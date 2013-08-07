#!/usr/bin/env ruby

require 'bio'


class VariantCaller
  
  def initialize(config)
    @config = config
    @factory = Bio::Blast.local('blastp', 'db/blastdb/influenza-a')
    @factory.blastall = "#{config['pwd']}/bin/blastall"
  end
  
  def cached(sequence, gene_name)
    VariantCallerCache.where(:sequence => sequence, :gene_name => gene_name).first_or_create
  end
  
  def blastp(sequence, gene_name)
    if (row = cached(sequence, gene_name)) && row.output
      YAML::load(row.output)
    else
      hits = @factory.query(sequence).iterations.last.hits
      top_hit = if block_given? then hits.find &Proc.new else hits.first end
      unless top_hit
        raise "No valid blastp matches for #{gene_name} (Other hits: #{hits.map(&:target_def).join(', ')})!"
      end
      calls = call_blastp_hit(top_hit)
      output = VariantCalls.new(gene_name, calls, top_hit.target_def)
      row.output = YAML::dump(output)
      row.save
      output
    end
  end
  
  def call_blastp_hit(hit)
    calls = {}
    pos = hit.target_start  # pos is a counter for where we are in the target sequence.
    qseq = hit.query_seq
    # I do not believe we need to consider top_hit.midline (conservative vs. nonconservative substitutions)
    hit.target_seq.split("").each_with_index do |char, i|
      case char
      when "-"  # Insertion
        # Append the query sequence's character to the call for the previous position
        calls[pos - 1] = (calls[pos - 1] || '.') + qseq[i] unless pos <= 1
        # Do NOT advance pos for an insertion
      when "X" # Ambigious reference sequence
        # In this case, we only care about deletions, as an AA difference here is not significant
        calls[pos] = "-" if qseq[i] == "-"
        pos += 1
      else # Substitution or deletion
        # Call the position if the query sequence is different and nonambiquous
        calls[pos] = qseq[i] unless char == qseq[i] || qseq[i] == "X"
        pos += 1
      end
    end
    calls
  end
  
end


class VariantCalls
  
  attr_reader :gene, :reference, :calls, :knockout
  
  include Enumerable
  extend Forwardable
  
  def initialize(gene, calls={}, reference=nil, knockout=false)
    @gene = gene
    @reference = reference
    @knockout = knockout
    @calls = {}
    if calls.is_a?(String)
      if calls == 'KO' then @knockout = true
      else @calls = Bio::Sequence::AA.calls_from_patch(calls); end
    elsif calls.is_a?(Hash)
      @calls = calls
    end
  end
  
  # Variant calls can be added.  Calls in the second operand overwrite calls in the first at the same position.
  def +(vcs)
    # The gene must be identical, as well as the reference, although a nil reference for one operand is OK.
    raise "Cannot add variant calls from different genes" if vcs.gene && @gene && vcs.gene != @gene
    raise "Cannot add variant calls from different references" if vcs.reference && @reference && vcs.reference != @reference
    
    # If gene has been knocked out, our work is done
    return self.class.new(@gene || vcs.gene, {}, @reference || vcs.reference, true) if vcs.knockout || @knockout
    
    # Clone so as to not modify the existing object
    new_calls = @calls.clone
    vcs.each do |pos, change|
      # Only one exception to the relatively simple process of overriding calls in the new VariantCalls object:
      # If there is an *addition* variant with the initial position unspecified (a period), any existing variant
      # for that one position is preserved.  This allows an addition variant to combine with a substitution variant
      # without losing information.
      prev_call = new_calls[pos]
      if change.size > 1 && change[0] == '.' && prev_call
        new_calls[pos] = prev_call[0] + change[1..-1]
      elsif prev_call && prev_call.size > 1 && prev_call[0] == '.' && change.size == 1
        new_calls[pos] = change + prev_call[1..-1]
      else
        new_calls[pos] = change
      end
    end
    self.class.new(@gene || vcs.gene, new_calls, @reference || vcs.reference)
  end
  
  # Allow iteration and indexing into the @calls attribute directly from this object
  def_delegators :@calls, :size, :[], :each, :count
  def_delegator :@calls, :keys, :positions
  
end