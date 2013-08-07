#!/usr/bin/env ruby

require 'bio'

class Bio::Sequence::NA
  # These describe a subset of HGVS-style variant specifications that we might want to apply
  # to a Bio::Sequence.  See http://www.hgvs.org/mutnomen/recs.html#general
  DNA_PATCH = /^([cg])\.(\d+)(_(\d+))?((del|ins)([ACTG]+)|([ACTG]+)>([ACTG]+))$/
  RNA_PATCH = /^([rn])\.(\d+)(_(\d+))?((del|ins)([acug]+)|([acug]+)>([acug]+))$/
  
  class << self
    def is_patch?(variant, type = nil)
      case type
      when :dna then variant =~ DNA_PATCH
      when :rna then variant =~ RNA_PATCH
      else variant =~ RNA_PATCH or variant =~ DNA_PATCH
      end
    end
  end
  
  def patch(variant, locs = nil)
    locs = Bio::Locations.new(locs) unless !locs || locs.is_a?(Bio::Locations)
    # The feature has to be complete (no range endpoints)
    raise "Inexact location cannot be patched" if locs && (locs.any?(&:gt) || locs.any?(&:lt))
    variant =~ DNA_PATCH or raise "Could not parse NA variant description"    # TODO: handle RNA?
    raise "DNA coding sequence variant '#{variant}' given without a gene location" if $1 != 'g' && !locs
    
    abs_origin = $1 == 'g' ? 0 : locs.absolute(1) - 1
    nil_seq = self.class.new('')
    from = $2.to_i + abs_origin
    to = $4 && $4.size > 0 ? $4.to_i + abs_origin : from
    orig = subseq(from, to)
    
    case $6
    when 'del'
      raise "DNA deletion variant says '#{variant}' but ref is: #{orig.upcase}" if orig != $7.downcase
      (from > 1 ? subseq(1, from - 1) : nil_seq) + (subseq(to + 1) || nil_seq)
    when 'ins'
      raise "DNA insertion variant '#{variant} should be 1 nt wide" if to != from + 1
      subseq(1, from) + self.class.new($7.downcase) + subseq(to)
    else
      raise "DNA substitution variant says '#{variant}' but ref is: #{orig.upcase}" if orig != $8.downcase
      replace = self.class.new $9.downcase
      (from > 1 ? subseq(1, from - 1) : nil_seq) + replace + (subseq(to + 1) || nil_seq)
    end
  end
  
end

class Bio::Sequence::AA
  AA_PATCH_SIMPLE = /^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$/
  # This describes a subset of HGVS-style variant specifications that we might want to apply
  # to a Bio::Sequence.  See http://www.hgvs.org/mutnomen/recs-prot.html#prot
  AA_PATCH = /^([p])\.(\d+)(_(\d+))?(del|ins)([ACDEFGHIKLMNPQRSTVWY]+)$/
  
  class << self
    def is_patch?(variant)
      variant =~ AA_PATCH or variant =~ AA_PATCH_SIMPLE
    end
    
    def calls_from_patch(variant)
      calls = {}
      if variant =~ AA_PATCH_SIMPLE
        calls[$2.to_i] = $3
      elsif variant =~ AA_PATCH
        from = $2.to_i
        to = $4 && $4.size > 0 ? $4.to_i : from
        case $5
        when 'del'
          (from..to).each {|i| calls[i] = '-' }
        when 'ins'
          calls[from] = '.' + $6.upcase
        end
      else
        raise "Could not parse protein-level variant description '#{variant}'"
      end
      calls
    end
  end
  
  def patch(variant)
    if variant =~ AA_PATCH_SIMPLE
      from = to = $2.to_i
      orig = subseq(from, to)
      raise "AA substitution variant says '#{variant}' but ref is: #{orig}" if orig != 'X' && orig != $1.upcase
      replace = self.class.new $3.upcase
      (from > 1 ? subseq(1, from - 1) : nil_seq) + replace + (subseq(to + 1) || nil_seq)
    elsif variant =~ AA_PATCH
      nil_seq = self.class.new('')
      from = $2.to_i
      to = $4 && $4.size > 0 ? $4.to_i : from
      orig = subseq(from, to)
      case $5
      when 'del'
        raise "AA deletion variant says '#{variant}' but ref is: #{orig}" if orig != $6.upcase
        (from > 1 ? subseq(1, from - 1) : nil_seq) + (subseq(to + 1) || nil_seq)
      when 'ins'
        raise "AA insertion variant '#{variant} should be 1 AA wide" if to != from + 1
        subseq(1, from) + self.class.new($6.upcase) + subseq(to)
      end
    else
      raise "Could not parse protein-level variant description '#{variant}'"
    end
  end
end
