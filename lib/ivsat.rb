#!/usr/bin/env ruby

# libraries provided by ruby
require 'rubygems'
require 'uri'
require 'net/http'
require 'yaml'

# gems you will need (run bundle install to get them)
require 'bundler/setup'
require 'json'
require 'nokogiri'
require 'htmlentities'

# libraries we define
require 'models'

class IVSATClient
  
  def initialize(config=nil)
    @config = config
  end
  
  def cached(sequence)
    IvsatCache.where(:sequence => sequence).first_or_create
  end
  
  def annotate(sequence)
    if (row = cached(sequence)) && row.output
      YAML::load(row.output)
    else
      uri = URI.parse(@config['endpoints']['ivsat'])
      http = Net::HTTP.new(uri.host, uri.port)
      request = Net::HTTP::Post.new(uri.path)
      request.set_form_data(:sequence => sequence, :SUBMIT => 'Annotate FASTA')
      doc = Nokogiri::HTML(http.request(request).body)
      output = {
        :genbank => doc.xpath("//input[@name='gbf_h']").first['value'],
        :asn1 => doc.xpath("//input[@name='sqn_h']").first['value'],
        :feature_table => doc.xpath("//input[@name='tbl_h']").first['value'],
        :xml => doc.xpath("//input[@name='xml_h']").first['value'],
        :fasta_protein => doc.xpath("//input[@name='faa_h']").first['value'],
        :fasta_na => doc.xpath("//input[@name='ffn_h']").first['value'],
        :alignment => doc.xpath("//input[@name='aln_h']").first['value']
      }
      row.output = YAML::dump(output)
      row.save
      output
    end
  end
  
end

module IVSATAnnotatable
  
  def ivsat_annotate(ivsat_client = nil)
    ivsat_client ||= @ivsat_client || @@ivsat_client
    unless ivsat_table && ivsat_genbank && ivsat_alignment
      annotations = ivsat_client.annotate(na_seq)
      self.ivsat_table = annotations[:feature_table]
      self.ivsat_genbank = annotations[:genbank]
      self.ivsat_alignment = annotations[:alignment]
      puts "Fetched new IVSAT annotations#{self.respond_to?(:acc) ? " for #{acc}" : ''}."
    end
  end
  
  # Public: Returns the completeness of the genes within this sequence according 
  # to the analysis in ivsat_table
  def complete?
    return nil unless ivsat_table
    line = ivsat_table.split("\n").find{|l| l.include? 'INFO: Sequence completeness:' }
    line && !line.include?('partial')
  end
  
end

class IVSATAnnotatableSequence
  
  attr_reader :na_seq 
  attr_accessor :ivsat_table, :ivsat_genbank, :ivsat_alignment
  include IVSATAnnotatable

  def initialize(na_seq, config=nil)
    @na_seq = na_seq
    @config = config
    @ivsat_client = IVSATClient.new config
    
    # Populates ivsat_table, ivsat_genbank, and ivsat_alignment
    ivsat_annotate
  end
    
end

class IVSATProteinSequence
  
  attr_accessor :aa_seq
  attr_reader :completeness, :gene, :acc
  
  def initialize(feature_assoc, entry)
    @aa_seq = Bio::Sequence::AA.new feature_assoc['translation']
    @completeness = entry.complete?
    @gene = feature_assoc['gene']
    @acc = entry.respond_to?(:acc) && entry.acc
  end
  
  def patch!(variant); @aa_seq = @aa_seq.patch variant; end
  
  def complete?; @completeness; end
  
end
