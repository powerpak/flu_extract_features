#!/usr/bin/env ruby

require 'bio'
require 'ivsat'

class GenbankEntry < Model
  self.table_name = 'genbank_entry'
  
  include IVSATAnnotatable
  
  class << self
    def cache(acc)
      row = self.where(:acc => acc).first_or_create
      row.genbank ||= download(acc)
      row.na_seq ||= Bio::GenBank.new(row.genbank).naseq
      row.ivsat_annotate ivsat_client
      if row.changed? then row.save && row else row end
    end
    
    alias_method :for_acc, :cache
    
    def download(acc)
      open(config["endpoints"]["genbank"] % [acc]).read
    end
    
    def ivsat_client
      @@ivsat_client ||= IVSATClient.new config
    end
  end
  
end