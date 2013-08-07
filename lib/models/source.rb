#!/usr/bin/env ruby

class Source < Model
  self.table_name = 'source'
  
  has_many :phenotypes
end