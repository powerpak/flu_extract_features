#!/usr/bin/env ruby

class Feature < Model
  self.table_name = 'feature'
  
  has_many :genotype_features
  has_many :genotypes, :through => :genotype_features
end