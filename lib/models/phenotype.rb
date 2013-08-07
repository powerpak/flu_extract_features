#!/usr/bin/env ruby

class Phenotype < Model
  self.table_name = 'phenotype'
  
  belongs_to :source
  
  has_many :strain_name_phenotypes
  has_many :strain_names, :through => :strain_name_phenotypes
  
  has_many :genotype_phenotypes
  has_many :genotypes, :through => :genotype_phenotypes
end