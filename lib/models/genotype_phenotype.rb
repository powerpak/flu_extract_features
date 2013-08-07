#!/usr/bin/env ruby

class GenotypePhenotype < Model
  self.table_name = 'genotype_phenotype'
  
  belongs_to :genotype
  belongs_to :phenotype
end