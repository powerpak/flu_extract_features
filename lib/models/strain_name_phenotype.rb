class StrainNamePhenotype < Model
  self.table_name = 'strain_name_phenotype'
  
  belongs_to :strain_name
  belongs_to :phenotype
end