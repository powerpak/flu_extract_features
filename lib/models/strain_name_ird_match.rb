class StrainNameIrdMatch < Model
  self.table_name = 'strain_name_ird_match'
  
  belongs_to :strain_name
  belongs_to :ird_strain
end