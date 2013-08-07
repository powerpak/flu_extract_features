#!/usr/bin/env ruby

class IrdStrain < Model
  self.table_name = 'ird_strain'
  
  has_many :strain_name_ird_matches
  has_many :strain_names, :through => :strain_name_ird_matches
end