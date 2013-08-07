class GenotypeFeature < Model
  self.table_name = 'genotype_feature'
  
  belongs_to :genotype
  belongs_to :feature
  
  class << self
    def create_for_feature(kind, name, value, weight=1.0)
      row = self.new
      row.feature = Feature.where(:kind => kind, :name => name).first_or_create
      row.value = value
      row.weight = weight
      row
    end
  end
end