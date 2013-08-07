#!/usr/bin/env ruby

require 'mysql2'
require 'active_record'
require 'activerecord-import/base'
ActiveRecord::Import.require_adapter('mysql2')

class Model < ActiveRecord::Base
  self.abstract_class = true
  class_attribute :config
  
  class << self
    def reset_auto_increment
      ActiveRecord::Base.connection.execute "ALTER TABLE `#{self.table_name}` AUTO_INCREMENT = 1;"
    end
  end
end

Dir[File.expand_path("../models/*.rb", __FILE__)].each {|file| require file.sub(/.rb$/, '') }