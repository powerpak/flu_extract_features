#!/usr/bin/env ruby

require 'models'

Dir[File.expand_path("../genome/*.rb", __FILE__)].each {|file| require file.sub(/.rb$/, '') }