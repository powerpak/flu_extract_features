#!/usr/bin/env ruby

# libraries provided by ruby
require 'rubygems'
require 'open-uri'

# gems you will need (run bundle install to get them)
require "bundler/setup"
require 'json'
require 'nokogiri'
require 'htmlentities'

class IRDClient
  
  def initialize(config=nil)
    @config = config
  end
  
end