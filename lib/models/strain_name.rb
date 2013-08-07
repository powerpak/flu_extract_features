#!/usr/bin/env ruby

require_relative './ird_strain'
require_relative './strain_name_ird_match'
require_relative './genbank_entry'
require 'set'

class StrainName < Model
  self.table_name = 'strain_name'
  
  has_many :strain_name_ird_matches
  has_many :ird_strains, :through => :strain_name_ird_matches
  
  has_many :strain_name_phenotypes
  has_many :phenotypes, :through => :strain_name_phenotypes
  
  @@synonym_dict = nil
  
  class << self
    def link_to_ird_strains!(noisy = false)
      error_phenotypes = Set.new
      StrainNameIrdMatch.delete_all and StrainNameIrdMatch.reset_auto_increment
      self.select("id, strain_name").find_each do |sn|
        # Link this StrainName to all possible combinations of part_synonyms
        if !sn.link_to_ird_strains!(noisy)
          ids = sn.phenotypes.pluck(:id).join ','
          puts "#{sn.strain_name} (phenotypes #{ids}) could not be linked to an IRD strain." if noisy
          error_phenotypes.merge sn.phenotypes
        end
      end
      puts "#{error_phenotypes.size} phenotypes (out of #{Phenotype.count}) have strain names that could not be linked to IRD." if noisy
    end
  end
  
  # Internal: Loads the synonym sets from the config file and pre-processes it into a dictionary.
  # This is only done once, after which the result is cached in a class variable.
  #
  # Returns the synonym dictionary.
  def synonym_dict
    return @@synonym_dict if @@synonym_dict
    @@synonym_dict = {}
    config["synonyms"].each do |syn_set|
      syn_set.each {|syn| @@synonym_dict[syn] ||= Set.new; @@synonym_dict[syn].merge syn_set }
    end
    @@synonym_dict
  end
  
  # Public: Links this StrainName to possible IRDStrains using a list of names to try.
  #
  # Returns the linked StrainName.
  def link_to_ird_strains!(noisy = false)
    cutoff_year = Time.now.year % 100
    
    parts = strain_name.split(/\//)
    part_synonyms = parts.each_with_index.map do |part, i|
      part.gsub! /^\s+|\s+$/, ''
      synonyms = case
      when i == 0
        # First part, the type of the influenza virus, we'll leave alone
        [part]
      when (i == parts.size - 1 and part =~ /^(\d{2}?)(\d{2})$/)
        # Try two digit AND four digit date, assuming that two digit dates past this year's are in the 1900's
        if $1.size > 0 then [part, $2] else [part, "#{$2.to_i > cutoff_year ? 19 : 20}#{$2}"] end
      when part =~ /^0*(\d+)$/
        # Allow some flexibility in leading zeros for middle segments containing numbers
        [$1, "0#{$1}", "00#{$1}", "000#{$1}"]
      when synonym_dict[part]
        synonym_dict[part].to_a
      else
        [part]
      end
      # Now add more possibilities: hyphens removed, underscores vs. spaces
      synonyms.map do |syn|
        more_syns = [syn]
        more_syns << syn.gsub(/-/, '') if syn =~ /-/
        more_syns << syn.gsub(/ /, '_') if syn =~ /_/
        more_syns << syn.gsub(/_/, ' ') if syn =~ / /
        more_syns
      end.flatten
    end
    
    part_synonyms.first.product(*part_synonyms.drop(1)).each do |parts|
      name = parts.join('/').gsub('_', '\\_').gsub('%', '\\%')
      IrdStrain.where("strain_name LIKE ?", name).find_each do |ird_strain|
        puts "#{strain_name} linked to IRD strain #{ird_strain.id}, name #{ird_strain.strain_name}." if noisy
        ird_strains << ird_strain unless ird_strains.include? ird_strain
      end
    end
    
    ird_strains.size > 0
  end
  
  # Public: Gets all Genbank entries according to IRD for this StrainName.
  #
  # seg_num - optionally, provide the number of the segment that you are interested in.
  #           This number is 1-indexed (the first segment is 1).
  #
  # Returns an array of eight arrays of GenbankEntries for every segment, or if a seg_num was
  # provided, just the array of GenbankEntries for that segment.
  # TODO: refactor this into returning a GenomeSegments class
  def ird_genbank_entries(seg_num = nil)
    seg_cols = seg_num ? [config["ird_strains_cols"][seg_num]] : config["ird_strains_cols"].drop(1)
    segments = seg_cols.map{ [] }
    ird_strains.each do |ird_strain|
      seg_cols.each_with_index do |col, i|
        ird_strain.send(col.to_sym).gsub(/[*\s]/, '').split(',').each do |acc|
          segments[i] << GenbankEntry.for_acc(acc)
        end
      end
    end
    seg_num ? segments[0] : segments
  end
  
end