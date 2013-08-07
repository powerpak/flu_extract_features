#!/usr/bin/env ruby

require_relative './feature'
require 'bio'

class Genotype < Model
  self.table_name = 'genotype'
  @@references = nil
  
  has_many :genotype_phenotypes
  has_many :phenotypes, :through => :genotype_phenotypes
  
  has_many :genotype_features
  has_many :features, :through => :genotype_features
  
  class << self
    # Retrieves the reference sequences we used to generate variant calls.
    def reference_protein_seqs
      if !@@references
        @@references = {}
        ff = Bio::FlatFile.new(Bio::FastaFormat, open("#{config['pwd']}/db/protein-a.fasta"))
        ff.each {|entry| @@references[entry.definition] = entry.seq }
      end
      @@references
    end
    
    def from_variant_list(var_list)
      row = self.new
      var_list.each do |gene, gene_var_calls|
        
        # We have to correct for the possibility that different var_calls will use different references.
        # In this case, we use most popular reference and discard the others.
        # TODO: perhaps, rebase the var_calls built on the unpopular reference onto the popular reference
        most_popular_ref = gene_var_calls.map(&:reference).group_by{|i| i}.map{|k, v| [k, v.count] }.max_by{|a| a[1] }
        gene_var_calls.reject!{|var_calls| var_calls.reference != most_popular_ref[0] }
        
        # The most popular reference is a GenotypeFeature.  All the variant_call features for that gene are built on it.
        gf_reference = GenotypeFeature.create_for_feature('reference', gene, most_popular_ref[0])
        row.genotype_features << gf_reference
        
        if gene_var_calls.any?(&:knockout)
          # If the variant for this gene was a knockout, save that as a feature
          raise "All variant calls should be a knockout for #{gene} if there was one!" unless gene_var_calls.all?(&:knockout)
          gf_reference = GenotypeFeature.create_for_feature('knockout', gene, "true")
          row.genotype_features << gf_reference
        else
          # This not a knocked out gene, save individual variant calls.
          # If there were multiple sequences for a segment in IRD, there may be conflicting variant calls.
          # We resolve this by adding a weight, which is the percentage each call was observed.
          reference_seq = Bio::Sequence::AA.new(reference_protein_seqs[most_popular_ref[0]])
          num_reads = gene_var_calls.size.to_f
          all_positions = gene_var_calls.map(&:positions).flatten.uniq
          var_calls_weighted = Hash[all_positions.map do |pos|
            # Also, we check all calls against the reference, and drop any that are identity variants.
            calls_at_pos = gene_var_calls.map{|vc| vc[pos] != reference_seq.subseq(pos, pos) ? vc[pos] : nil }.compact
            [pos, Hash[calls_at_pos.group_by{|i| i}.map{|k, v| [k, v.count / num_reads] }]]
          end]
          var_calls_weighted.each do |pos, weighted_calls|
            weighted_calls.each do |change, weight|
              gf_reference = GenotypeFeature.create_for_feature('variant_call', "#{gene}:#{pos}", change, weight)
              row.genotype_features << gf_reference
            end
          end
        end
        
      end
      row.save && row
    end
  end
  
end