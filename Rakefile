$:.unshift File.expand_path("../lib", __FILE__), File.expand_path("../ext", __FILE__)

# From gems or Ruby
require 'rubygems'
require 'yaml'
require 'shellwords'
require 'uri'
require 'andand'
require 'set'
require 'pp'

# Our libraries
require 'bio_ext'
require 'models'
require 'genome'

config = nil
REQUIRED_BINS = {"wget"=>"wget", "makeblastdb"=>"BLAST+"}

task :check do |t, args|
  if missing = REQUIRED_BINS.keys.find{|b| `which #{b}`.strip.size == 0 }
    fail "FAIL: Could not find \`#{missing}\` in your $PATH; please ensure #{REQUIRED_BINS[missing]} is installed."
  end
end

task :config => :check do |t, args|
  unless File.exist? "config.yaml"
    fail "Copy config.dist.yaml to config.yaml and edit before proceeding."
  end
  config = YAML.load_file("config.yaml")
  config["segments_by_gene"] = {}
  config["prefixes_by_gene"] = {}
  config["genes_by_segment"].each_with_index do |genes, i|
    genes.each do |gene, prefixes| 
      config["prefixes_by_gene"][gene] = [*prefixes]
      config["segments_by_gene"][gene] = i + 1
    end
  end
  config["pwd"] = pwd
  database_cfg = config["database"].each_with_object({}){|(k,v), h| h[k.to_sym] = v }
  ActiveRecord::Base.establish_connection(database_cfg)
  Model.config = config
end


desc "Start a console with the database models loaded."
task :console => :config do
  require 'irb'
  ARGV.clear
  IRB.start
end


desc "Loads db/ird_strains.tsv into the database (see db/README.md for details)"
task :load_ird_strains => :config do |t, args|
  Rake::Task["db/ird_strains.tsv"].invoke
  KEYS = config["ird_strains_cols"].map &:to_sym
  File.open("db/ird_strains.tsv") do |f|
    first_line = true
    strains = []
    f.each_line do |tsv|
      if first_line then first_line = false; next; end
      fields = tsv.chomp.split(/\t/).map{|v| v =='-N/A-' ? '' : v }
      strains << IrdStrain.new(Hash[KEYS.zip(fields)])
      if strains.length >= 1000
        puts "Importing #{strains.size} rows, now at: #{fields[0]}"
        IrdStrain.import strains
        strains = []
      end
    end
    if strains.length > 0
      puts "Importing final #{strains.size} rows."
      IrdStrain.import strains
    end
    puts "Done."
  end
end


desc "[1] Updates the mapping of strain names to IRD strains"
task :update_strain_names => :config do |t, args|
  MOD_COLS = (1..config["number_mod_cols"]).map{|n| "mod_#{n}"}
  COLS = "id, strain_name, #{MOD_COLS.join ', '}"
  
  # Clear all the existing strain names and strain name to phenotype relationships
  StrainNamePhenotype.delete_all and StrainNamePhenotype.reset_auto_increment
  StrainName.delete_all and StrainName.reset_auto_increment
    
  # Scan phenotypes for strain names, create entries for them, and link them to the phenotypes
  Phenotype.select(COLS).find_each do |pheno|
    sn = StrainName.where(:strain_name => pheno.strain_name).first_or_create
    sn.phenotypes << pheno unless sn.phenotypes.include? pheno
    puts "#{pheno.strain_name} linked to phenotype #{pheno.id}."
    MOD_COLS.each do |col|
      if pheno.send(col.to_sym) =~ /^\w+:(A\/.+)$/
        sn = StrainName.where(:strain_name => $1).first_or_create
        sn.phenotypes << pheno unless sn.phenotypes.include? pheno
        puts "#{$1} linked to phenotype #{pheno.id}."
      end
    end
  end
  
  # Link all strain names to their likely corresponding IRD strains based on synonyms
  StrainName.link_to_ird_strains! true
end


desc "[2] Caches and pre-processes all GenBank entries relevant to phenotypes"
task :update_genbank_entries => :config do |t, args|
  StrainName.select("id, strain_name").find_each do |sn|
    sn.ird_strains.find_each do |ird_strain|
      config["ird_strains_cols"].drop(1).each do |col|
        ird_strain.send(col.to_sym).gsub(/[*\s]/, '').split(',').each do |acc|
          puts "Caching #{acc}..."
          GenbankEntry.cache(acc)
        end
      end
    end
  end
end


directory "db/PROTEIN-A"
file "db/PROTEIN-A" => :config do
  uri = URI.parse config["endpoints"]["annotation"]
  cut = uri.path.split('/').size
  here = pwd
  cd 'db/PROTEIN-A'
  sh "wget -r -nH --cut-dirs=#{cut} #{Shellwords.escape uri.to_s}"
  cd here
end

file "db/protein-a.fasta" => "db/PROTEIN-A" do
  sh "find 'db/PROTEIN-A' -type f | xargs cat > db/protein-a.fasta"
end

directory "db/blastdb"
file "db/blastdb" => "db/protein-a.fasta" do
  sh "makeblastdb -dbtype prot -in db/protein-a.fasta -title influenza-a -out db/blastdb/influenza-a"
end


desc "[3] Compiles GenBank entries, applies alterations, and links phenotypes to sequence features"
task :update_sequence_features, [:reset] => [:config, "db/protein-a.fasta", "db/blastdb"] do |t, args|
  MOD_COLS = (1..config["number_mod_cols"]).map{|n| "mod_#{n}"}
  COLS = "id, strain_name, #{MOD_COLS.join ', '}"
  variant_caller = VariantCaller.new(config)
  
  # We could clear all existing genotypes, features, and their relationships to be neat...
  if !args.reset.nil?
    Feature.delete_all and Feature.reset_auto_increment
    GenotypeFeature.delete_all and GenotypeFeature.reset_auto_increment
    Genotype.delete_all and Genotype.reset_auto_increment
    GenotypePhenotype.delete_all and GenotypePhenotype.reset_auto_increment
  end
  
  complete = 0
  at_least_one = 0
  aa_successes = 0; aa_attempts = 0; pheno_killed = Set.new
  Phenotype.select(COLS).find_each do |pheno|
    puts "Working on phenotype #{pheno.id}."
    
    # Get the base strain's Genbank entries according to IRD, arranged as an array of 8 segments
    base_strain = StrainName.find_by(:strain_name => pheno.strain_name)
    # TODO: refactor base_segs into a GenomeSegments class
    base_segs = base_strain.ird_genbank_entries 
    
    # Gather all the alterations into a list of e.g. {:gene => 'PB1', :alt => 'U123V'}
    # TODO: refactor into AlterationList.new; AlterationList#push
    alterations = []
    MOD_COLS.map(&:to_sym).each do |col| 
      mod = pheno.send(col)
      if mod && mod =~ /^(\w+):(.*)$/
        $2.split(',').each {|alt| alterations << {:gene => $1, :alt => alt} }
      end
    end
    
    # Apply reassortment alterations, removing them from the list as we go
    # TODO: refactor into AlterationList#apply_reassortments base_segs
    alterations.reject! do |alt|
      if alt[:alt] =~ /^[ABC]\/.+$/
        strain_name = StrainName.find_by(:strain_name => alt[:alt])
        seg_num = config["segments_by_gene"][alt[:gene]]
        base_segs[seg_num - 1] = strain_name.ird_genbank_entries(seg_num)
        true
      end
    end
    # Track completeness: whether we have enough complete Genbank entries to work with
    completeness = base_segs.map {|entries| entries.any? &:complete? }
    complete += 1 if completeness.all?
    at_least_one += 1 if completeness.any?
    
    # Apply nucleotide-level alterations, removing them from the list as we go
    # TODO: refactor into AlterationList#apply_nt_alterations base_segs
    alterations.reject! do |alt|
      # Only allow coding sequence changes: genomic coordinates are too unreliable for influenza
      if alt[:alt] =~ /^c\./ && Bio::Sequence::NA.is_patch?(alt[:alt], :dna)
        seg_num = config["segments_by_gene"][alt[:gene]]
        base_segs[seg_num - 1].map! do |entry|
          ivsat_gb = Bio::GenBank.new(entry.ivsat_genbank)
          gene = ivsat_gb.features.find{|f| f.assoc['gene'] == alt[:gene] }
          next unless gene  # The GenBank file has to include the gene we're altering, otherwise throw it out
          
          # #patch is defined in bio_ext: applies the alteration to the naseq
          # This is wrapped immediately into an object that re-annotates the sequence with IVSAT
          IVSATAnnotatableSequence.new(ivsat_gb.naseq.patch(alt[:alt], gene.position), config)
        end.compact!
        true
      end
    end
    
    # At this point we want to collect all predicted protein sequences by gene name
    # TODO: refactor into GenomeSegments#translate --> ProteinSequences
    prot_seqs = Hash.new
    base_segs.each_with_index do |seg, i|
      seg.each do |entry|
        ivsat_gb = Bio::GenBank.new(entry.ivsat_genbank)
        gene = ivsat_gb.features.map(&:assoc).each do |f| 
          if f['translation'] && f['gene']
            gene = config['gene_synonyms'][f['gene']] || f['gene']
            (prot_seqs[gene] ||= []) << IVSATProteinSequence.new(f, entry)
          end
        end
      end
    end
    # If we have a complete sequence for a gene, throw out the partials
    # Why?  Otherwise, the logic for the following step becomes much more complicated.
    prot_seqs.each do |gene, seqs|
      seqs.reject! {|seq| !seq.complete? } if seqs.any?(&:complete?)
    end
    
    # Apply AA-level alterations for protein sequences that are complete, removing them from the list as we go
    # TODO: refactor into AlterationList#apply_aa_alterations prot_seqs --> (boolean) killed
    alterations.reject! do |alt|
      if Bio::Sequence::AA.is_patch?(alt[:alt])
        if prot_seqs[alt[:gene]].andand.first.andand.complete?
          prot_seqs[alt[:gene]].reject! do |seq|
            aa_attempts += 1
            begin
              seq.patch! alt[:alt]
              aa_successes += 1
              false
            rescue Exception => e
              puts "#{alt[:gene]}, #{seq.acc}: #{e.message}"
              true
            end
          end
          if prot_seqs[alt[:gene]].size == 0
            puts "Killed by #{alt[:gene]}:#{alt[:alt]}"; pheno_killed << pheno.id
          end
          true
        end
      end
    end
    
    # TODO: refactor into ProteinSequences#call_variants --> VariantList
    #       (right now it's just a hash of genes -> lists of VariantCalls)
    # If we weren't killed by problems with applying AA-level alterations, it's time for blastp
    next if pheno_killed.include? pheno.id
    var_list = {}
    prot_seqs.each do |gene, seqs|
      seg_num = config["segments_by_gene"][gene]
      prefixes = config["prefixes_by_gene"][gene]
      var_list[gene] = seqs.map do |seq|
        variant_caller.blastp(seq.aa_seq, gene) do |hit|
          refname = hit.target_def
          refname =~ /^[Ss]eg#{seg_num}\D/ && !refname.include?('mature') && prefixes.any?{|p| refname =~ /#{p}(\D|$)/ } 
        end
      end
    end
    
    # Apply remaining AA-level and knockout alterations, removing them from the list as we go
    alterations.reject! do |alt|
      if Bio::Sequence::AA.is_patch?(alt[:alt]) || alt[:alt] == 'KO'
        puts "A straggler alteration appeared! #{alt[:gene]}:#{alt[:alt]}"
        simulated_calls = VariantCalls.new(alt[:gene], alt[:alt])
        puts "It could not be applied." if !var_list[alt[:gene]]
        # TODO: Should we remove variants on genes for which we have no GenBank sequences from the list, 
        # or should they kill the genotype-in-progress?  Hmmm, not sure
        var_list[alt[:gene]].andand.map! do |var_calls|
          var_calls + simulated_calls
        end
        true
      end
    end
    
    # Are there any alterations left in the list?  Hopefully not.
    if alterations.size > 0
      alterations.each {|alt| puts "Alteration could not be applied! #{alt[:gene]}:#{alt[:alt]}" }
      pheno_killed << pheno.id
      puts "Killed."
      next
    end
    
    # Save all the VariantCalls to the feature and genotype tables.
    # WARNING: This deletes and overwrites the existing genotype for the phenotype!
    # As they are saved, VariantCalls should be checked against their reference, and redundant calls removed
    #   (e.g., if a call is made that is identical to the reference don't save it.)
    # TODO: Refactor into VariantList#save_as_genotype!(pheno)
    pheno.genotypes.each {|geno| geno.delete }
    pheno.reload
    pheno.genotypes = [Genotype.from_variant_list(var_list)]
  end
  puts "#{complete} complete genomes and #{at_least_one} with >1 complete segment (out of #{Phenotype.count})."
  puts "#{aa_attempts} AA substitutions attempted, #{aa_successes} succeeded."
  puts "#{pheno_killed.size} phenotypes were killed by failed AA substitutions."
end

desc "[4] Extracts sequence features and phenotypes from the database into an ARFF file"
task :make_arff, [:file] => :config do |t, args|
  $stdout = File.open(args.file || "influenza-a.arff", 'w') unless args.file == '-'
  
  PHENOTYPE_FILTER = "pathogenicity_qual IN ('low', 'high')"
  PHENOTYPE_ATTRIBUTES = [
    ["pathogenicity_qual", "{low,high}"]
  ]
  COLS = "id," + PHENOTYPE_ATTRIBUTES.map(&:first).join(',')
  
  # Print the header @RELATION line
  puts "@RELATION influenza-a-virulence\n\n"
    
  # First, we need to tabulate all possible features, which become separate @ATTRIBUTE's in the ARFF.
  # Each attribute is nominal and all the possible values also need to be tabulated.
  # We should collect the variant_call features by gene, as well.
  attr_tree = {:knockout=>[], :reference=>[], :variant_call=>{}}
  Feature.select("id, kind, name").find_each do |feat|
    possible_values = feat.genotype_features.pluck(:value).uniq
    if feat.kind == 'variant_call'
      gene, pos = feat.name.split(':')
      (attr_tree[:variant_call][gene] ||= []) << [feat.name, possible_values, pos.to_i]
    else
      attr_tree[feat.kind.to_sym] << [feat.name, possible_values]
    end
  end
  
  # Sort the attributes in a sensible manner and print the @ATTRIBUTE lines
  attr_cols = []
  attr_tree.each do |kind, sub_attrs|
    if kind == :variant_call
      # For each gene in the variant_call subtree, sort the attributes by their numerical position (prettier output)
      sub_attrs.each {|gene, attrs| sub_attrs[gene] = attrs.sort_by {|attr| attr[2] }}
      sub_attrs.keys.sort.each do |gene|
        sub_attrs[gene].each do |attr|
          attr_name = "variant_call:#{attr[0]}"
          puts "@ATTRIBUTE \"#{attr_name}\" {\"#{(attr[1] <<= '.').join('","')}\"}"
          attr_cols << attr_name
          attr[3] = attr_cols.size - 1
        end
      end
    else
      attr_tree[kind] = sub_attrs.sort_by {|attr| attr[0] }
      sub_attrs.each do |attr|
        attr_name = "#{kind}:#{attr[0]}"
        puts "@ATTRIBUTE \"#{attr_name}\" {#{attr[1].join(',')}}"
        attr_cols << attr_name
        attr[3] = attr_cols.size - 1
      end
    end
  end
  # Print @ATTRIBUTE lines for the phenotype columns
  PHENOTYPE_ATTRIBUTES.each do |pa|
    puts "@ATTRIBUTE \"#{pa[0]}\" #{pa[1]}"
    attr_cols << pa[0]
    pa[2] = attr_cols.size - 1
  end
  $stderr.puts "Wrote attributes."
  
  # And now for the data...
  puts "\n@DATA\n"
  # Iterate over phenotypes, branching on related genotypes (although there should be one per at this point)
  # Every phenotype-genotype will be one row of data for the ARFF file.
  Phenotype.select(COLS).where(PHENOTYPE_FILTER).find_each do |pheno|
    # Although we're currently using each genotype only once, we could conceivably use the weights on genotype_features
    # to create multiple observations for combinations of conflicting variant calls.
    # For now, we just take the most heavily weighted observation at each position.
    pheno.genotypes.each do |g|
      # If we have a reference feature for a gene, then we should include all available variant call features for that gene.
      # Variant call features that are not in the DB are set to ".", indicating it is the same as reference sequence.
      # For the genes with no reference feature, that and all variant call features are set to "?" (the unknown value)
      # Don't forget the knockout feature.
      data_vals = Array.new(attr_cols.size)
      
      attr_tree[:reference].each do |attr|
        gf = g.genotype_features.joins(:feature).where('feature.kind' => 'reference', 'feature.name'=>attr[0]).first
        next unless gf
        
        # We have a reference feature for this gene
        data_vals[attr[3]] = gf.value
        filter = "feature.kind = 'variant_call' AND feature.name LIKE '#{attr[0].gsub(/[^a-zA-Z0-9-]/, '')}:%'"
        var_calls = {}
        g.genotype_features.joins(:feature).where(filter).each do |gf|
          # This is how the most heavily weighted observation at each position is picked
          if !var_calls[gf.feature.name] or var_calls[gf.feature.name][1] < gf.weight
            var_calls[gf.feature.name] = [gf.value, gf.weight]
          end
        end
        attr_tree[:variant_call][attr[0]].each do |vc_attr|
          data_vals[vc_attr[3]] = var_calls[vc_attr[0]] && var_calls[vc_attr[0]][1] >= 0.5 ? var_calls[vc_attr[0]][0] : '.'
        end
      end
      data_vals.map!{|v| v.nil? ? '?' : v }
      
      # Finally add the phenotype data values
      PHENOTYPE_ATTRIBUTES.each do |pa|
        pheno_val = pheno.send(pa[0].to_sym)
        data_vals[pa[2]] = pheno_val && pheno_val != '' ? pheno_val : '?'
      end
      puts data_vals.join ','
      $stderr.puts "Finished pheno #{pheno.id}."
    end
  end
  $stdout.flush
end