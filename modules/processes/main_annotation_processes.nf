#!/usr/bin/env nextflow


/*
 * Process 1: Get all the parameters
 */

process getParams {

    output:
    val params

    exec:
    if (params.species == 'human' && params.genome_version == 'hg38') {
        params = [
            CpG_species: params.CpG_species ?: 'hg38',
            TSS_species: params.TSS_species ?: "Homo' 'sapiens",
            vep_species: params.vep_species ?: 'homo_sapiens',
            vep_cache_species: params.vep_cache_species ?: 'GRCh38',
			fasta_start: params.fasta_start ?: 'Homo_sapiens'
        ]
    } else if (params.species == 'human' && params.genome_version == 'hg19') {
        params = [
            CpG_species: params.CpG_species ?: 'hg19',
            TSS_species: params.TSS_species ?: "Homo' 'sapiens",
            vep_species: params.vep_species ?: 'homo_sapiens',
            vep_cache_species: params.vep_cache_species ?: 'GRCh37',
			fasta_start: params.fasta_start ?: 'Homo_sapiens'
        ]
    } else if (params.species == 'cattle' && params.genome_version == 'bosTau9') {
        params = [
            CpG_species: params.CpG_species ?: 'bosTau9',
            TSS_species: params.TSS_species ?: "Bos' 'taurus",
            vep_species: params.vep_species ?: 'bos_taurus',
            vep_cache_species: params.vep_cache_species ?: 'ARS-UCD1.2',
			fasta_start: params.fasta_start ?: 'Bos_taurus'
        ]
    } else if (params.species == 'pig' && params.genome_version == 'Sscrofa11.1') {
        params = [
            CpG_species: params.CpG_species ?: 'Sscrofa11.1',
            TSS_species: params.TSS_species ?: "Sus' 'scrofa",
            vep_species: params.vep_species ?: 'sus_scrofa',
            vep_cache_species: params.vep_cache_species ?: 'Sscrofa11.1',
			fasta_start: params.fasta_start ?: 'Sus_scrofa'
        ]
    } else if (params.species == 'dog' && params.genome_version == 'CanFam3.1') {
        params = [
            CpG_species: params.CpG_species ?: 'CanFam3.1',
            TSS_species: params.TSS_species ?: "Canis' 'lupus' 'familiaris",
            vep_species: params.vep_species ?: 'canis_lupus_familiaris',
            vep_cache_species: params.vep_cache_species ?: 'CanFam3.1',
			fasta_start: params.fasta_start ?: 'Canis_lupus_familiaris'
        ]
    } else if (params.species == 'sheep' && params.genome_version == 'Oar_v4.0') {
        params = [
            CpG_species: params.CpG_species ?: 'Oar_v4.0',
            TSS_species: params.TSS_species ?: "Ovis' 'aries",
            vep_species: params.vep_species ?: 'ovis_aries',
            vep_cache_species: params.vep_cache_species ?: 'Oar_v3.1',
			fasta_start: params.fasta_start ?: 'Ovis_aries'
        ]
    } else {
        throw new IllegalArgumentException("Unknown species and genome version. Please specify valid species and genome version.")
    }

    // Check if any of the required parameters are still null
    if (!params.CpG_species || !params.TSS_species || !params.vep_species) {
        throw new IllegalArgumentException(
            "Missing required parameters: TSS_species, vep_species, or CpG_species. " +
            "Please specify these parameters in your Nextflow config file or via command line."
        )
    }
}

process data_downloading {

    input:
    val par

    output:
    path "${par['CpG_species']}_download_complete.txt", emit: download_complete_ch

    script:
        // Define variables
        def fastaUrl = "https://ftp.ensembl.org/pub/release-${params.version}/fasta/${par['vep_species']}/dna/${par['fasta_start']}.${par['vep_cache_species']}.dna.toplevel.fa.gz"
        def vepCacheUrl = "https://ftp.ensembl.org/pub/release-${params.version}/variation/vep/${par['vep_species']}_vep_${params.version}_${par['vep_cache_species']}.tar.gz"
        def vepCacheFilename = vepCacheUrl.tokenize('/')[-1]
        def cache_dir = "$baseDir/data/vep_cache/${par['vep_species']}/${params.version}_${par['vep_cache_species']}"

        """
        # Create the data directory if it doesn't exist
        if [ ! -d "${baseDir}/data" ]; then
            echo "data directory not found, creating..."
            mkdir -p ${baseDir}/data
        fi

        # Download and prepare the genome file
        
        if [ ! -f "${baseDir}/data/${par['CpG_species']}.fa.gz" ]; then
            echo "Reference file not found, downloading..."
            wget -O ${baseDir}/data/${par['CpG_species']}.fa.gz ${fastaUrl}
        else
            echo "Reference file already exists, skipping download."
        fi

        if [ -f "${baseDir}/data/${par['CpG_species']}.fa.gz" ] && [ ! -f "${baseDir}/data/${par['CpG_species']}.fa.gz.fai" ]; then
            echo "Recompressing genome.fa.gz with bgzip..."
            gunzip -c ${baseDir}/data/${par['CpG_species']}.fa.gz | bgzip -c > ${baseDir}/data/${par['CpG_species']}.fa.bgz
            echo "Indexing the genome.fa.bgz..."
            samtools faidx ${baseDir}/data/${par['CpG_species']}.fa.bgz
            mv ${baseDir}/data/${par['CpG_species']}.fa.bgz ${baseDir}/data/${par['CpG_species']}.fa.gz
            mv ${baseDir}/data/${par['CpG_species']}.fa.bgz.fai ${baseDir}/data/${par['CpG_species']}.fa.gz.fai
        else
            echo "Index file already exists, skipping indexing."
        fi

        # Download the VEP cache if it doesn't exist
        if [ ! -d ${cache_dir} ]; then
            echo "vep cache version ${params.version} not found, downloading..."
            mkdir -p $baseDir/data/vep_cache
            curl -O ${vepCacheUrl}
            tar -xzf ${vepCacheFilename} -C ${baseDir}/data/vep_cache
            rm ${vepCacheFilename}
        else
            echo "vep cache already exists, skipping download."
        fi

        # Download conservation files based on species
        if [ "${params.species}" == "human" ]; then
            if [ ! -f "${baseDir}/data/hg38.phastCons100way.bw" ]; then
                wget -O ${baseDir}/data/hg38.phastCons100way.bw ${params.phastCons100wayUrl}
            fi
            if [ ! -f "${baseDir}/data/hg38.phastCons30way.bw" ]; then
                wget -O ${baseDir}/data/hg38.phastCons30way.bw ${params.phastCons30wayUrl}
            fi
            if [ ! -f "${baseDir}/data/hg38.phyloP100way.bw" ]; then
                wget -O ${baseDir}/data/hg38.phyloP100way.bw ${params.phyloP100wayUrl}
            fi
            if [ ! -f "$baseDir/data/hg38.phyloP30way.bw" ]; then
                wget -O ${baseDir}/data/hg38.phyloP30way.bw ${params.phyloP30wayUrl}
            fi
        elif [ "${params.species}" == "cattle" ]; then
			if [ ! -d "${baseDir}/data/cattle_conservation_scores" ]; then
            	wget ${params.cattle_conservation_scores}
				mv 'cattle_conservation_scores.zip?download=1' cattle_conservation_scores.zip
				unzip -o cattle_conservation_scores.zip -d ${baseDir}/data -x "__MACOSX/*" "*.DS_Store"
				rm cattle_conservation_scores.zip
			fi
        else
            echo "Conservation bigwig files unavailable for this species."
            exit 1
        fi

        touch ${par['CpG_species']}_download_complete.txt
		"""
}

/*
 * Process 2: Convert the input file to the sorted BED format file and vep input format.
 */

process formatting {

	input:
	path raw_data

	output:
	path 'data_bed_sorted.bed', emit: bed_ch
	path 'data_sorted.vep', emit: vep_ch
	path 'data_enformer_input', emit: enformer_ch

	script:

	if(raw_data.name =~ /.*[.]vcf$/ || raw_data.name =~ /.*[.]vcf[.]gz$/ || raw_data.name =~ /.*[.]bcf$/ || raw_data.name =~ /.*[.]bcf[.]gz$/)
	  """
	  bcftools query -f'%CHROM\t%POS0\t%END\t%CHROM_%POS_%REF_%ALT\n' ${raw_data} | sort -k1,1 -k2,2n -o data_bed_sorted.bed
	  python $baseDir/modules/script/getVepEnformerFormat.py data_bed_sorted.bed data_sorted.vep data_enformer_input
	  """

	else if(raw_data.name =~ /.*[.]bed$/ || raw_data.name =~ /.*[.]txt$/)
	  """
      sort -k1,1 -k2,2n ${raw_data} -o data_bed_sorted.bed
	  python $baseDir/modules/script/getVepEnformerFormat.py ${raw_data} data_sorted.vep data_enformer_input
	  """

	else
	  throw new IllegalArgumentException("Unknown format, please use bed, vcf, vcf.gz, bcf or bcf.gz file as input.")
	
	stub:
	"""
	touch data_bed_sorted.bed
	"""
}

/*
 * Process 3: Get the flanking sequence of the varaints.
 */
 
process flankingSequence {

	input:
	val par
	file data_bed_sorted
	path download_signal

	output:
	file 'flanking_annotation_result.txt'

	script:
	"""
	python $baseDir/modules/script/getFlankingSequence.py ${data_bed_sorted} ${params.flanking_len} ${par['CpG_species']}.fa.gz flanking_annotation_result.txt
	"""

	stub:
	"""
	touch flanking_annotation_result.txt
	"""
}

/*
 * Process 4: Get distances to different classes of genes (TSS)
 */

process distanceToTSS {

	input:
	val par
	file data_bed_sorted 

	output:
	file 'TSS_annotation_result.txt'

	script:
	"""
	Rscript $baseDir/modules/script/distanceToGenes.R ${par['TSS_species']} ${params.dataProvider} ${params.version} ${data_bed_sorted} TSS_annotation_result.txt
	"""

	stub:
	"""
	touch TSS_annotation_result.txt
	"""
}

/*
 * Process 5: Get distance to CpG island
 */

process distanceToCpG {

	input:
	val par 
	file data_bed_sorted 

	output:
	file 'CpG_annotation_result.txt'

	script:
	"""
    wget -O- http://hgdownload.cse.ucsc.edu/goldenpath/${par['CpG_species']}/database/cpgIslandExt.txt.gz | gunzip -c | awk 'BEGIN{ OFS="\t"; }{ print \$2, \$3, \$4, \$5\$6 }' | sort -k1,1 -k2,2n > CpG.bed
    bedtools closest -a ${data_bed_sorted} -b CpG.bed -t first -D b | awk ' BEGIN{ print "chrom", "start", "end", "variant_id", "CpG_id", "distance_to_CpG"} NF && NF-1 { print \$1, \$2, \$3, \$4, \$(NF - 1), \$NF}' > CpG_annotation_result.txt 
    """

	stub:
	"""
	touch CpG_annotation_result.txt
	"""
}

/*
 * Process 6: Conservation score (human and cattle)
 */

process conservationScore {
 	
	input:
	file data_bed_sorted
	path download_signal

	output:
	file 'conservation_result.txt'

	when:
	params.species == 'human' || params.species == 'cattle'

	script:
	if (params.species == 'human')
		"""
		python $baseDir/modules/script/conservationScore.py human ${data_bed_sorted} conservation_result.txt --phastCons100way_path $baseDir/data/hg38.phastCons100way.bw --phyloP100way_path $baseDir/data/hg38.phyloP100way.bw --phastCons30way_path $baseDir/data/hg38.phastCons30way.bw --phyloP30way_path $baseDir/data/hg38.phyloP30way.bw
		"""
	else if (params.species == 'cattle')
		"""
		python $baseDir/modules/script/conservationScore.py cattle ${data_bed_sorted} conservation_result.txt --phastCons241way_path $baseDir/data/cattle_conservation_scores/phastCons_bigwig --phyloP241way_path $baseDir/data/cattle_conservation_scores/phyloP_bigwig
		"""
	else 
		error "Conservation bigwig file unavailable for this species."

	stub:
	"""
	touch conservation_result.txt
	"""
}

/*
 * Process 7: Distance to Chromatin data (for human only)
 */

process distanceToChromatin {

 	input:
 	file data_bed_sorted

 	output:
 	file 'chromatin.txt'

 	when:
 	params.species == 'human'

 	script:
 	"""
 	awk '{print \$NF}' $baseDir/data/bedFiles.txt > path.txt
 	cat ${data_bed_sorted} | awk ' BEGIN{ print "chrom", "start", "end", "variant_id"}{ print \$1, \$2, \$3, \$4 }' > chromatin.txt
 	cat path.txt | while read LINE
 	do
 	wget -O- \$LINE | gunzip -c | sort -k1,1 -k2,2n > chromatin.bed
 	bedtools closest -a ${data_bed_sorted} -b chromatin.bed -t first -D b | awk -v array=\$LINE '{OFS="\\t"} BEGIN{OFS = "_"; split(array, a, "/"); print a[9], a[8], "enriched_site"}{ print \$NF}' > temp.txt
 	awk 'NR == FNR {a[NR] = \$0;next} {print a[FNR],\$0}' chromatin.txt temp.txt > result.txt
 	cat result.txt > chromatin.txt
 	done
 	"""

	stub:
	"""
	touch chromatin.txt
	"""
}

/*
 * Process 8: Gene density
 */

process geneDensity {

	input:
	val par 
	file data_bed_sorted 

	output:
	file 'density.txt'

	script:
	"""
	Rscript $baseDir/modules/script/geneDensity.R ${par['TSS_species']} ${params.dataProvider} ${params.version} ${data_bed_sorted} density.txt
	"""

	stub:
	"""
	touch density.txt
	"""
}


/*
 * Process 9: Distance to regulatory features (enhancers, promoters, CTCF binding sites, TF binding sites, human only)
 */

process distanceToRegulatoryFeatures {

  	input:
  	file data_bed_sorted

  	output:
  	file 'regulatory_features.txt'

  	when:
  	params.species == 'human'

  	script:
  	"""
  	Rscript $baseDir/modules/script/distanceToregulatoryFeatures.R ${data_bed_sorted} regulatory_features.txt $baseDir/data/enhancer.bed $baseDir/data/TF.bed $baseDir/data/CTCF.bed $baseDir/data/promoter.bed
  	"""

	stub:
	"""
	touch regulatory_features.txt
	"""
}

/*
 * Process 10: Get VEP annotations: consequence 
 */
def inputFileName = params.path.split('/')[-1]
def suffix = inputFileName.replaceFirst(/\.[^.]+$/, '')
def outputFileName = "vep_${suffix}.txt"
process getVep {
	publishDir 'Annotation_result', mode: 'copy'

	input:
	val par
	file data_sorted_vep
	path download_signal

	output:
	file "vep_${suffix}.txt"

	script:
	
	if(params.species == 'human')
	  """
	  vep -i ${data_sorted_vep} -o ${outputFileName} --cache --dir_cache $baseDir/data/vep_cache --no_stats --pick --af --af_1kg --regulatory --biotype --symbol --fork ${task.cpus} --no_check_variants_order --species ${par['vep_species']}
	  """
	else
	  """
      vep -i ${data_sorted_vep} -o ${outputFileName} --cache --dir_cache $baseDir/data/vep_cache --no_stats --pick --regulatory --biotype --symbol --fork ${task.cpus} --no_check_variants_order --species ${par['vep_species']}
	  """
}


