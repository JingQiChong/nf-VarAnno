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

/*
 * Process 2: Download data
 */
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
        # Ensure necessary tools are available
        command -v wget >/dev/null 2>&1 || { echo >&2 "wget is required but not installed. Aborting."; exit 1; }
        command -v curl >/dev/null 2>&1 || { echo >&2 "curl is required but not installed. Aborting."; exit 1; }
        command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools is required but not installed. Aborting."; exit 1; }
        command -v gunzip >/dev/null 2>&1 || { echo >&2 "gunzip is required but not installed. Aborting."; exit 1; }
        command -v bgzip >/dev/null 2>&1 || { echo >&2 "bgzip is required but not installed. Aborting."; exit 1; }

        # Create the data directory if it doesn't exist
        if [ ! -d "$baseDir/data" ]; then
            echo "Data directory not found, creating..."
            mkdir -p $baseDir/data
        fi

        # Download and prepare the genome file
        if [ ! -f $baseDir/data/${par['CpG_species']}.fa.gz ]; then
            echo "Reference file not found, downloading..."
            wget --continue --tries=5 --no-verbose -O $baseDir/data/${par['CpG_species']}.fa.gz ${fastaUrl}
            gunzip -c $baseDir/data/${par['CpG_species']}.fa.gz | bgzip -c > $baseDir/data/${par['CpG_species']}.fa.bgz
            samtools faidx $baseDir/data/${par['CpG_species']}.fa.bgz
            mv $baseDir/data/${par['CpG_species']}.fa.bgz $baseDir/data/${par['CpG_species']}.fa.gz
            mv $baseDir/data/${par['CpG_species']}.fa.bgz.fai $baseDir/data/${par['CpG_species']}.fa.gz.fai
        else
            echo "Reference file already exists, skipping download."
        fi

        # Download the VEP cache if it doesn't exist
        if [ ! -d ${cache_dir} ]; then
            echo "VEP cache version ${params.version} not found, downloading..."
            mkdir -p $baseDir/data/vep_cache
            wget --continue --no-verbose -O $baseDir/data/vep_cache/${vepCacheFilename} ${vepCacheUrl}
            tar -xzf $baseDir/data/vep_cache/${vepCacheFilename} -C $baseDir/data/vep_cache
            rm $baseDir/data/vep_cache/${vepCacheFilename}
        else
            echo "VEP cache already exists, skipping download."
        fi

        # Download conservation files based on species
        if [ "${params.species}" == "human" ]; then
            if [ ! -f $baseDir/data/hg38.phastCons100way.bw ]; then
                wget --continue --no-verbose -O $baseDir/data/hg38.phastCons100way.bw ${params.phastCons100wayUrl}
            fi
            if [ ! -f $baseDir/data/hg38.phastCons30way.bw ]; then
                wget --continue --no-verbose -O $baseDir/data/hg38.phastCons30way.bw ${params.phastCons30wayUrl}
            fi
            if [ ! -f $baseDir/data/hg38.phyloP100way.bw ]; then
                wget --continue --no-verbose -O $baseDir/data/hg38.phyloP100way.bw ${params.phyloP100wayUrl}
            fi
            if [ ! -f $baseDir/data/hg38.phyloP30way.bw ]; then
                wget --continue --no-verbose -O $baseDir/data/hg38.phyloP30way.bw ${params.phyloP30wayUrl}
            fi
        elif [ "${params.species}" == "cattle" ]; then
            if [ ! -d $baseDir/data/cattle_conservation_scores ]; then
                wget --continue --no-verbose -O $baseDir/data/cattle_conservation_scores.zip ${params.cattle_conservation_scores}

                unzip -o $baseDir/data/cattle_conservation_scores.zip -d $baseDir/data/
                if [ ! -d $baseDir/data/cattle_conservation_scores ]; then
                    echo "Error: Unzipping conservation scores failed."
                    exit 1
                fi
                rm $baseDir/data/cattle_conservation_scores.zip
            fi
        else
            echo "Conservation bigwig files unavailable for this species."
            exit 1
        fi

        # Create a signal file to indicate the process has completed successfully
        touch ${par['CpG_species']}_download_complete.txt
        """
}


/*
 * Process 3: Convert the input file to the sorted BED format file and vep input format.
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
