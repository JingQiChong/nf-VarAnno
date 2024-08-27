#!/usr/bin/env nextflow

/*
 * Process 1: Get the 5-mer flanking sequence of the varaints.
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
 * Process 2: Get distances to different classes of genes (TSS)
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
 * Process 3: Get distance to CpG island
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
 * Process 4: Conservation score (human and cattle)
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
 * Process 5: Distance to Chromatin data (for human only)
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
 * Process 6: Gene density
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
 * Process 7: Distance to regulatory features (enhancers, promoters, CTCF binding sites, TF binding sites, human only)
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
 * Process 8: Get VEP annotations: consequence 
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


