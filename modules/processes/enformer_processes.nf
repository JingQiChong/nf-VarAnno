#!/usr/bin/env nextflow

process getVariantScore {
    maxForks 1
    publishDir 'Annotation_result/Enformer'
    input:
    val prefix
    file data_enformer_input
    val par

    output:
    file "${prefix}_enformer_variant_score.txt"
    file "${prefix}_enformer_variant_score_abs.txt"
    file "${prefix}_enformer_variant_score_prop.txt"

    script:
    """
    wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/${par['CpG_species']}/bigZips/${par['CpG_species']}.fa.gz | gunzip -c > $baseDir/data/${par['CpG_species']}.fa
    python $baseDir/modules/script/getVariantScore.py ${data_enformer_input} ${prefix}_enformer_variant_score.txt ${prefix}_enformer_variant_score_abs.txt ${prefix}_enformer_variant_score_prop.txt $baseDir/data/${par['CpG_species']}.fa
    """
}

