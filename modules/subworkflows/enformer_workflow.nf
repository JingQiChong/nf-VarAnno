#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
 
nextflow.enable.dsl = 2

//params.path = "$baseDir/*.vcf"
//params.enformer_prefix = 'temp'

include{
    getVariantScore
} from '../processes/enformer_processes.nf'

include{
    formatting
} from '../processes/main_annotation_processes.nf'

workflow enformer {
    formatting(params.path)
    getVariantScore(params.enformer_prefix, formatting.out.enformer_ch)
}
