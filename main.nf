#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.path) {exit 1 "Please provide the target file path."}

//params.function = 'general_annotation'

include {main_annotation} from './modules/subworkflows/main_annotation_workflow' params(params)
include {enformer} from './modules/subworkflows/enformer_workflow' params(params)

workflow {
    if(params.function == 'general_annotation'){
        main_annotation()
    } else if (params.function == 'enformer_scores'){
        enformer()
    }
}

