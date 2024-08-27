#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info "============================================================"
log.info "          __     __      __                                 "
log.info "        / _|    \\ \\    / /         /\\                      "
log.info "  _ __ | |_ _____\\ \\  / /_ _ _ __ /  \\   _ __  _ __   ___  "
log.info " | '_ \\|  _|______\\ \\/ / _` | '__/ /\\ \\ | '_ \\| '_ \\ / _ \\ "
log.info " | | | | |         \\  / (_| | | / ____ \\| | | | | | | (_) |"
log.info " |_| |_|_|          \\/ \\__,_|_|/_/    \\_\\_| |_|_| |_|\\___/ "
log.info "                                                           "
log.info "============================================================"

if ( !params.path ) { 
    log.error "Please provide the target file path."
    exit 1 
}
if ( !params.ref_file ) { 
    log.error "Please provide the reference genome."
    exit 1 
}
if ( !params.species ) { 
    log.error "Please provide the species."
    exit 1 
}
if (params.species !in ['human', 'cattle', 'pig', 'dog']) {
    log.error "Please provide a valid species. Supported species are: human, cattle, pig, dog."
    exit 1 
}
if ( !params.genome_version ) { 
    log.error "Please provide the genome assembly version."
    exit 1 
}
if ( !params.function ) { 
    log.error "Please provide the type of annotations."
    exit 1 
}

log.info "Species: $params.species"
log.info "Genome version: $params.genome_version"
log.info "Function: $params.function"

include {main_annotation} from './modules/subworkflows/main_annotation_workflow' params(params)
include {enformer} from './modules/subworkflows/enformer_workflow' params(params)

workflow {
    if(params.function == 'general_annotation'){
        main_annotation()
    } else if (params.function == 'enformer_scores'){
        enformer()
    }
}

