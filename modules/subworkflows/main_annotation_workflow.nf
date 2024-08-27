#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
 
nextflow.enable.dsl = 2


/*
 * Import modules
 */

include {
 	flankingSequence;
 	distanceToTSS;
 	distanceToCpG;
 	conservationScore;
 	distanceToChromatin;
 	geneDensity;
    distanceToRegulatoryFeatures;
	getVep
} from '../processes/main_annotation_processes.nf'

include {
	getParams;
	data_downloading;
	formatting;
} from '../processes/pre_processing.nf'

/*
 * Main pipeline
 */

workflow main_annotation {
 	getParams()
	getParams.out.view()
	data_downloading(getParams.out)
 	formatting(params.path)
	formatting.out.bed_ch.splitText(by: params.split_len, file: true).set{chunks_ch}
 	flankingSequence(getParams.out, chunks_ch, data_downloading.out.download_complete_ch)
 	flankingSequence.out.collectFile(name: 'flanking_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
 	distanceToTSS(getParams.out, chunks_ch)
 	distanceToTSS.out.collectFile(name: 'TSS_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
 	distanceToCpG(getParams.out, chunks_ch)
 	distanceToCpG.out.collectFile(name: 'CpG_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
 	conservationScore(chunks_ch, data_downloading.out.download_complete_ch)
 	conservationScore.out.collectFile(name: 'conservation_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
 	distanceToChromatin(chunks_ch)
    distanceToChromatin.out.collectFile(name: 'chromatin_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
	geneDensity(getParams.out, chunks_ch)
 	geneDensity.out.collectFile(name: 'density_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
	distanceToRegulatoryFeatures(chunks_ch)
 	distanceToRegulatoryFeatures.out.collectFile(name: 'regulatory_' + params.path.split('/')[-1], keepHeader: true, skip: 1, sort: true, storeDir: 'Annotation_result')
	getVep(getParams.out, formatting.out.vep_ch, data_downloading.out.download_complete_ch)
}




