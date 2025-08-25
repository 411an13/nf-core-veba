/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQPREPROCESSOR }      from '../subworkflows/local/fastqpreprocessor/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VEBA {

    take:
      ch_samplesheet // channel: samplesheet read in from --input

    main:
      FASTQPREPROCESSOR(
        reads                = ch_samplesheet,
        adapter_fasta        = params.adapter_fasta,
        discard_trimmed_pass = params.discard_trimmed_pass,
        save_trimmed_fail    = params.save_trimmed_fail,
        save_merged          = params.save_merged
      )

    emit:
      // forward useful outputs (extend as needed)
      reads_out       = FASTQPREPROCESSOR.out.reads_out
      aligned_bam     = FASTQPREPROCESSOR.out.aligned_bam
      multiqc_report  = FASTQPREPROCESSOR.out.multiqc_report
      versions        = FASTQPREPROCESSOR.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
