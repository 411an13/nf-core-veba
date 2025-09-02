// ─────────────────────────────────────────────
// MODULE IMPORTS
// ─────────────────────────────────────────────
include { SEQKIT_STATS }  from '../../../modules/nf-core/seqkit/stats/main'
include { FASTP }         from '../../../modules/nf-core/fastp/main'
include { BOWTIE2_ALIGN } from '../../../modules/nf-core/bowtie2/align/main'
include { BBMAP_BBDUK }   from '../../../modules/nf-core/bbmap/bbduk/main'
include { FASTQC }        from '../../../modules/nf-core/fastqc/main'
include { MULTIQC }       from '../../../modules/nf-core/multiqc/main'

include { paramsSummaryMap } from 'plugin/nf-schema'

// NOTE: separate names with semicolons and use the correct path
include { paramsSummaryMultiqc; softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'

// include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_fastqpreprocessor_pipeline'

// ─────────────────────────────────────────────
// SUBWORKFLOW
// ─────────────────────────────────────────────
workflow FASTQPREPROCESSOR {

  take:
    reads                       // channel of tuples: usually [meta, reads]
    adapter_fasta
    discard_trimmed_pass
    save_trimmed_fail
    save_merged

  main:
    // Channels used to collect tool versions and MultiQC input files
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ─────────────────────
    // MODULE: FASTP
    // ─────────────────────
    FASTP(
      reads,
      adapter_fasta,
      discard_trimmed_pass,
      save_trimmed_fail,
      save_merged
    )

    // Add FASTP output to MultiQC inputs and version tracker
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })
    ch_versions      = ch_versions.mix(FASTP.out.versions.first())

    // SEQKIT_STATS on FASTP trimmed reads
    SEQKIT_STATS(FASTP.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats.collect { it[1] })
    ch_versions      = ch_versions.mix(SEQKIT_STATS.out.versions.first())

    // ─────────────────────
    // MODULE: BBMAP_BBDUK (Contaminant/adapter removal)
    // ─────────────────────
    // Decide if we run BBDUK (requires --bbduk_contaminants file)
def run_bbduk = params.bbduk_contaminants && params.bbduk_contaminants.toString().trim()

if (run_bbduk) {
    BBMAP_BBDUK(
        FASTP.out.reads,
        file(params.bbduk_contaminants, checkIfExists: true)
    )
    ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect { it[1] })
    ch_versions      = ch_versions.mix(BBMAP_BBDUK.out.versions.first())

    cleaned_reads = BBMAP_BBDUK.out.reads

    // Stats on cleaned reads
    SEQKIT_STATS(cleaned_reads)
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats.collect { it[1] })
    ch_versions      = ch_versions.mix(SEQKIT_STATS.out.versions.first())
} else {
    log.info "[FASTQPREPROCESSOR] Skipping BBDUK: --bbduk_contaminants not set."
    cleaned_reads = FASTP.out.reads
}

// ─────────────────────
// Optional: BOWTIE2_ALIGN (Host genome alignment)
// Requires --host_index at minimum
// ─────────────────────
def run_bt2 = params.host_index && params.host_index.toString().trim()
if (run_bt2) {
    BOWTIE2_ALIGN(
        cleaned_reads,        // use BBDUK-cleaned reads if run; else FASTP reads
        params.host_index,
        params.host_fasta,
        false,                // save_unaligned
        true                  // sort_bam
    )
    ch_versions      = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect { it[1] })
    aligned_bam = BOWTIE2_ALIGN.out.bam
} else {
    log.info "[FASTQPREPROCESSOR] Skipping Bowtie2: --host_index not set."
    aligned_bam = Channel.empty()
}

    // ─────────────────────
    // MODULE: FASTQC
    // ─────────────────────
    FASTQC(FASTP.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions      = ch_versions.mix(FASTQC.out.versions.first())

    // ─────────────────────
    // COLLECT SOFTWARE VERSIONS
    // ─────────────────────
    softwareVersionsToYAML(ch_versions)
      .collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_fastqpreprocessor_software_mqc_versions.yml',
        sort: true,
        newLine: true
      )
      .set { ch_collated_versions }

    // ─────────────────────
    // MULTIQC SETUP
    // ─────────────────────
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files    = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    // Use default or custom methods description
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
      file(params.multiqc_methods_description, checkIfExists: true) :
      file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    // Add versions and methods to MultiQC input
    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    //ch_multiqc_files = ch_multiqc_files.mix(
    //ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
    //)

    // ─────────────────────
    // MODULE: MULTIQC
    // ─────────────────────
    MULTIQC(
      ch_multiqc_files.collect(),
      ch_multiqc_config.toList(),
      ch_multiqc_custom_config.toList(),
      ch_multiqc_logo.toList(),
      [],
      []
    )

  emit:
    // FASTP outputs
    reads_out       = cleaned_reads           // post-BBDUK cleaned reads
    reads_fail      = FASTP.out.reads_fail
    reads_merged    = FASTP.out.reads_merged
    fastp_html      = FASTP.out.html
    fastp_json      = FASTP.out.json
    // log             = FASTP.out.log

    // Alignment output (optional)
    aligned_bam     = aligned_bam

    // Versions and final MultiQC report
    versions        = ch_versions
    multiqc_report  = MULTIQC.out.report.toList()
}
