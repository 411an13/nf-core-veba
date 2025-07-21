// ─────────────────────────────────────────────
// MODULE IMPORTS
// ─────────────────────────────────────────────
include { FASTP  } from '../../modules/nf-core/fastp'                  // Replaces the fastq_fastp subworkflow
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { MULTIQC } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap, softwareVersionsToYAML } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_fastqpreprocessor_pipeline'

// ─────────────────────────────────────────────
// SUBWORKFLOW
// ─────────────────────────────────────────────
workflow {

  take:
    reads                       // renamed from ch_samplesheet to match FASTP's input format
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
    // Direct call to FASTP replaces the fastq_fastp subworkflow
    FASTP(
      reads                = reads,
      adapter_fasta        = adapter_fasta,
      discard_trimmed_pass = discard_trimmed_pass,
      save_trimmed_fail    = save_trimmed_fail,
      save_merged          = save_merged
    )

    // Add FASTP output to MultiQC inputs and version tracker
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })
    ch_versions      = ch_versions.mix(FASTP.out.versions.first())

    // ─────────────────────
    // MODULE: FASTQC
    // ─────────────────────
    // Run FastQC on trimmed reads from FASTP
    FASTQC(FASTP.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions      = ch_versions.mix(FASTQC.out.versions.first())

    // ─────────────────────
    // COLLECT SOFTWARE VERSIONS
    // ─────────────────────
    // Combine all version info into a YAML file for MultiQC and reproducibility
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
    // Load config, custom logo, and summary files for MultiQC report
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

    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    // Add versions and methods to MultiQC input
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
      ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
    )

    // ─────────────────────
    // MODULE: MULTIQC
    // ─────────────────────
    // Build final MultiQC report
    MULTIQC(
      ch_multiqc_files.collect(),
      ch_multiqc_config.toList(),
      ch_multiqc_custom_config.toList(),
      ch_multiqc_logo.toList(),
      [],
      []
    )

  emit:
    // Outputs from FASTP
    reads_out       = FASTP.out.reads
    reads_fail      = FASTP.out.reads_fail
    reads_merged    = FASTP.out.reads_merged
    fastp_html      = FASTP.out.html
    fastp_json      = FASTP.out.json
    log             = FASTP.out.log

    // Versions and final MultiQC report
    versions        = ch_versions
    multiqc_report  = MULTIQC.out.report.toList()
}
