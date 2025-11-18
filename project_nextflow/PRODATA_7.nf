params.reads = null
params.outdir = "results"
params.fastqc_threads = 4
params.multiqc_title = 'Analiza MultiQC'
params.trimmed = "trimmed"

process FASTQC {
    tag 'fastqc'
    publishDir "$params.outdir/quality_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path '*_fastqc.zip', emit: fastqc_zip
    path '*_fastqc.html', emit: fastqc_html

    script:
    """
    fastqc ${read1} ${read2} -o . -t $params.fastqc_threads
    """
}

process MULTIQC {
    tag 'multiqc'
    container 'multiqc/multiqc:v1.29'
    publishDir "$params.outdir/quality_reports_multiqc", mode: 'copy'

    input:
    path fastqc_zip

    output:
    path '*multiqc_report.html', emit: multiqc_html
    path '*multiqc_report_data', emit: multiqc_data

    script:
    """
    multiqc . --title '$params.multiqc_title'
    """
}
process TRIMMED {
    tag 'trimmed'
    publishDir "$params.outdir/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple path("trimmed_${read1.name}"), path("trimmed_${read2.name}")

    script:
    """
    trimmomatic PE -threads $params.fastqc_threads \
      ${read1} ${read2} \
      trimmed_${read1.name} unpaired_${read1.name} \
      trimmed_${read2.name} unpaired_${read2.name} \
      SLIDINGWINDOW:4:20 MINLEN:50
    """
}

workflow {
    log.info """\
    ====================================================
     F A S T Q C   +   M U L T I Q C   +   T R I M M E D
    ====================================================
    input           : ${params.input}
    outdir          : ${params.outdir}
    fastqc_threads  : ${params.fastqc_threads}
    multiqc_title   : ${params.multiqc_title}
    """
    .stripIndent()

    if (!params.input) {
        log.error "Input files are not specified. Please provide input files using the --input parameter."
        exit 1
    }

    input_ch = Channel.fromFilePairs(params.input, flat: true)

    FASTQC(input_ch)
    MULTIQC(FASTQC.out.fastqc_zip)
    TRIMMED(input_ch)
}
