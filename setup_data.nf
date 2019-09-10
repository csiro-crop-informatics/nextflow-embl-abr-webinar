#!/usr/bin/env nextflow

//Take accessions defined in nextflow.config.
//Use --take N to process first N accessions or --take all to process all
accessionsChannel = Channel.from(params.accessions).take( params.download == 'all' ? -1 : params.download )

region = "${params.chr}_${params.start}-${params.end}"

//Reference fasta.gz specified in nextflow.config, override with --reference path_or_url.fasta.gz
reference_url = params.reference


process get_reference { //alternative: referencesChannel = Channel.fromPath(params.reference)
  tag { "${region}"}
  publishDir "data/reference", mode: 'copy'

  input:
    reference_url

  output:
    file('*')

  script:
  """
  wget ${reference_url}
  """
}

process get_reads {
  tag { "${accession} @ ${region}"}
  publishDir "data/raw_reads", mode: 'copy'

  input:
    val accession from accessionsChannel
    //e.g. ACBarrie

  output:
    set val(accession), file("${accession}_R?.fastq.gz")
    //e.g. ACBarrie, [ACBarrie_R1.fastq.gz, ACBarrie_R2.fastq.gz]

  script:
  URL_BASE = [params.reads_base_url, region, accession].join('/')
  """
  wget ${URL_BASE}_${params.r1_suffix} ${URL_BASE}_${params.r2_suffix} \
  && zcat ${accession}_R1.fastq.gz | head | awk 'END{exit(NR<4)}' \
  && zcat ${accession}_R2.fastq.gz | head | awk 'END{exit(NR<4)}'
  """
}

