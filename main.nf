#!/usr/bin/env nextflow

//Take accessions defined in nextflow.config.
//Use --take N to process first N accessions or --take all to process all
accessionsChannel = Channel.from(params.accessions).take( params.take == 'all' ? -1 : params.take )

region = "${params.chr}_${params.start}-${params.end}"

//Reference fasta.gz specified in nextflow.config, override with --reference path_or_url.fasta.gz
reference = Channel.from(params.reference)

//fetch adapters file - either local or remote
adaptersChannel = Channel.fromPath(params.adapters) //NF will download if remote



process get_reference {
  tag { "${region}"}
  input:
    val(ref_url) from reference

  output:
    file('*') into referencesChannel

  script:
  """
  wget ${ref_url}
  """
}

process bwa_index {
  tag { ref }
  input:
    file(ref) from referencesChannel

  output:
    set val(ref), file("*") into indexChannel //also valid: set val("${ref}"), file("*") into indexChannel

  script:
  """
  bwa index -a bwtsw ${ref}
  """
}


process get_reads {
  tag { "${accession} @ ${region}"}

  input:
    val accession from accessionsChannel
    //e.g. ACBarrie

  output:
    set val(accession), file("${accession}_R?.fastq.gz") into (extractedReadsChannelA, extractedReadsChannelB)
    //e.g. ACBarrie, [ACBarrie_R1.fastq.gz, ACBarrie_R2.fastq.gz]

  script:
  URL_BASE = [params.reads_base_url, region, accession].join('/')
  """
  wget ${URL_BASE}_${params.r1_suffix} ${URL_BASE}_${params.r2_suffix} \
  && zcat ${accession}_R1.fastq.gz | head | awk 'END{exit(NR<4)}' \
  && zcat ${accession}_R2.fastq.gz | head | awk 'END{exit(NR<4)}'
  """
}


process fastqc_raw {
  tag { accession }

  input:
    set val(accession), file('*') from extractedReadsChannelA

  output:
    file('*') into fastqcRawResultsChannel

  script:
  """
  fastqc  --quiet --threads ${task.cpus} *
  """
}

process multiqc_raw {
  input:
    file('*') from fastqcRawResultsChannel.collect()

  output:
    file('*') into multiqcRawResultsChannel

  script:
  """
  multiqc .
  """
}

process trimmomatic_pe {
  echo true
  tag {accession}

  input:
    set file(adapters), val(accession), file('*') from adaptersChannel.combine(extractedReadsChannelB)

  output:
    set val(accession), file('*.paired.fastq.gz') into (trimmedReadsChannelA, trimmedReadsChannelB)

  script:
  """
  trimmomatic PE \
    *.fastq.gz \
    ${accession}_R1.paired.fastq.gz \
    ${accession}_R1.unpaired.fastq.gz \
    ${accession}_R2.paired.fastq.gz \
    ${accession}_R2.unpaired.fastq.gz \
    ILLUMINACLIP:${adapters}:2:30:10:3:true \
    LEADING:2 \
    TRAILING:2 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
  """
}

process fastqc_trimmed {
  tag { accession }

  input:
    set val(accession), file('*') from trimmedReadsChannelB

  output:
    file('*') into fastqcTrimmedResultsChannel

  script:
  """
  fastqc --quiet --threads ${task.cpus} *
  """
}

process multiqc_trimmed {
  input:
    file('*') from fastqcTrimmedResultsChannel.collect()

  output:
    file('*') into multiqcTrimmedResultsChannel

  script:
  """
  multiqc .
  """
}

// indexChannel.combine(trimmedReadsChannelA).view()
// [
//  reference.fasta.gz,
//  [
//    /home/rad/repos/nextflow-embl-abr-webinar/work/c0/c01aea6edd42b4f981d2e97718d0ed/reference.fasta.gz.amb,
//    /home/rad/repos/nextflow-embl-abr-webinar/work/c0/c01aea6edd42b4f981d2e97718d0ed/reference.fasta.gz.ann,
//    /home/rad/repos/nextflow-embl-abr-webinar/work/c0/c01aea6edd42b4f981d2e97718d0ed/reference.fasta.gz.bwt,
//    /home/rad/repos/nextflow-embl-abr-webinar/work/c0/c01aea6edd42b4f981d2e97718d0ed/reference.fasta.gz.pac,
//    /home/rad/repos/nextflow-embl-abr-webinar/work/c0/c01aea6edd42b4f981d2e97718d0ed/reference.fasta.gz.sa
//  ],
//  ACBarrie,
//  [
//    /home/rad/repos/nextflow-embl-abr-webinar/work/8f/f03aa5c708c5b39d27e6066d3a8338/ACBarrie_R1.paired.fastq.gz,
//    /home/rad/repos/nextflow-embl-abr-webinar/work/8f/f03aa5c708c5b39d27e6066d3a8338/ACBarrie_R2.paired.fastq.gz]
// ]

process bwa_mem {
  tag { accession }

  input:
    // set val(ref), file('*'), val(accession), file(reads) from indexChannel.combine(trimmedReadsChannelA)
    // set val(ref), file('*'), val(accession), file(reads) from indexChannel.combine(trimmedReadsChannelA)
    set val(accession), file(reads) from trimmedReadsChannelA
    set val(ref), file(index) from indexChannel

	// output:
	// 	file('*.bam') into alignedReadsChannel

  script:
  """
  bwa mem -t ${task.cpus} -R '@RG\\tID:${accession}\\tSM:${accession}' ${ref} ${reads} | samtools view -b > ${accession}.bam
  """
}