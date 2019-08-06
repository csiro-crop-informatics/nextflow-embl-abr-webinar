#!/usr/bin/env nextflow

region = "${params.chr}_${params.start}-${params.end}"

referencesChannel = Channel.fromPath("data/${region}/*.fasta.gz")

//fetch adapters file - either local or remote
adaptersChannel = Channel.fromPath(params.adapters) //NF will download if remote

Channel.fromFilePairs("data/${region}/*_R{1,2}.fastq.gz")
.take( params.take == 'all' ? -1 : params.take ) //Use --take N to process first N accessions or --take all to process all
.into { extractedReadsChannelA; extractedReadsChannelB }

process bwa_index {
  tag { ref }
  input:
    file(ref) from referencesChannel

  output:
    set val(ref.name), file("*") into indexChannel //also valid: set val("${ref}"), file("*") into indexChannel

  script:
  """
  bwa index -a bwtsw ${ref}
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

process bwa_mem {
  tag { accession }

  input:
    set val(ref), file('*'), val(accession), file(reads) from indexChannel.combine(trimmedReadsChannelA)

	output:
		file('*.bam') into alignedReadsChannel

  script:
  """
  bwa mem -t ${task.cpus} -R '@RG\\tID:${accession}\\tSM:${accession}' ${ref} ${reads} | samtools view -b > ${accession}.bam
  """
}