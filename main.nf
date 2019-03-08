#!/usr/bin/env nextflow

referenceLink = params.ref.base_url + params.ref.chr + ".fsa.zip"
// reference = Channel.fromPath(params.ref.base_url + params.ref.chr + ".fsa.zip")
accessionsChannel = Channel.from(params.accessions).take( params.take == 'all' ? -1 : params.take )

process download_chromosome {
  tag { params.ref.chr }
  input:
    referenceLink

  output:
    file('*') into references

  script:
  """
  wget ${referenceLink}
  """
}

process bgzip_chromosome {
  cpus '4'
  tag { ref }
  input:
    file ref from references

  output:
    file('*') into chromosomesChannel

  script:
  """
  unzip -p ${ref} \
    | bgzip --threads ${task.cpus} \
    > ${ref}.gz
  """
}

process bgzip_chromosome_subregion {
  input:
    file(chr) from chromosomesChannel

  output:
    file('subregion') into subregionsChannel

  script:
  """
  samtools faidx ${chr} ${params.ref.chr}:${params.ref.start}-${params.ref.end} \
    | bgzip --threads ${task.cpus} \
    > subregion
  """
}

process extract_reads {
  tag { accession }
  input:
    val accession from accessionsChannel

  output:
    set val(accession), file('*.fastq.gz') into (extractedReadsChannel1, extractedReadsChannel2)
    // set val(accession), file('*.fastq.gz') into (extractedReadsChannel1, extractedReadsChannel2)
    // ACBarrie.realigned.bam.bai, raw_reads/ACBarrie_R1.fastq.gz, raw_reads/ACBarrie_R2.fastq.gz

  script:
  """
  samtools view -hu "${params.bam.base_url}/chr4A_part2/${accession}.realigned.bam" ${params.bam.chr}:${params.bam.start}-${params.bam.end} \
  | samtools collate -uO - \
  | samtools fastq -F 0x900 -1 ${accession}_R1.fastq.gz -2 ${accession}_R2.fastq.gz -s /dev/null -0 /dev/null -
  """

}

process fastqc_raw {
  echo true
  time 1.m
  tag { accession }
  input:
    set val(accession), file('*') from extractedReadsChannel1

  output:
    file('*') into fastqcRawResultsChannel

  script:
  """
  #ls -l
  #echo
  fastqc  --quiet --threads ${task.cpus} *
  #zcat * | head -4
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

adaptersSingletonChannel = Channel.fromPath(params.adapters)

process trimmomatic_pe {
  tag {accession}

  input:
    set val(accession), file('*') from extractedReadsChannel2
    file(adapters) from adaptersSingletonChannel

  output:
    set val(accession), file('*.paired.fastq.gz') into (trimmedReadsChannel1, trimmedReadsChannel2)

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
    set val(accession), file('*') from trimmedReadsChannel2

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

process bwa_index {
  input:
    file(ref) from subregionsChannel

  output:
    set val(ref.name), file("*") into indexChannel
    // set val("${ref}"), file("*") into indexChannel

  script:
  """
  bwa index -a bwtsw ${ref}
  """
}


process bwa_mem {
  tag { accession }
	input:
    set val(ref), file('*'), val(accession), file(reads) from indexChannel.combine(trimmedReadsChannel1)

	output:
		file('*.bam')

  script:
  """
  bwa mem -t ${task.cpus} -R '@RG\\tID:${accession}\\tSM:${accession}' ${ref} ${reads} | samtools view -b > ${accession}.bam
  """
}