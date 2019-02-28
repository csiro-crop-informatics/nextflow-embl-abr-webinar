#!/usr/bin/env nextflow

reference = Channel.fromPath(params.ref.base_url + params.ref.chr. + ".fsa.zip")

process bgzip_chromosome {
	input:
		file ref from reference

  output:
    file('*') into chromosomesChannel

	script:
  """
		unzip -p ${ref} \
		  | bgzip --threads ${task.cpus} > chr
	"""
	// > ${ref.take(filename.lastIndexOf('.')}.gz// """
}

process bgzip_chromosome_subregion {
	input:
		file(chr) from chromosomesChannel

	output:
		file(subregion) into subregionsChannel

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
    accession from Channel.from(params.accessions)

  output:
    set val(accession), file('*.fastq.gz') into extractedReadsChannel
    // ACBarrie.realigned.bam.bai, raw_reads/ACBarrie_R1.fastq.gz, raw_reads/ACBarrie_R2.fastq.gz

    script:
    """
      samtools view -hu "${params.bam.base_url}/chr4A_part2/${accession}.realigned.bam" ${params.bam.chr}:${params.bam.start}-${params.bam.end} \
		  | samtools collate -uO - \
		  | samtools fastq -F 0x900 -1 ${accession}_R1.fastq.gz -2 ${accession}_R2.fastq.gz -s /dev/null -0 /dev/null -
    """
    // | samtools fastq -F 0x900 -1 R1.fastq.gz -2 R2.fastq.gz -s /dev/null -0 /dev/null -
}

