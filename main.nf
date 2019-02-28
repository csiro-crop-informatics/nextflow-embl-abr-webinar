#!/usr/bin/env nextflow

reference = Channel.fromPath(params.ref.base_url + params.ref.chr + ".fsa.zip")
accessionsChannel = Channel.from(params.accessions)

process bgzip_chromosome {
	input:
		file ref from reference

  output:
    file('*') into chromosomesChannel

	script:
  """
		unzip -p ${ref}
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
    // | samtools fastq -F 0x900 -1 R1.fastq.gz -2 R2.fastq.gz -s /dev/null -0 /dev/null -
}

process fastqc_raw {
  tag { accession }
	input:
		set val(accession), file('*') from extractedReadsChannel1

	output:
		file('*') into fastqcRawResultsChannel

  script:
   	"""
		fastqc --threads ${task.cpus} *
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


adaptersSingletonChannel = Channel.fromPath(params.adpaters)

process trimmomatic_pe {
	input:
		file('*') from
		file(adapaters) from adaptersSingletonChannel

	// output:
	// 	r1          = "qc_reads/{accession}_R1.fastq.gz",
	// 	r2          = "qc_reads/{accession}_R2.fastq.gz",
	// 	# reads where trimming entirely removed the mate
	// 	r1_unpaired = "qc_reads/{accession}_R1.unpaired.fastq.gz",
	// 	r2_unpaired = "qc_reads/{accession}_R2.unpaired.fastq.gz",
	script:
	"""
	  trimmomatic \
    * \
    r1_paired.gz \
		r1_unpaired.gz \
		r2_paired.gz \
		r2_unpaired.gz \
		ILLUMINACLIP:${adapters}:2:30:10:3:true \
	  LEADING:2 \
	  TRAILING:2 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36 \
	"""
}

// rule fastqc_trimmed:
// 	input:
// 		"qc_reads/{prefix}.fastq.gz",
// 	output:
// 		zip  = "reports/qc_reads/{prefix}_fastqc.zip",
// 		html = "reports/qc_reads/{prefix}_fastqc.html",
// 	conda:
// 		"envs/tutorial.yml"
// 	benchmark:
// 		repeat("benchmarks/fastqc_trimmed/{prefix}.txt", N_BENCHMARKS),
// 	wrapper:
// 		"0.31.1/bio/fastqc"
//   script:
// #		"""
// #		fastqc {input}
// #		"""

// rule multiqc_trimmed:
// 	input:
// 		expand("reports/qc_reads/{accession}_R{read}_fastqc.zip", accession=ACCESSIONS, read=[1,2]),
// 	output:
// 		"reports/qc_reads_multiqc.html",
// 	conda:
// 		"envs/tutorial.yml"
// 	benchmark:
// 		repeat("benchmarks/multiqc_trimmed/benchmark.txt", N_BENCHMARKS),
// 	wrapper:
// 		"0.31.1/bio/multiqc"
//   script:
// #		"""
// #		multiqc --filename {output} {input}
// #		"""


process bwa_index {
	input:
		file(ref) from subregionsChannel

	output:
		file("*") into indexChannel

  script:
	"""
		bwa index -a bwtsw ${ref}
	"""
}



// rule bwa_mem:
// 	input:
// 		reference = expand(REFERENCE + ".{ext}", ext=["amb","ann","bwt","pac","sa"]),
// 		reads     = [ "qc_reads/{sample}_R1.fastq.gz", "qc_reads/{sample}_R2.fastq.gz"],
// #		r1        = "qc_reads/{sample}_R1.fastq.gz",
// #		r2        = "qc_reads/{sample}_R2.fastq.gz",
// 	output:
// 		"mapped/{sample}.bam"
// 	conda:
// 		"envs/tutorial.yml"
// 	params:
// 		index      = REFERENCE,
// 		extra      = r"-R '@RG\tID:{sample}\tSM:{sample}'",
// 		sort       = "none",
// 		sort_order = "queryname",
// 		sort_extra = ""
// 	threads:
// 		MAX_THREADS
// 	benchmark:
// 		repeat("benchmarks/bwa_mem/{sample}.txt", N_BENCHMARKS),
// 	wrapper:
// 		"0.31.1/bio/bwa/mem"
