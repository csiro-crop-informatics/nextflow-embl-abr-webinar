#!/usr/bin/env nextflow

//Build link to reference
referenceLink = params.ref.base_url + params.ref.chr + ".fsa.zip"

//Take accessions defined in nextflow.config.
//Use --take N to process first N accessions or --take all to process all
accessionsChannel = Channel.from(params.accessions).take( params.take == 'all' ? -1 : params.take )

//fetch adapters file - either local or remote
adaptersChannel = Channel.fromPath(params.adapters)

process download_chromosome {
  tag { params.ref.chr }

  //Prevent re-downloading of large files
  // storeDir { executor == 'awsbatch' ? null : "${params.outdir}/downloaded" }  //use with care, caching will not work as normal so changes to input may not take effect
  storeDir { "${params.outdir}/downloaded" }  //use with care, caching will not work as normal so changes to input may not take effect
  scratch false //must be false otherwise storeDir ignored

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
  cpus '2' //consider defining in conf/requirements.config based on process name or label
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
    file chr from chromosomesChannel

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
    //e.g. ACBarrie

  output:
    set val(accession), file('*.fastq.gz') into (extractedReadsChannelA, extractedReadsChannelB)
    //e.g. ACBarrie.realigned.bam.bai, ACBarrie_R1.fastq.gz, ACBarrie_R2.fastq.gz

  script:
  """
  samtools view -hu "${params.bam.base_url}/chr4A_part2/${accession}.realigned.bam" \
    ${params.bam.chr}:${params.bam.start}-${params.bam.end} \
  | samtools collate -uO - \
  | samtools fastq -F 0x900 -1 ${accession}_R1.fastq.gz -2 ${accession}_R2.fastq.gz \
    -s /dev/null -0 /dev/null -
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

process bwa_index {
  input:
    file(ref) from subregionsChannel

  output:
    set val(ref.name), file("*") into indexChannel //also valid: set val("${ref}"), file("*") into indexChannel

  script:
  """
  bwa index -a bwtsw ${ref}
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