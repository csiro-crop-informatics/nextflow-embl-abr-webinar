#!/usr/bin/env nextflow

/*
Locate the reference file, emit it through the referencesChannel.
We expect only a single reference but this syntax
allows to have muliple if needed.
*/
Channel.fromPath("data/references/reference.fasta.gz").set{ referencesChannel }

process bwa_index {
  input:
    file(ref) from referencesChannel

  output:
    set val("${ref}"), file("*") into indexChannel //also valid: set val(ref.name), file("*") into indexChannel

  script:
  """
  bwa index -a bwtsw ${ref}
  """
}

/*
Locate FATSQ files, emit each read file seperately
- The .take(n) operator allows us to limit how many files
  are emitted during a given run of the workflow.
- The .first() and .last() operatorsa are redundant given .take(1)
*/
Channel.fromPath("data/raw_reads/*.fastq.gz")
  .take(1)
  .first()
  .last()
  .set { readsForQcChannel }

process fastqc {
  tag { reads.name }
  input:
    file(reads) from readsForQcChannel

  output:
    file('*') into fastqcRawResultsChannel

  script:
  """
  fastqc  --quiet --threads ${task.cpus} ${reads}
  """
}