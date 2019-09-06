#!/usr/bin/env nextflow

/*
Locate the reference file, emit it through the referencesChannel.
We expect only a single reference but this syntax
allows to have muliple if needed.
*/
Channel.fromPath("data/**.fasta.gz").set{ referencesChannel }

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
Locate paired FATSQ files, emit each pair seperately
in the form [common_prefix_string, [R1_file, R2_file]]
- The .take() operator allows us to limit how many pairs
  are emitted during a given run of the workflow.
- The .into{ } operator allows us to fork the newly created channel
*/
Channel.fromFilePairs("data/**_R{1,2}.fastq.gz")
  .take ( {
    if(params.take == 'all') {
      return -1 //.take() treats -1 as "let all emissions through"
    } else {
      return params.take // "let first n emissions through" - as specified using --take n runtime param
    }
  }() )
  // .take ( params.take == 'all' ? -1 : params.take ) //Alternative (ternary) syntax
  .into { readPairsForQcChannel; readPairsForTrimmingChannel } //send each item into two separate channels


process fastqc {
  tag { accession }

  input:
    set val(accession), file(reads) from readPairsForQcChannel

  output:
    file('*') into fastqcRawResultsChannel

  script:
  """
  fastqc  --quiet --threads ${task.cpus} ${reads}
  """
}