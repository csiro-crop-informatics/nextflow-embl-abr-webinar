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