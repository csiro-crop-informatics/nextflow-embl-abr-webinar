#!/usr/bin/env nextflow


process get_versions {
  echo true

  script:
  """
  fastqc --version
  multiqc --version
  echo -n "bwa " &&  bwa 2>&1 | grep 'Version'
  echo -n "samtools "&&  samtools 2>&1 | grep 'Version'
  """
}