#!/usr/bin/env nextflow


process get_versions {
  echo true

  script:
  """
  echo -n "bwa ";  bwa 2>&1 | grep 'Version'
  echo -n "samtools ";  samtools 2>&1 | grep 'Version'
  fastqc --version
  multiqc --version
  """
}