#!/usr/bin/env nextflow
// Create fastq channel from samples.tsv

samples = file("samples.tsv")
aligned_reads = Channel
  .from(samples.readLines())
  .map {line ->
    lsplit      = line.split()
    sampleID    = lsplit[0]
    bam  = file(lsplit[1])
    bai  = file(lsplit[2])
    [ sampleID, bam, bai ]
}

aligned_reads.into{bam_files; bai_files}
bam_files = bam_files.map{sampleID, bam_file, bai_file -> bam_file}
bai_files = bai_files.map{sampleID, bam_file, bai_file -> bai_file}

all_bam_files = bam_files.toList()
all_bai_files = bai_files.toList()

all_bam_files.into{ all_bam_files_freebayes; all_bam_files_indexsplit}
all_bai_files.into{ all_bai_files_freebayes; all_bai_files_indexsplit}

process IndexSplit {
  publishDir "outputs/regions"

  cpus { 1 }
  memory { 32.GB }
//  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 8.GB}
  time { 2.d }
  errorStrategy { 'retry' }
  maxRetries 1
  maxErrors '-1'

  input:
  file(bams) from all_bam_files_indexsplit
  file(bais) from all_bai_files_indexsplit

  output:
  file("regions.bed") into regions_bed

  """
  /usr/local/bin/goleft indexsplit -n 500 *.bam > regions.bed
  """
}
counter=0
regions_channel = regions_bed.splitText{ line -> 
 lsplit = line.split()
 if (lsplit[4] != "0") {
   region = "--region " + lsplit[0] + ":" + lsplit[1] + "-" + lsplit[2]
   task = counter
   counter += 1
   return [task, region]
  }
 }.filter{ it != null}



process Freebayes1 {
publishDir "outputs/regionVCFs"
tag "$regionTask"

cpus { 1 }
memory { 8.GB }
//  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 8.GB}
time { 2.d }
errorStrategy { 'retry' }
maxRetries 0
maxErrors '-1'

input:
file(bams) from all_bam_files_freebayes
file(bais) from all_bai_files_freebayes
set regionTask, region from regions_channel
file params.genome
file params.genomeIndex

output:
set regionTask, file("region_${regionTask}.vcf") into region_VCFs

script:
input_bams = bams.collect{"-b $it"}.join(' ')

"""
freebayes ${params.freebayes_options} \
  -f $params.genome \
  $input_bams \
  $region > region_${regionTask}.vcf
"""
}

region_vcf_files = region_VCFs.map{id, file -> file}
all_VCFs = region_vcf_files.toList()

process ConcatenateVCFs {
publishDir "outputs"

memory { 16.GB * task.attempt}
time { 2.d }
errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
maxRetries 1
maxErrors '-1'

input:
file vcfs from all_VCFs

output:
file "${params.output_prefix}.raw.snps.indels.vcf"  into finalVCF

script:
input_vcfs = vcfs.collect{"$it"}.join(' ')

"""
vcf-concat *.vcf > ${params.output_prefix}.raw.snps.indels.vcf
"""
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
