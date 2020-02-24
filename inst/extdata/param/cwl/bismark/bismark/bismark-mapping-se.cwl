################################################################
##                       Bismark-mapping                      ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[Bismark]()"
label: Last updated 02/2020
hints:
  SoftwareRequirement:
    packages:
    - package: bismark
      version: [ v0.22.3 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [bismark]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --bowtie2
  - prefix: -N
    valueFrom: '0'
  - prefix: -L
    valueFrom: '20'
  - valueFrom: $(inputs.idx_basedir.path)
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).fastq.gz.out

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      prefix: 
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  bismark_mapping:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).fastq.gz.out
