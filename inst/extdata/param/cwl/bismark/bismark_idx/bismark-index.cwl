################################################################
##                        Bismark-Index                       ##
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

baseCommand: [ bismark_genome_preparation ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.idx_basedir) ]
    
arguments:
  - prefix: --bowtie2
  - valueFrom: $(inputs.idx_basedir.path)
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Basename of the bismark index files"
    type: Directory

outputs:
  index_files:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/Bisulfite_Genome/
