# hapTEins
A pipeline for transposable elements(TE) insertion detection at haplotype-level
## Usage
- clone this repository, then run script in your working directory
- change the `Python` interpreter path of python scripts
## The pipeline comprise three steps to detect TE insertions at haplotype-level
1. align NGS data to both haplotypes and perform reads phasing
2. detect TE insertions utilizing TE insertion signals for each haplotype 
3. generate IGV snapshot and manual check
## 1. align and reads phasing, `pipeline.sh`
This procedure comprises three steps:
1. align NGS reads to reference genome with `align-M.sh` and `align-P.sh`
   
   - **substitute the ${GATK} and ${java} as your software path**
   
   - **substitute the ${ref} and ${genome} as your reference genome path**
   
   - **generate *.fai* index by `samtools` and *.dict* index by `picard`**
3. extract reads alignment with `getReadsMateTag`
   
   - **a C script implemented with htslib**
5. reads phasing by comparing alignment score(AS) and mismatch number(NM) with `readsCompare.py`
6. align the phased reads to haplotypes appended with your target TEs, `align-M_TE.sh` and `align-P_TE.sh`
   
   - **same to step1, but substitute the ${ref} and ${genome} as your reference genome appended with target TEs, with each TE as individual sequences**

## 2. detect TE insertion for each haplotype, `pipeline-dectTE.sh`

This procedure comprises two steps:

1. get candidate insertion positions, `getPos.py`
   
   - **substitute variable TEs as your TE seqeuence IDs, same with the target TE IDs in step1.6**

   - **substitute variable mL as your TE sequence length**

2. filter candidate insertion position to get high condidant insertion, `fltDiscordant.py`

## 3. generate IGV snapshot

  - **substitute variables bamDir, bam, disAln and TEAlns as corresponding path**
