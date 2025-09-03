#Completed on 08/28/2025
#to run this script:
#cd /private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/
#conda activate snakemake
# snakemake --executor slurm --default-resources slurm_partition=medium runtime=720 mem_mb=1000000 -j 10 -s Snakefile

#Global Variables:

import os  
import glob

samples = ['20250828_Nanopore']
conda: '/private/groups/russelllab/jodie/bootcamp2024/scripts/read_filtering.yaml'   

rule all:
    input:
        expand('data/polished/{sample}_wRi_M23.assembly.fasta', sample=samples)

# ==============================================================================
# SHARED PREPROCESSING STEPS
# ==============================================================================

rule basecalling:
    input:
        pod5 = 'data/pod5/'
    output:
        bam = 'data/aligned/aligned_basecalled.bam'
    params:
        genome = 'data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz',
        outdir = 'data/aligned/'
    threads: 25
    resources:
        slurm_partition='gpu',
        mem_mb=100000,  # 100GB in MB
        runtime=360,    # 6 hours in minutes
        slurm_extra='--gpus-per-node=4 --nodes=1 --exclude=phoenix-09 --mail-user=jomojaco@ucsc.edu --mail-type=ALL --output=logs/%x.%j.log'
    shell:
        """
        # Dorado Basecalling:
        /private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado basecaller hac {input.pod5} \
            --device cuda:all \
            --reference {params.genome} \
            --kit-name SQK-NBD114-24 > {output.bam}
        """

rule sort_bam:
    input:
        bam = 'data/aligned/aligned_basecalled.bam'
    output:
        sorted_bam = 'data/aligned/aligned_basecalled.sorted.bam',
        index = 'data/aligned/aligned_basecalled.sorted.bam.bai'
    resources: 
        mem_mb=100000,
        runtime=200
    threads: 16  
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly

        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
        samtools index -@ {threads} {output.sorted_bam}  
        '''

rule separate_bam:
    input:
        bam = 'data/aligned/aligned_basecalled.sorted.bam',
    output:
        host_bam = 'data/aligned/{sample}.DsimM23.bam',
        wolbachia_bam = 'data/aligned/{sample}.wRiM23.bam'
    resources: 
        mem_mb=100000,
        runtime=200
    threads: 8  
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly

        samtools view -@ {threads} -b -f 4 {input.bam} > {output.wolbachia_bam}
        samtools view -@ {threads} -b -F 4 {input.bam} > {output.host_bam}
        '''

# ==============================================================================
# WRI (WOLBACHIA) PROCESSING  
# ==============================================================================

rule wri_bam2fastq:
    input:
        wolbachia_bam = 'data/aligned/{sample}.wRiM23.bam'
    output:
        wolbachia_fastq = 'data/basecalled/{sample}.wRiM23.fastq.gz'
    resources: 
        mem_mb=25000,
        runtime=60
    threads: 8  
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        samtools bam2fq -@ {threads} {input.wolbachia_bam} | seqkit seq -m 3000 | gzip > {output.wolbachia_fastq}
        '''

rule wri_assembly:
    input:
        wolbachia_fastq = 'data/basecalled/{sample}.wRiM23.fastq.gz'
    output:
        wolbachia_assembly = 'data/flye/{sample}/wRi/assembly.fasta'
    params:
        wolbachia_dir = 'data/flye/{sample}/wRi/'
    resources: 
        mem_mb=50000,
        runtime=180  # Shorter runtime for smaller genome
    threads: 16  
    shell:
        ''' 
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        mkdir -p {params.wolbachia_dir}
        flye --nano-hq {input.wolbachia_fastq} -t {threads} --out-dir {params.wolbachia_dir} --genome-size 1.3m
        '''

rule wri_prepare_short_reads:
    input:
        wri_r1 = 'data/short_reads/wRi_R1.fastq.gz',
        wri_r2 = 'data/short_reads/wRi_R2.fastq.gz'
    output:
        wri_trimmed = 'data/short_reads/{sample}.wri.trimmed.filtered.fastq.gz'
    threads: 4
    resources:
        mem_mb=10000,
        runtime=30
    shell:
        """
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly

        # Create output directory
        mkdir -p data/short_reads/
        # Simply combine R1 and R2 for wRi without trimming
        zcat {input.wri_r1} {input.wri_r2} | gzip > {output.wri_trimmed}
        """

rule wri_polish:
    input:
        wolbachia_assembly = 'data/flye/{sample}/wRi/assembly.fasta',
        wolbachia_short_reads = 'data/short_reads/{sample}.wri.trimmed.filtered.fastq.gz'
    output:    
        wolbachia_polished = 'data/polished/{sample}_wRi_M23.assembly.fasta'
    params:
        wolbachia_dir = 'data/polished/{sample}/wRi/'
    resources: 
        mem_mb=50000,
        runtime=180  # Shorter runtime for smaller genome
    threads: 16
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        # Create output directory
        mkdir -p {params.wolbachia_dir}
        
        # Index assembly for BWA
        bwa index {input.wolbachia_assembly}
        
        # Align Wolbachia reads to Wolbachia assembly
        bwa mem -t {threads} {input.wolbachia_assembly} {input.wolbachia_short_reads} | \
            sambamba view -f bam -S /dev/stdin | \
            sambamba sort -t {threads} -o {params.wolbachia_dir}/wolbachia_sorted.bam /dev/stdin

        # Index the sorted BAM
        sambamba index -t {threads} {params.wolbachia_dir}/wolbachia_sorted.bam

        # Mark duplicates for Wolbachia
        sambamba markdup -t {threads} {params.wolbachia_dir}/wolbachia_sorted.bam {params.wolbachia_dir}/wolbachia_dedup.bam
        
        # Index the deduplicated BAM
        sambamba index -t {threads} {params.wolbachia_dir}/wolbachia_dedup.bam
        
        export _JAVA_OPTIONS="-Xmx32g"

        # Polish Wolbachia assembly with Pilon
        pilon --genome {input.wolbachia_assembly} \
              --bam {params.wolbachia_dir}/wolbachia_dedup.bam \
              --output {params.wolbachia_dir}/wolbachia_polished \
              --threads {threads}
              
        # Copy polished assembly to final output location
        cp {params.wolbachia_dir}/wolbachia_polished.fasta {output.wolbachia_polished}
        '''

# ==============================================================================
# OPTIONAL: QUALITY ASSESSMENT
# ==============================================================================

rule wri_busco:
    input:
        wolbachia_assembly = 'data/flye/{sample}/wRi/assembly.fasta',
        wolbachia_polished = 'data/polished/{sample}_wRi_M23.assembly.fasta'
    output:
        wolbachia_assembly = '/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/busco/{sample}/wRi/assembly/short_summary.specific.rickettsiales_odb10.txt',
        wolbachia_polished = '/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/busco/{sample}/wRi/polished/short_summary.specific.rickettsiales_odb10.txt'
    params:
        wolbachia_assembly_dir = '/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/busco/{sample}/wRi/assembly/',
        wolbachia_polished_dir = '/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/busco/{sample}/wRi/polished/'
    resources: 
        mem_mb=25000,
        runtime=120
    threads: 16  
    shell:
        ''' 
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate busco

        mkdir -p {params.wolbachia_assembly_dir}
        mkdir -p {params.wolbachia_polished_dir}

        busco -i {input.wolbachia_assembly} -o {params.wolbachia_assembly_dir} -l rickettsiales_odb10 -m genome --cpu {threads}
        busco -i {input.wolbachia_polished} -o {params.wolbachia_polished_dir} -l rickettsiales_odb10 -m genome --cpu {threads}
        '''