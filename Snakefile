#Completed on 08/28/2025
#to run this script:
#cd /private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_Dsim_assembly/
#conda activate snakemake
# snakemake --executor slurm --default-resources slurm_partition=medium runtime=720 mem_mb=100000 -j 10 -s Snakefile

#Global Variables:

import os  
import glob

samples = ['20250828_Nanopore']
conda: '/private/groups/russelllab/jodie/bootcamp2024/scripts/read_filtering.yaml'   

rule all:
    input:
        expand('data/integrated/{sample}_Dsim_integrated.assembly.fasta', sample=samples),
        expand('data/comparison/{sample}_stats_summary.txt', sample=samples),
        expand('data/comparison/{sample}_vs_ref_dnadiff.txt', sample=samples)

# ==============================================================================
# SHARED PREPROCESSING STEPS
# ==============================================================================

rule basecalling:
    input:
        pod5 = 'data/pod5'
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
        slurm_extra='--gpus-per-node=4'
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

        samtools view -@ {threads} -b -F 4 {input.bam} > {output.host_bam}
        samtools view -@ {threads} -b -f 4 {input.bam} > {output.wolbachia_bam}
        '''

# ==============================================================================
# D. SIMULANS HOST GENOME PROCESSING  
# ==============================================================================

rule host_bam2fastq:
    input:
        host_bam = 'data/aligned/{sample}.DsimM23.bam'
    output:
        host_fastq = 'data/basecalled/{sample}.DsimM23.fastq.gz'
    resources: 
        mem_mb=50000,
        runtime=120
    threads: 8  
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        samtools bam2fq -@ {threads} {input.host_bam} | seqkit seq -m 1000 | gzip > {output.host_fastq}
        '''

rule host_flye_assembly:
    input:
        host_fastq = 'data/basecalled/{sample}.DsimM23.fastq.gz'
    output:
        host_assembly = 'data/flye/{sample}/Dsim/assembly.fasta'
    params:
        host_dir = 'data/flye/{sample}/Dsim/'
    resources: 
        slurm_partition='long',
        mem_mb=200000,  # 200GB for large genome
        runtime=1440    # 24 hours for D. simulans assembly
    threads: 24  
    shell:
        ''' 
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        mkdir -p {params.host_dir}
        flye --nano-hq {input.host_fastq} -t {threads} --out-dir {params.host_dir} --genome-size 140m
        '''
rule host_medaka_polish:
    input:
        assembly="data/flye/{sample}/Dsim/assembly.fasta",
        reads="data/basecalled/{sample}.DsimM23.fastq.gz"
    output:
        "data/polished/{sample}_Dsim_flye_polished.fasta"
    threads: 20
    resources:
        mem_mb=150000,
        runtime=720,
        slurm_partition="medium"
    shell:
        """
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        mkdir -p data/polished/{wildcards.sample}/flye_polish/
        
        # Map reads to assembly
        minimap2 -ax map-ont -t {threads} {input.assembly} {input.reads} | \
            samtools sort -@ {threads} > data/polished/{wildcards.sample}/flye_polish/mapped.bam
        samtools index data/polished/{wildcards.sample}/flye_polish/mapped.bam
        
        # Step 1: Generate features from BAM
        medaka inference data/polished/{wildcards.sample}/flye_polish/mapped.bam \
            data/polished/{wildcards.sample}/flye_polish/consensus.hdf \
            --model r1041_e82_400bps_hac_v4.3.0 --threads {threads}
        
        # Step 2: Generate consensus FASTA from features
        medaka consensus_from_features data/polished/{wildcards.sample}/flye_polish/consensus.hdf \
            {output}
        """

rule reference_based_racon_polish:
    input:
        reads="data/basecalled/{sample}.DsimM23.fastq.gz"
    output:
        "data/polished/{sample}_Dsim_ref_polished.fasta"
    threads: 20
    resources:
        mem_mb=150000,
        runtime=720,
        slurm_partition="medium"
    shell:
        """
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        mkdir -p data/polished/{wildcards.sample}/ref_polish/
        
        # Decompress reference if needed
        if [[ data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz == *.gz ]]; then
            gunzip -c data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz > data/polished/{wildcards.sample}/ref_polish/reference.fasta
        else
            cp data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz data/polished/{wildcards.sample}/ref_polish/reference.fasta
        fi
        
        # Map nanopore reads to reference
        minimap2 -ax map-ont -t {threads} data/polished/{wildcards.sample}/ref_polish/reference.fasta {input.reads} | \
            samtools sort -@ {threads} > data/polished/{wildcards.sample}/ref_polish/mapped.bam
        samtools index data/polished/{wildcards.sample}/ref_polish/mapped.bam
        
        # Step 1: Generate features from BAM
        medaka inference data/polished/{wildcards.sample}/ref_polish/mapped.bam \
            data/polished/{wildcards.sample}/ref_polish/consensus.hdf \
            --model r1041_e82_400bps_hac_v4.3.0 --threads {threads}
        
        # Step 2: Generate consensus FASTA from features
        medaka consensus_from_features data/polished/{wildcards.sample}/ref_polish/consensus.hdf \
            {output}
        """
        
rule reference_scaffold:
    input:
        assembly="data/flye/{sample}/Dsim/assembly.fasta",
        reference="data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz"
    output:
        "data/scaffolded/{sample}_Dsim_scaffolded.fasta"
    shell:
        """
        ragtag.py scaffold {input.reference} {input.assembly} -o data/scaffolded/{wildcards.sample}_Dsim_ragtag/
        cp data/scaffolded/{wildcards.sample}_Dsim_ragtag/ragtag.scaffold.fasta {output}
        """

rule integrate_assemblies:
    input:
        flye_polished = 'data/polished/{sample}_Dsim_flye_polished.fasta',
        ref_polished = 'data/polished/{sample}_Dsim_ref_polished.fasta'
    output:
        integrated = 'data/integrated/{sample}_Dsim_integrated.assembly.fasta'
    params:
        work_dir = 'data/integrated/{sample}/'
    resources: 
        mem_mb=100000,
        runtime=360  # 6 hours
    threads: 16
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        mkdir -p {params.work_dir}
        
        # Use Quickmerge to integrate assemblies
        # Reference-polished as "reference", Flye as "query" 
        quickmerge -d {input.ref_polished} -q {input.flye_polished} \
            -hco 5.0 -c 1.5 -l 100000 -ml 5000 \
            -pre {params.work_dir}/integrated
            
        # Copy final assembly
        cp {params.work_dir}/integrated.fasta {output.integrated}
        '''

# ==============================================================================
# QUALITY ASSESSMENT
# ==============================================================================

rule assembly_stats:
    input:
        flye_assembly = 'data/flye/{sample}/Dsim/assembly.fasta',
        flye_polished = 'data/polished/{sample}_Dsim_flye_polished.fasta',
        ref_polished = 'data/polished/{sample}_Dsim_ref_polished.fasta',
        integrated = 'data/integrated/{sample}_Dsim_integrated.assembly.fasta'
    output:
        stats = 'data/comparison/{sample}_stats_summary.txt'
    resources: 
        mem_mb=25000,
        runtime=60
    threads: 4
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        echo "Assembly Statistics Comparison" > {output.stats}
        echo "===============================" >> {output.stats}
        echo "" >> {output.stats}
        
        echo "Flye De Novo Assembly:" >> {output.stats}
        assembly-stats {input.flye_assembly} >> {output.stats}
        echo "" >> {output.stats}
        
        echo "Flye Polished Assembly:" >> {output.stats}
        assembly-stats {input.flye_polished} >> {output.stats}
        echo "" >> {output.stats}
        
        echo "Reference Polished Assembly:" >> {output.stats}
        assembly-stats {input.ref_polished} >> {output.stats}
        echo "" >> {output.stats}
        
        echo "Integrated Assembly:" >> {output.stats}
        assembly-stats {input.integrated} >> {output.stats}
        '''

rule busco_assessment:
    input:
        flye_polished = 'data/polished/{sample}_Dsim_flye_polished.fasta',
        ref_polished = 'data/polished/{sample}_Dsim_ref_polished.fasta',
        integrated = 'data/integrated/{sample}_Dsim_integrated.assembly.fasta'
    output:
        flye_busco = 'busco/{sample}/flye_polished/short_summary.specific.diptera_odb10.txt',
        ref_busco = 'busco/{sample}/ref_polished/short_summary.specific.diptera_odb10.txt', 
        integrated_busco = 'busco/{sample}/integrated/short_summary.specific.diptera_odb10.txt'
    params:
        flye_dir = 'busco/{sample}/flye_polished/',
        ref_dir = 'busco/{sample}/ref_polished/',
        integrated_dir = 'busco/{sample}/integrated/'
    resources: 
        mem_mb=50000,
        runtime=300  # 5 hours for larger genome
    threads: 16  
    shell:
        ''' 
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate busco

        mkdir -p {params.flye_dir}
        mkdir -p {params.ref_dir}
        mkdir -p {params.integrated_dir}

        # BUSCO for Flye polished assembly
        busco -i {input.flye_polished} -o {params.flye_dir} -l diptera_odb10 -m genome --cpu {threads}
        
        # BUSCO for reference polished assembly  
        busco -i {input.ref_polished} -o {params.ref_dir} -l diptera_odb10 -m genome --cpu {threads}
        
        # BUSCO for integrated assembly
        busco -i {input.integrated} -o {params.integrated_dir} -l diptera_odb10 -m genome --cpu {threads}
        '''

# ==============================================================================
# OPTIONAL: HYBRID POLISHING WITH SHORT READS (if available)
# ==============================================================================

rule host_prepare_short_reads:
    input:
        host_r1 = 'data/short_reads/Dsim_R1.fastq.gz',
        host_r2 = 'data/short_reads/Dsim_R2.fastq.gz'
    output:
        host_trimmed_r1 = 'data/short_reads/{sample}.dsim.trimmed.R1.fastq.gz',
        host_trimmed_r2 = 'data/short_reads/{sample}.dsim.trimmed.R2.fastq.gz'
    threads: 8
    resources:
        mem_mb=25000,
        runtime=120
    shell:
        """
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly

        # Create output directory
        mkdir -p data/short_reads/
        
        # Trim and filter short reads with fastp
        fastp -i {input.host_r1} -I {input.host_r2} \
              -o {output.host_trimmed_r1} -O {output.host_trimmed_r2} \
              --thread {threads} --length_required 50 --cut_tail
        """

rule hybrid_polish_integrated:
    input:
        integrated = 'data/integrated/{sample}_Dsim_integrated.assembly.fasta',
        host_r1 = 'data/short_reads/{sample}.dsim.trimmed.R1.fastq.gz',
        host_r2 = 'data/short_reads/{sample}.dsim.trimmed.R2.fastq.gz'
    output:    
        final_polished = 'data/final/{sample}_Dsim_final.assembly.fasta'
    params:
        work_dir = 'data/final/{sample}/'
    resources: 
        mem_mb=150000,
        runtime=480  # 8 hours
    threads: 20
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        # Create output directory
        mkdir -p {params.work_dir}
        
        # Index assembly for BWA
        bwa index {input.integrated}
        
        # Align short reads to integrated assembly
        bwa mem -t {threads} {input.integrated} {input.host_r1} {input.host_r2} | \
            sambamba view -f bam -S /dev/stdin | \
            sambamba sort -t {threads} -o {params.work_dir}/host_sorted.bam /dev/stdin

        # Index the sorted BAM
        sambamba index -t {threads} {params.work_dir}/host_sorted.bam

        # Mark duplicates
        sambamba markdup -t {threads} {params.work_dir}/host_sorted.bam {params.work_dir}/host_dedup.bam
        
        # Index the deduplicated BAM
        sambamba index -t {threads} {params.work_dir}/host_dedup.bam
        
        export _JAVA_OPTIONS="-Xmx128g"

        # Final polish with Pilon
        pilon --genome {input.integrated} \
              --bam {params.work_dir}/host_dedup.bam \
              --output {params.work_dir}/final_polished \
              --threads {threads} --fix snps,indels,gaps
              
        # Copy final assembly 
        cp {params.work_dir}/final_polished.fasta {output.final_polished}
        '''

# ==============================================================================
# ASSEMBLY COMPARISON AND VALIDATION
# ==============================================================================

rule compare_to_reference:
    input:
        integrated = 'data/integrated/{sample}_Dsim_integrated.assembly.fasta',
        reference = 'data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz'
    output:
        comparison = 'data/comparison/{sample}_vs_ref_dnadiff.txt'
    params:
        work_dir = 'data/comparison/{sample}/'
    resources: 
        mem_mb=100000,
        runtime=240  # 4 hours
    threads: 16
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        mkdir -p {params.work_dir}
        
        # Decompress reference if needed
        if [[ {input.reference} == *.gz ]]; then
            gunzip -c {input.reference} > {params.work_dir}/reference.fasta
        else
            cp {input.reference} {params.work_dir}/reference.fasta
        fi
        
        # Compare assemblies with dnadiff
        dnadiff -p {params.work_dir}/comparison {params.work_dir}/reference.fasta {input.integrated}
        
        # Generate summary report
        echo "Assembly vs Reference Comparison" > {output.comparison}
        echo "================================" >> {output.comparison}
        cat {params.work_dir}/comparison.report >> {output.comparison}
        '''

# ==============================================================================
# ALTERNATIVE: RAGTAG SCAFFOLDING (if Quickmerge not available)
# ==============================================================================

rule ragtag_scaffold:
    input:
        flye_polished = 'data/polished/{sample}_Dsim_flye_polished.fasta',
        reference = 'data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz'
    output:
        scaffolded = 'data/ragtag/{sample}_Dsim_scaffolded.assembly.fasta'
    params:
        work_dir = 'data/ragtag/{sample}/',
        reference_unzipped = 'data/ragtag/{sample}/reference.fasta'
    resources: 
        mem_mb=100000,
        runtime=360  # 6 hours
    threads: 16
    shell:
        '''
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate assembly
        
        mkdir -p {params.work_dir}
        
        # Decompress reference
        if [[ {input.reference} == *.gz ]]; then
            gunzip -c {input.reference} > {params.reference_unzipped}
        else
            cp {input.reference} {params.reference_unzipped}
        fi
        
        # Scaffold with RagTag
        ragtag.py scaffold {params.reference_unzipped} {input.flye_polished} -o {params.work_dir} -t {threads}
        
        # Copy final scaffolded assembly
        cp {params.work_dir}/ragtag.scaffold.fasta {output.scaffolded}
        '''