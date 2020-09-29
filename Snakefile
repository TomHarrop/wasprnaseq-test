#!/usr/bin/env python3


#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    return {
        'l2r1': pep.get_sample(wildcards.sample).l2r1,
        'l2r2': pep.get_sample(wildcards.sample).l2r2,
        'l3r1': pep.get_sample(wildcards.sample).l3r1,
        'l3r2': pep.get_sample(wildcards.sample).l3r2}


###########
# GLOBALS #
###########

# samples
pepfile: 'data/config.yaml'

# references
ref = 'data/ref/Vvulg_24_August_2020.assembly.annotated.fna'
gff = 'data/ref/Vvulg_23_June_2020.gff3'
mrna = 'output/000_ref/vvulg.mrna.fa'


# containers
bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.86'
gffread = 'shub://TomHarrop/assembly-utils:gffread_0.12.3'
salmon = 'docker://combinelab/salmon:1.3.0'
salmontools = 'shub://TomHarrop/align-utils:salmontools_23eac84'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'

#########
# RULES #
#########

rule all:
    input:
        expand('output/020_salmon/{sample}/quant.sf',
               sample=pep.sample_table['sample_name']),
        # 'output/000_ref/gentrome.fa'


# quantify
rule salmon:
    input:
        'output/005_index/seq.bin',
        'output/005_index/pos.bin',
        r1 = 'output/010_process/{sample}.r1.fastq',
        r2 = 'output/010_process/{sample}.r2.fastq'
    output:
        'output/020_salmon/{sample}/quant.sf'
    params:
        index = 'output/005_index',
        outdir = 'output/020_salmon/{sample}'
    log:
        'output/logs/salmon.{sample}.log'
    threads:
        workflow.cores
    singularity:
        salmon
    shell:
        'salmon quant '
        '--libType ISR '
        '--index {params.index} '
        '--mates1 {input.r1} '
        '--mates2 {input.r2} '
        '--output {params.outdir} '
        '--threads {threads} '
        '--validateMappings '
        '--gcBias '
        '&> {log}'


# process the reads
rule trim:
    input:
        'output/010_process/{sample}.repair.fastq'
    output:
        r1 = 'output/010_process/{sample}.r1.fastq',
        r2 = 'output/010_process/{sample}.r2.fastq'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/trim.{sample}.log'
    threads:
        1
    container:
        bbmap
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'


rule check_pairing:
    input:
        r1 = 'output/010_process/{sample}.joined.r1.fastq',
        r2 = 'output/010_process/{sample}.joined.r2.fastq',
    output:
        pipe = pipe('output/010_process/{sample}.repair.fastq')
        # pipe = 'output/010_process/{sample}.repair.fastq'
    log:
        'output/logs/{sample}_repair.txt'
    threads:
        1
    container:
        bbmap
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '>> {output.pipe} '
        '2> {log}'

rule join_reads:
    input:
        unpack(get_reads)
    output:
        r1 = pipe('output/010_process/{sample}.joined.r1.fastq'),
        r2 = pipe('output/010_process/{sample}.joined.r2.fastq'),
    shell:
        'zcat {input.l2r1} {input.l3r1} >> {output.r1} & '
        'zcat {input.l2r2} {input.l3r2} >> {output.r2} & '
        'wait'


# process the reference
rule generate_index:
    input:
        transcriptome = 'output/000_ref/gentrome.fa',
        decoys = 'output/000_ref/decoys.txt'
    output:
        'output/005_index/seq.bin',
        'output/005_index/pos.bin'
    params:
        outdir = 'output/005_index'
    log:
        'output/logs/generate_index.log'
    threads:
        workflow.cores
    singularity:
        salmon
    shell:
        'salmon index '
        '--transcripts {input.transcriptome} '
        '--index {params.outdir} '
        '--threads {threads} '
        '--decoys {input.decoys} '
        '&> {log}'


rule generate_decoy_trancriptome:
    input:
        fasta = ref,
        transcriptome = mrna,
        gff = gff
    output:
        'output/000_ref/gentrome.fa',
        'output/000_ref/decoys.txt'
    params:
        outdir = 'output/000_ref'
    log:
        'output/logs/generate_decoy_trancriptome.log'
    threads:
        workflow.cores
    singularity:
        salmontools
    shell:
        'generateDecoyTranscriptome.sh '
        '-j {threads} '
        '-b /usr/bin/bedtools '
        '-m /usr/local/bin/mashmap '
        '-a {input.gff} '
        '-g {input.fasta} '
        '-t {input.transcriptome} '
        '-o {params.outdir} '
        '&> {log}'

rule gffread:
    input:
        ref = ref,
        gff = gff,
        fai = f'{ref}.fai'
    output:
        mrna = mrna
    log:
        'output/logs/gffread.log'
    container:
        gffread
    shell:
        'gffread '
        '{input.gff} '
        '-w {output.mrna} '
        '-g {input.ref} '
        '&> {log}'


rule faidx:
    input:
        '{path}/{file}.{ext}'
    output:
        '{path}/{file}.{ext}.fai'
    wildcard_constraints:
        ext = 'fasta|fa|fna'
    singularity:
        samtools
    shell:
        'samtools faidx {input}'