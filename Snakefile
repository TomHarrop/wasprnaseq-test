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

pepfile: 'data/config.yaml'

bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.86'


#########
# RULES #
#########

rule all:
    input:
        expand('output/010_process/{sample}.r1.fastq',
               sample=pep.sample_table['sample_name'])



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