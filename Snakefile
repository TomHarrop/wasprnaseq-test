pepfile: 'data/config.yaml'

def get_reads(wildcards):
    return {
        'r1': pep.get_sample(wildcards.sample).read1,
        'r2': pep.get_sample(wildcards.sample).read2}


rule all:
    input:
        expand('output/010_process/{sample}.fastq.gz',
               sample=pep.sample_table['sample_name'])


# WORKING, STILL NEED TO MAP THE MULTPILE LANES BACK TO SAMPLE
rule a:
    input:
        unpack(get_reads)
    output:
        'output/010_process/{sample}.fastq.gz'
    shell:
        'bbmap stuff '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output}'