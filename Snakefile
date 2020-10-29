#!/usr/bin/env python3


#############
# FUNCTIONS #
#############

def get_fastqc_reads(wildcards):
    if wildcards.type == 'raw':
        return {'r1': 'output/010_process/{sample}.joined.r1.fastq',
                'r2': 'output/010_process/{sample}.joined.r2.fastq'}
    elif wildcards.type == 'processed':
        return {'r1': 'output/010_process/{sample}.r1.fastq',
                'r2': 'output/010_process/{sample}.r2.fastq'}


def get_reads(wildcards):
    input_keys = ['l2r1', 'l2r2', 'l3r1', 'l3r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}


###########
# GLOBALS #
###########

# samples
pepfile: 'data/config.yaml'
all_samples = pep.sample_table['sample_name']

# references
ref = 'data/ref/GCA_014466185.1_ASM1446618v1_genomic.fna'
gff = 'data/ref/GCA_014466185.1_ASM1446618v1_genomic.gff'
mrna = 'output/000_ref/vvulg.mrna.fa'
# mrna = 'data/ref/Vespula_vulgaris.transcripts.fa'

pipelines = ['salmon', 'star']

# containers
bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.86'
bioconductor = ('shub://TomHarrop/r-containers:bioconductor_3.11'
                '@ae3e49fbdb6c7a9a05fc5b88cc55ac3663b40036')    # has tximeta
fastqc = 'docker://biocontainers/fastqc:v0.11.9_cv7'
gffread = 'shub://TomHarrop/assembly-utils:gffread_0.12.3'
multiqc = 'docker://ewels/multiqc:1.9'
salmon = 'docker://combinelab/salmon:1.3.0'
salmontools = 'shub://TomHarrop/align-utils:salmontools_23eac84'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
star = 'shub://TomHarrop/align-utils:star_2.7.6a'


#########
# RULES #
#########

wildcard_constraints:
    sample = '|'.join(all_samples),
    pl = '|'.join(pipelines)

rule all:
    input:
        expand('output/030_deseq/wald/res.annot.{pl}.csv',
               pl=pipelines),
        # 'output/017_multiqc/multiqc_report.html'

# DE analysis
rule de_wald:
    input:
        dds = 'output/030_deseq/dds.{pl}.Rds',
    params:
        alpha = 0.1,
        lfc_threshold = 0.585   # log(1.5, 2)
    output:
        ma = 'output/030_deseq/wald/ma.{pl}.pdf',
        res = 'output/030_deseq/wald/res.{pl}.csv'
    log:
        'output/logs/de_wald.{pl}.log'
    threads:
        min(16, workflow.cores)
    container:
        bioconductor
    script:
        'src/de_wald.R'

def pick_quant_files(wildcards):
    if wildcards.pl == 'salmon':
        my_files = expand('output/020_salmon/{sample}/quant.sf',
                          sample=all_samples)
    elif wildcards.pl == 'star':
        my_files = expand('output/025_star/pass2/{sample}.ReadsPerGene.out.tab',
                          sample=all_samples)
    else:
        raise ValueError(f'wtf {wildcards.pl}')
    return(my_files)

rule generate_deseq_object:
    input:
        # quant_files = expand('output/020_salmon/{sample}/quant.sf',
        #                      sample=all_samples),
        quant_file = pick_quant_files,
        gff = gff,
        mrna = mrna
    output:
        dds = 'output/030_deseq/dds.{pl}.Rds'
    params:
        index = 'output/005_index',
        # script = 'src/generate_deseq_object.{pl}.R'
        script = lambda wildcards: print(wildcards.pl)
    log:
        'output/logs/generate_deseq_object.{pl}.log'
    singularity:
        bioconductor
    script:
        '{params.script}'

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
        r1 = temp('output/010_process/{sample}.joined.r1.fastq'),
        r2 = temp('output/010_process/{sample}.joined.r2.fastq'),
    shell:
        'zcat {input.l2r1} {input.l3r1} >> {output.r1} & '
        'zcat {input.l2r2} {input.l3r2} >> {output.r2} & '
        'wait'

# generic annotation rule
rule annot_res:
    input:
        res = '{path}/{file}.{pl}.csv',
        annot = 'output/000_ref/annot.csv'
    output:
        res_annot = '{path}/{file}.annot.{pl}.csv'
    log:
        'output/logs/annot_res.{path}.{file}.{pl}.log'
    container:
        bioconductor
    script:
        'src/annot_res.R'

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


rule generate_gentrome:
    input:
        fasta = ref,
        transcriptome = mrna
    output:
        'output/000_ref/gentrome.fa',
    container:
        salmon
    shell:
        'cat {input.transcriptome} {input.fasta} > {output}'

rule generate_decoys:
    input:
        f'{ref}.fai'
    output:
        temp('output/000_ref/decoys.txt')
    container:
        salmon
    shell:
        'cut -f1 {input} > {output}'

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

rule parse_annotations:
    input:
        gff = gff
    output:
        annot = 'output/000_ref/annot.csv'
    log:
        'output/logs/parse_annotations.log'
    container:
        bioconductor
    script:
        'src/parse_annotations.R'

# fastqc
rule multiqc:
    input:
        expand('output/015_fastqc/{type}/{sample}.fastqc',
               type=[
                    'raw',
                    'processed'
                    ],
               sample=all_samples),
        expand('output/020_salmon/{sample}/quant.sf',
               sample=all_samples)
    output:
        'output/017_multiqc/multiqc_report.html'
    params:
        outdir = 'output/017_multiqc'
    log:
        'output/logs/multiqc.log'
    container:
        multiqc
    shell:
        'multiqc '
        '-o {params.outdir} '
        'output '
        '2> {log}'

rule fastqc:
    input:
        unpack(get_fastqc_reads)
    output:
        'output/015_fastqc/{type}/{sample}.fastqc'
    params:
        outdir = 'output/015_fastqc/{type}'
    log:
        'output/logs/fastqc.{sample}.{type}.log'
    threads:
        2
    container:
        fastqc
    shell:
        'fastqc '
        '--threads {threads} '
        '-o {params.outdir} '
        '{input.r1} {input.r2} '
        '&> {log} '
        '; touch {output}'


# test star mapping
rule star_target:
    input:
        expand('output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
               sample=all_samples)

rule star_second_pass:
    input:
        r1 = 'output/010_process/{sample}.r1.fastq',
        r2 = 'output/010_process/{sample}.r2.fastq',
        star_reference = 'output/007_star-index/SA',
        junctions = expand('output/025_star/pass1/{sample}.SJ.out.tab',
                           sample=all_samples)
    output:
        bam = 'output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
        counts = 'output/025_star/pass2/{sample}.ReadsPerGene.out.tab'
    threads:
        workflow.cores
    params:
        genome_dir = 'output/007_star-index',
        prefix = 'output/025_star/pass2/{sample}.'
    log:
        'output/logs/star_second_pass.{sample}.log'
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outReadsUnmapped Fastx '
        '--quantMode GeneCounts '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_first_pass:
    input:
        r1 = 'output/010_process/{sample}.r1.fastq',
        r2 = 'output/010_process/{sample}.r2.fastq',
        star_reference = 'output/007_star-index/SA'
    output:
        sjdb = 'output/025_star/pass1/{sample}.SJ.out.tab'
    threads:
        workflow.cores
    params:
        genome_dir = 'output/007_star-index',
        prefix = 'output/025_star/pass1/{sample}.'
    log:
        'output/logs/star_first_pass.{sample}.log'
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '          # troubleshoot gtf
        # '--outSAMtype SAM '               # troubleshoot gtf
        # '--quantMode GeneCounts '       # troubleshoot gtf
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_index:
    input:
        fasta = ref,
        gff = gff
    output:
        'output/007_star-index/SA'
    params:
        outdir = 'output/007_star-index'
    log:
        'output/logs/star_index.log'
    threads:
        workflow.cores
    container:
        star
    shell:
        'STAR '
        'runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gff} '
        '--genomeSAindexNbases 12 '
        '--sjdbGTFtagExonParentTranscript Parent '
        '--sjdbGTFtagExonParentGene locus_tag '
        # '--sjdbGTFtagExonParentGeneName Name '
        '&> {log}'

