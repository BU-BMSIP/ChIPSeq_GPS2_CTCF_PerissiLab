# Snakefile with conda envs

# GLOBAL config
WORKDIR = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/iMargi"
SRA_ID = "SRR8206679"

rule all:
    input:
        expand(WORKDIR + "/fastqc/{sample}_fastqc.html", sample=[f"{SRA_ID}_1", f"{SRA_ID}_2"]),
        WORKDIR + "/multiqc/multiqc_report.html",
        WORKDIR + "/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        WORKDIR + "/output/final_test_sample.pairs.gz",
        WORKDIR + "/output/pipelineStats_test_sample.log"

# download fastq.gz
rule download_fastq:
    output:
        fq1 = WORKDIR + "/fastq/SRR8206679_1.fastq.gz",
        fq2 = WORKDIR + "/fastq/SRR8206679_2.fastq.gz"
    shell:
        """
        mkdir -p {WORKDIR}/fastq
        wget -O {output.fq1} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR820/009/SRR8206679/SRR8206679_1.fastq.gz
        wget -O {output.fq2} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR820/009/SRR8206679/SRR8206679_2.fastq.gz
        """

# FASTQC with your env
rule fastqc:
    input:
        fq=WORKDIR + "/fastq/{sample}.fastq.gz"
    output:
        html=WORKDIR + "/fastqc/{sample}_fastqc.html",
        zip=WORKDIR + "/fastqc/{sample}_fastqc.zip"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/fastqc_env.yml"
    threads: 8
    shell:
        """
        mkdir -p {WORKDIR}/fastqc
        fastqc {input.fq} -o {WORKDIR}/fastqc
        """

# MULTIQC with your env
rule multiqc:
    input:
        expand(WORKDIR + "/fastqc/{sample}_fastqc.zip", sample=[f"{SRA_ID}_1", f"{SRA_ID}_2"])
    output:
        html=WORKDIR + "/multiqc/multiqc_report.html"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/multiqc_env.yml"
    threads: 1
    shell:
        """
        mkdir -p {WORKDIR}/multiqc
        multiqc {WORKDIR}/fastqc -o {WORKDIR}/multiqc
        """

#download reference
rule download_reference:
    output:
        ref = WORKDIR + "/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    shell:
        """
        mkdir -p {WORKDIR}/ref
        wget -O {output.ref}.gz \
          https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
        gunzip -c {output.ref}.gz > {output.ref}
        rm {output.ref}.gz
        """

# run iMARGI with singularity
rule run_imargi:
    input:
        r1 = WORKDIR + "/fastq/SRR8206679_1.fastq.gz",
        r2 = WORKDIR + "/fastq/SRR8206679_2.fastq.gz",
        ref = WORKDIR + "/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    output:
        pairs = WORKDIR + "/output/final_test_sample.pairs.gz",
        log = WORKDIR + "/output/pipelineStats_test_sample.log"
    singularity:
        "docker://zhonglab/imargi"
    threads: 16
    shell:
        """
        mkdir -p {WORKDIR}/output
        imargi_wrapper.sh \
        -r hg38 \
        -N test_sample \
        -t {threads} \
        -i {WORKDIR}/ref/bwa_index/bwa_index_hg38 \
        -g {WORKDIR}/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
        -1 {WORKDIR}/fastq/SRR8206679_1.fastq.gz \
        -2 {WORKDIR}/fastq/SRR8206679_2.fastq.gz \
        -o {WORKDIR}/output
        """