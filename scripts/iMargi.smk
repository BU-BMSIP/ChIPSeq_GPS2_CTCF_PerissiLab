# GLOBAL config
WORKDIR = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/iMargi"
SRA_ID = ["SRR9900120", "SRR9900121", "SRR9900122"]
SRA_PATH_MAP = {
    "SRR9900120": "SRR990/000",
    "SRR9900121": "SRR990/001",
    "SRR9900122": "SRR990/002"
}

# rule all
rule all:
    input:
        expand(WORKDIR + "/fastqc/{sample}_fastqc.html", sample=[f"{sra}_1" for sra in SRA_ID] + [f"{sra}_2" for sra in SRA_ID]),
        WORKDIR + "/multiqc/multiqc_report.html",
        WORKDIR + "/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        expand(WORKDIR + "/output/final_{sra}.pairs.gz", sra=SRA_ID),
        expand(WORKDIR + "/output/pipelineStats_{sra}.log", sra=SRA_ID)

# download fastq
rule download_fastq:
    output:
        fq1 = WORKDIR + "/fastq/{sra}_1.fastq.gz",
        fq2 = WORKDIR + "/fastq/{sra}_2.fastq.gz"
    params:
        path = lambda wildcards: SRA_PATH_MAP[wildcards.sra]
    shell:
        """
        mkdir -p {WORKDIR}/fastq
        wget -nc -O {output.fq1} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{params.path}/{wildcards.sra}/{wildcards.sra}_1.fastq.gz
        wget -nc -O {output.fq2} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{params.path}/{wildcards.sra}/{wildcards.sra}_2.fastq.gz
        """

# FASTQC
rule fastqc:
    input:
        fq = WORKDIR + "/fastq/{sample}.fastq.gz"
    output:
        html = WORKDIR + "/fastqc/{sample}_fastqc.html",
        zip = WORKDIR + "/fastqc/{sample}_fastqc.zip"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/fastqc_env.yml"
    threads: 8
    shell:
        """
        mkdir -p {WORKDIR}/fastqc
        fastqc {input.fq} -o {WORKDIR}/fastqc
        """

# MULTIQC
rule multiqc:
    input:
        expand(WORKDIR + "/fastqc/{sample}_fastqc.zip", sample=[f"{sra}_1" for sra in SRA_ID] + [f"{sra}_2" for sra in SRA_ID])
    output:
        html = WORKDIR + "/multiqc/multiqc_report.html"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/multiqc_env.yml"
    threads: 1
    shell:
        """
        mkdir -p {WORKDIR}/multiqc
        multiqc {WORKDIR}/fastqc -o {WORKDIR}/multiqc
        """

# download reference
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

# bwa index
rule bwa_index:
    input:
        fasta = WORKDIR + "/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    output:
        expand(WORKDIR + "/ref/bwa_index/{prefix}.{ext}",
            prefix="bwa_index_hg38",
            ext=["amb","ann","bwt","pac","sa"]
        )
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bwa_env.yml"
    threads: 4
    shell:
        """
        mkdir -p {WORKDIR}/ref/bwa_index
        bwa index -p {WORKDIR}/ref/bwa_index/bwa_index_hg38 {input.fasta}
        """

# run iMARGI
rule run_imargi:
    input:
        r1 = WORKDIR + "/fastq/{sra}_1.fastq.gz",
        r2 = WORKDIR + "/fastq/{sra}_2.fastq.gz",
        fasta = WORKDIR + "/ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        bwa_index = expand(WORKDIR + "/ref/bwa_index/bwa_index_hg38.{ext}", 
                           ext=["amb","ann","bwt","pac","sa"])
    output:
        pairs = WORKDIR + "/output/final_{sra}.pairs.gz",
        log = WORKDIR + "/output/pipelineStats_{sra}.log"
    singularity:
        "docker://zhonglab/imargi"
    threads: 16
    shell:
        """
        mkdir -p {WORKDIR}/output
        iMargi/imargi_wrapper.sh \
            -r hg38 \
            -N {wildcards.sra} \
            -t {threads} \
            -g {input.fasta} \
            -i {WORKDIR}/ref/bwa_index/bwa_index_hg38 \
            -c {WORKDIR}/ref/chromsize.hg38.txt \
            -R {WORKDIR}/ref/AluI_frags.bed.gz \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {WORKDIR}/output
        """
