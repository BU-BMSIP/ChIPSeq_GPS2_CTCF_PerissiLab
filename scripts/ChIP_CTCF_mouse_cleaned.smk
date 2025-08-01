wildcard_constraints:
    orig_file="[^/]+",
    file="[^/]+"

import pandas as pd
import os

sample_csv = pd.read_csv('Sample_Data_for_CTCF.csv')
updated_samples_sheet = pd.read_csv('Updated_Sample_Sheet.csv')
peak_calls_with_control = pd.read_csv("Peak_Calls_with_Control.csv")
peak_calls_without_control = pd.read_csv("Peak_Calls_without_Control.csv")

# Add a 'Control' column to the without_control DataFrame
peak_calls_without_control['Control'] = 'None'

# Combine both DataFrames
all_peak_calls = pd.concat([peak_calls_with_control, peak_calls_without_control], ignore_index=True)

EXTENSIONS = [1, 2, 3, 4, 'rev.1', 'rev.2']
DIRECTORY = ['CTCF_3T3L1']
samples_to_concat = updated_samples_sheet.query('Concat == "yes"')["Name"].tolist()

def get_input_files(wildcards, column, dataframe):
    files = dataframe.loc[dataframe['Peak_file'] == wildcards.peak_file, column].values[0]
    return [f"{wildcards.directory}/results/sorted/{file.strip()}.sorted.bam" for file in files.split(';')]

# SAMPLES_TO_COMPARE/REGIONS_FOR_COMPARISON
SAMPLES_TO_COMPARE = {
    "CTCF_all_T_samples": ["CTCF_t1", "CTCF_t2", "CTCF_t3", "CTCF_t4"],
}
CTCF_TIMEPOINTS = ["t1", "t2", "t3", "t4"]

REGIONS_OF_INTEREST = {
    "CTCF_union_peaks": "CTCF_3T3L1/results/filtered/CTCF_all_timepoints_union.bed",
    "GPS2_d6": "Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_bf.narrowPeak",
    "GPS2_siCTL_peaks": "Silencing/GPS2/results/filtered/sictl_bf.narrowPeak"
}

SIZES = [100, 200, 300, 400, 500]

rule all:
    input:
        # 1. Download raw samples
        expand(
            "{orig_dir}/raw_samples/{orig_file}.fastq.gz",
            zip,
            orig_dir=sample_csv['Target Directory'],
            orig_file=sample_csv['File name']
        ),

        # 2. Downstream FASTQ
        expand(
            "{dir}/samples/{file}.fastq.gz",
            zip,
            dir=updated_samples_sheet['Directory'],
            file=updated_samples_sheet['Name']
        ),

        # 3. QC outputs
        expand(
            "{dir}/results/fastqc/{file}_fastqc.html",
            zip,
            dir=updated_samples_sheet['Directory'],
            file=updated_samples_sheet['Name']
        ),
        expand(
            "{dir}/results/trimqc/{file}.trim_fastqc.html",
            zip,
            dir=updated_samples_sheet['Directory'],
            file=updated_samples_sheet['Name']
        ),
        expand(
            "{dir}/results/flagstat/{file}.txt",
            zip,
            dir=updated_samples_sheet['Directory'],
            file=updated_samples_sheet['Name']
        ),

        # 4. MultiQC report
        expand(
            "{directory}/multiqc/multiqc_report.html",
            directory=DIRECTORY
        ),

        # 5. Peak calls & motifs
        expand(
            "{directory}/results/filtered/{peak_file}_bf.narrowPeak",
            zip,
            directory=all_peak_calls['Directory'],
            peak_file=all_peak_calls['Peak_file']
        ),
        expand(
            "{directory}/results/motifs/{peak_file}",
            zip,
            directory=all_peak_calls['Directory'],
            peak_file=all_peak_calls['Peak_file']
        ),
        expand(
            "{directory}/results/promoter_motifs/{peak_file}",
            zip,
            directory=all_peak_calls['Directory'],
            peak_file=all_peak_calls['Peak_file']
        ),

        # 6. BigWig
        expand(
            "{dir}/results/bigwig/{file}.bw",
            zip,
            dir=updated_samples_sheet['Directory'],
            file=updated_samples_sheet['Name']
        ),

        # 7. CTCF–GPS2 intersections
        expand(
            "CTCF_3T3L1/results/intersect/CTCF_{tp}_GPS2_overlap.bed",
            tp=CTCF_TIMEPOINTS
        ),

        # 8. Signal matrices & profiles
        "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_promoters.gz",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_promoters_profile.png",
        "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_peaks.gz",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_peaks_profile.png",

        # 9. Additional analyses
        # 9a. GPS2 silencing vs CTCF promoters
        "Silencing/GPS2/results/plots/GPS2_signal_on_CTCF_t1_promoters_profile.png",

        # 9b. Adipocyte differentiation: GPS2 on CTCF t1 promoters
        "Adipocyte_differentiation/GPS2/results/matrix/gps2_day0_on_CTCF_t1_promoters.gz",
        "Adipocyte_differentiation/GPS2/results/plots/gps2_day0_on_CTCF_t1_promoters_profile.png",

        # 9c. Merged peak BEDs
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_preadipocyte_merged.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_adipocyte_merged.bed",

        # 9d. Pairwise intersections
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/GPS2_siCTL_vs_CTCF_preadipocyte.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/GPS2_d6_vs_CTCF_adipocyte.bed",

        # 9e. Union of overlaps
        "CTCF_3T3L1/results/intersect/CTCF_GPS2_overlap_union_t1to4_3T3L1_mm39.bed",

        # 9f. CommonGeneFiltered GPS2 peaks analysis
        "CTCF_3T3L1/results/matrix/CTCF_signal_over_CommonGeneFiltered_GPS2_peaks.gz",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_peaks_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_promoter_only_peaks_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_non_promoter_peaks_profile.png",

        # 9g. Motif directories
        directory('/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/motifs/GPS2_CTCF_Intersect'),
        directory('/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/motifs/Filtered_GPS2_Promoter'),
        directory('/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/motifs/Filtered_GPS2_Non_Promoter'),

        # 9h. Motif flank analyses
        expand(
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/motifs/Filtered_GPS2_flanks_{size}bp",
            size=SIZES
        ),
        expand(
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/motifs/Filtered_GPS2d6_flanks_{size}bp",
            size=SIZES
        ),
        expand(
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/motifs/Filtered_ATF4_flanks_{size}bp",
            size=SIZES
        ),

        # 9i. Motif center analyses
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/motifs/GPS2_d6_promoter_center",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/motifs/ATF4_promoter_center",

        # 9j. Extended CTCF profiles
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_500bpextended.png",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_500bpextended_heatmap.png",

        # 9k. d6 GPS2 peaks profiles
        "CTCF_3T3L1/results/plots/CTCF_signal_on_d6_GPS2_peaks_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_d6_GPS2_peaks_heatmap.png",

        # 9l. Common genes d6 profiles
        "CTCF_3T3L1/results/plots/GPS2_CTCF_signal_on_common_genes_d6_profile.png",
        "CTCF_3T3L1/results/plots/GPS2_CTCF_signal_on_common_genes_d6_heatmap.png",

        # 9m. Filtered GPS2 d6 promoter profiles
        "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_GPS2d6_promoter_heatmap.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_GPS2d6_promoter_profile.png",

        # 9n. ATF4-related profiles & intersections
        "CTCF_3T3L1/results/plots/CTCF_signal_on_ATF4_promoters_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_ATF4_profile.png",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/intersect/ctcf_atf4_overlap.bed",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_ATF4_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_promoter_only_ATF4_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_non_promoter_ATF4_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_3filtered_ATF4_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_3filtered_promoter_only_ATF4_profile.png",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_3filtered_non_promoter_ATF4_profile.png",
        "ATF4_3T3L1/results/plots/GPS2d6_CTCF_on_ATF4_promoter_heatmap.png",
        "ATF4_3T3L1/results/plots/GPS2d6_CTCF_on_ATF4_promoter_profile.png"


# Typical ChIPseq analysis workflow----------------------------------------------------------------------------------------------------------------  
# rule wget_files
rule wget_files:
    output:
        "{orig_dir}/raw_samples/{orig_file}.fastq.gz"
    params:
        link = lambda wildcards: sample_csv.loc[
            sample_csv['File name'] == wildcards.orig_file, 'FTP Link'
        ].iloc[0]
    threads: 8
    shell:
        '''
        mkdir -p $(dirname {output})
        wget -O {output} {params.link}
        '''


rule download_gtf:
    output:
        'Adapters_and_Annotations/GRCm39_annotation.fa'
    params:
        fa = 'Adapters_and_Annotations/GRCm39_annotation.fa.gz'
    shell:
        '''
        wget -O {params.fa} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.primary_assembly.genome.fa.gz
        gunzip {params.fa}
        '''
        
rule concatenate_fastq:
    input:
        lambda wildcards: [
            f"{wildcards.dir}/raw_samples/{fname.strip()}.fastq.gz"
            for fname in updated_samples_sheet.query('Name == @wildcards.sample')["Original_Files"].values[0].split(";")
        ]
    output:
        fastq = "{dir}/samples/{sample}.fastq.gz"
    shell:
        """
        mkdir -p $(dirname {output.fastq})
        cat {input} > {output.fastq}
        """

rule fastqc:
    input: 
        lambda wildcards: f"{wildcards.dir}/samples/{wildcards.file}.fastq.gz"
    output:
        html = '{dir}/results/fastqc/{file}_fastqc.html',
        zip = '{dir}/results/fastqc/{file}_fastqc.zip'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.html)
    threads: 8
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/fastqc_env.yml'
    shell:
        '''
        mkdir -p {params.outdir}
        fastqc {input} -o {params.outdir}
        '''

rule trimomatic:
    input:
        fastq = '{dir}/samples/{file}.fastq.gz'
    output:
        trimmed = '{dir}/results/trimmed/{file}.trim.fastq.gz'
    params:
        ad = 'Adapters_and_Annotations/Adapters.fa',
        outdir = lambda wildcards, output: os.path.dirname(output.trimmed)
    threads: 8
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/trimmomatic_env.yml'
    shell:
        '''
        mkdir -p {params.outdir}
        trimmomatic SE -threads {threads} -phred33 {input.fastq} {output.trimmed} ILLUMINACLIP:{params.ad}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        '''

rule bowtie2_build_gencode:
    input:
        fa = 'Adapters_and_Annotations/GRCm39_annotation.fa'
    output:
        expand('results/GRCm39.{ext}.bt2', ext=EXTENSIONS)
    params:
        outdir = 'results/GRCm39'
    threads: 16
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bowtie2_env.yml'
    shell:
        '''
        bowtie2-build {input.fa} {params.outdir}
        '''

rule unzip_fq:
    input:
        samples = '{dir}/results/trimmed/{file}.trim.fastq.gz'
    output:
        unzip = '{dir}/results/trimmed/{file}.trim.fastq'
    threads: 8
    shell:
        '''
        gunzip -c {input.samples} > {output.unzip}
        '''

rule trimqc:
    input: 
        fastq = '{dir}/results/trimmed/{file}.trim.fastq'
    output:
        fastqc = '{dir}/results/trimqc/{file}.trim_fastqc.html'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.fastqc)
    threads: 4
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/fastqc_env.yml'
    shell:
        '''
        mkdir -p {params.outdir}
        fastqc {input.fastq} -o {params.outdir}
        '''

rule bowtie2_align:
    input:
        samples = '{dir}/results/trimmed/{file}.trim.fastq',
        bt2 = expand('results/GRCm39.{ext}.bt2', ext=EXTENSIONS)
    output:
        bam = '{dir}/results/aligned/{file}.bam'
    threads: 16
    params:
        index = 'results/GRCm39',
        outdir = lambda wildcards, output: os.path.dirname(output.bam)
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bowtie2_env.yml'
    shell:
        '''
        mkdir -p {params.outdir}
        bowtie2 -x {params.index} -U {input.samples} | samtools view -bS - > {output.bam}
        '''

rule samtools_sort:
    input:
        bam = '{dir}/results/aligned/{file}.bam'
    output:
        sorted = '{dir}/results/sorted/{file}.sorted.bam'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.sorted)
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/samtools_env.yml'
    threads: 8
    shell:
        '''
        mkdir -p {params.outdir}
        samtools sort -o {output.sorted} {input.bam} 
        '''

rule samtools_idxstats:
    input:
        sorted = '{dir}/results/sorted/{file}.sorted.bam'
    output:
        index = '{dir}/results/view/{file}.txt'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.index)
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/samtools_env.yml'
    threads: 8
    shell:
        '''
        mkdir -p {params.outdir}
        samtools idxstats {input.sorted} > {output.index}
        '''

rule samtools_idx:
    input:
        filter = '{dir}/results/sorted/{file}.sorted.bam'
    output:
        index = '{dir}/results/sorted/{file}.sorted.bam.bai'
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/samtools_env.yml'
    threads: 8
    shell:
        '''
        samtools index {input.filter} {output.index}
        '''

rule samtools_flagstats:
    input:
        filter = '{dir}/results/sorted/{file}.sorted.bam'
    output:
        flagstat = '{dir}/results/flagstat/{file}.txt'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.flagstat)
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/samtools_env.yml'
    threads: 8
    
    shell:
        '''
        mkdir -p {params.outdir}
        samtools flagstat {input.filter} > {output.flagstat}
        '''

rule multiqc:
    input:
        expand("{dir}/results/flagstat/{file}.txt",
               zip, dir=updated_samples_sheet['Directory'], file=updated_samples_sheet['Name'])
    output:
        multiqc = expand("{directory}/multiqc/multiqc_report.html", directory=DIRECTORY)
    params:
        dirs = DIRECTORY
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/multiqc_env.yml"
    shell:
        """
        for dir in {params.dirs}; do 
            multiqc ${{dir}} -o ${{dir}}/multiqc
        done
        """

rule bamCoverage:
    input:
        bam = '{dir}/results/sorted/{file}.sorted.bam',
        bai = '{dir}/results/sorted/{file}.sorted.bam.bai',
        blacklist = 'Adapters_and_Annotations/mm39.excluderanges.bed'
    output:
        bigwig = '{dir}/results/bigwig/{file}.bw'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.bigwig)
    threads: 8
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml'
    shell:
        '''
        mkdir -p {params.outdir}
        bamCoverage -b {input.bam} -o {output.bigwig} -bl {input.blacklist} --effectiveGenomeSize 2654621783 --normalizeUsing RPKM -p {threads}
        '''

rule call_peaks_with_control:
    input:
        treatment = lambda wildcards: get_input_files(wildcards, 'Treatment', peak_calls_with_control),
        control = lambda wildcards: get_input_files(wildcards, 'Control', peak_calls_with_control)
    output:
        peaks = "{directory}/results/macs3/{peak_file}_peaks.narrowPeak"
    params:
        genome_size = "mm"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/macs3_env.yml'
    shell:
        """
        macs3 callpeak \
            -t {input.treatment} \
            -c {input.control} \
            -f BAM \
            -g {params.genome_size} \
            -n {wildcards.peak_file} \
            --outdir {wildcards.directory}/results/macs3
        """

rule call_peaks_without_control:
    input:
        treatment = lambda wildcards: get_input_files(wildcards, 'Treatment', peak_calls_without_control)
    output:
        peaks = "{directory}/results/macs3/{peak_file}_peaks.narrowPeak"
    params:
        genome_size = "mm"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/macs3_env.yml'
    shell:
        """
        macs3 callpeak \
            -t {input.treatment} \
            -f BAM \
            -g {params.genome_size} \
            -n {wildcards.peak_file} \
            --outdir {wildcards.directory}/results/macs3 \
            --nolambda
        """

rule filter_blacklist:
    input:
        peaks = '{directory}/results/macs3/{peak_file}_peaks.narrowPeak',
        blacklist = 'Adapters_and_Annotations/mm39.excluderanges.bed'
    output:
        filtered = '{directory}/results/filtered/{peak_file}_bf.narrowPeak'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.filtered)
    threads: 4
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml'
    shell:
        '''
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} > {output.filtered}
        '''

#----------------------------------------------------------------------------------------------------------------------------------------
#motif searh part
rule motifs:
    input:
        filtered = '{directory}/results/filtered/{peak_file}_bf.narrowPeak',
        fasta = 'Adapters_and_Annotations/GRCm39_annotation.fa'
    output:
        motifs = directory('{directory}/results/motifs/{peak_file}')
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml'
    threads: 8
    shell:
        '''
        findMotifsGenome.pl {input.filtered} {input.fasta} {output.motifs} -p {threads}
        '''

rule motifs_on_promoter:
    input:
        promoter_bed = '{directory}/results/annotation/{peak_file}_promoter_only.bed',
        fasta = 'Adapters_and_Annotations/GRCm39_annotation.fa'
    output:
        motifs = directory('{directory}/results/promoter_motifs/{peak_file}')
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml'
    threads: 8
    shell:
        '''
        findMotifsGenome.pl {input.promoter_bed} {input.fasta} {output.motifs} -size given -p {threads}
        '''

rule motifs_on_gps2_ctcf_intersect:
    input:
        bed = '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/CTCF_GPS2_overlap_union_t1to4_3T3L1_mm39.bed',
        fasta = 'Adapters_and_Annotations/GRCm39_annotation.fa'
    output:
        motifs = directory('/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/motifs/GPS2_CTCF_Intersect/')
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml'
    threads: 8
    shell:
        '''
        findMotifsGenome.pl {input.bed} {input.fasta} {output.motifs} -size given -p {threads}
        '''

rule motifs_on_CommonGeneFiltered_GPS2_non_promoter_peaks_and_promoter_only_peaks:
    input:
        bed_promoter = '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_filtered_by_GPS2_CTCF_common_promoter_only.bed',
        bed_nonpromoter = '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_filtered_by_GPS2_CTCF_common_non_promoter.bed',
        fasta = 'Adapters_and_Annotations/GRCm39_annotation.fa'
    output:
        motifs_promoter = directory('/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/motifs/Filtered_GPS2_Promoter'),
        motifs_nonpromoter= directory('/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/motifs/Filtered_GPS2_Non_Promoter')
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml'
    threads: 16
    shell:
        '''
        findMotifsGenome.pl {input.bed_promoter} {input.fasta} {output.motifs_promoter} -size given -p {threads}
        findMotifsGenome.pl {input.bed_nonpromoter} {input.fasta} {output.motifs_nonpromoter} -size given -p {threads}
        '''

rule motifs_on_filtered_gps2_flanks:
    input:
        bed = lambda wildcards: f"/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_filtered_by_GPS2_CTCF_common_slop{wildcards.size}bpflanks_only.bed",
        fasta = "Adapters_and_Annotations/GRCm39_annotation.fa"
    output:
        motifs = directory("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/motifs/Filtered_GPS2_flanks_{size}bp")
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml"
    threads: 16
    shell:
        """
        findMotifsGenome.pl {input.bed} {input.fasta} {output.motifs} -size given -p {threads}
        """

rule motifs_on_filtered_atf4_flanks:
    input:
        bed = lambda wildcards: f"/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_common_slop{wildcards.size}bpflanks_only.bed",
        fasta = "Adapters_and_Annotations/GRCm39_annotation.fa"
    output:
        motifs = directory("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/motifs/Filtered_ATF4_flanks_{size}bp")
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml"
    threads: 16
    shell:
        """
        findMotifsGenome.pl {input.bed} {input.fasta} {output.motifs} -size given -p {threads}
        """

rule motifs_on_filtered_gps2d6_flanks:
    input:
        bed = lambda wildcards: f"/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_filtered_by_CTCF_GPS2_d6_common_genes_slop{wildcards.size}bpflanks_only.bed",
        fasta = "Adapters_and_Annotations/GRCm39_annotation.fa"
    output:
        motifs = directory("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/motifs/Filtered_GPS2d6_flanks_{size}bp")
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml"
    threads: 16
    shell:
        """
        findMotifsGenome.pl {input.bed} {input.fasta} {output.motifs} -size given -p {threads}
        """

rule motifs_on_gps2_d6_center:
    input:
        bed = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_annotated_promoter_only.bed",
        fasta = "Adapters_and_Annotations/GRCm39_annotation.fa"
    output:
        motifs = directory("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/motifs/GPS2_d6_promoter_center")
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml"
    threads: 16
    shell:
        """
        findMotifsGenome.pl {input.bed} {input.fasta} {output.motifs} -size given -p {threads}
        """

rule motifs_on_atf4_center:
    input:
        bed = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_common_genes_promoter_only.bed",
        fasta = "Adapters_and_Annotations/GRCm39_annotation.fa"
    output:
        motifs = directory("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/motifs/ATF4_promoter_center")
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/homer_env.yml"
    threads: 16
    shell:
        """
        findMotifsGenome.pl {input.bed} {input.fasta} {output.motifs} -size given -p {threads}
        """

#-----------------------------------------------------------------------------------------------------------------------------

rule create_union_peakset:
    input:
        expand("CTCF_3T3L1/results/filtered/CTCF_t{n}_bf.narrowPeak", n=[1, 2, 3, 4])
    output:
        union_bed = "CTCF_3T3L1/results/filtered/CTCF_all_timepoints_union.bed"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml' 
    shell:
        '''
        cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output.union_bed}
        '''

rule intersect_ctcf_gps2:
    input:
        ctcf = "CTCF_3T3L1/results/filtered/CTCF_{tp}_bf.narrowPeak",
        gps2 = "Silencing/GPS2/results/filtered/sictl_bf.narrowPeak"
    output:
        overlap = "CTCF_3T3L1/results/intersect/CTCF_{tp}_GPS2_overlap.bed"
    conda: 
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml"
    shell:
        '''
        mkdir -p results/intersect
        bedtools intersect -a <(cut -f1-3 {input.ctcf}) \
                           -b <(cut -f1-3 {input.gps2}) -u \
                           > {output.overlap}
        '''

rule merge_ctcf_gps2_overlap_union:
    input:
        expand("CTCF_3T3L1/results/intersect/CTCF_t{tp}_GPS2_overlap.bed", tp=[1,2,3,4])
    output:
        "CTCF_3T3L1/results/intersect/CTCF_GPS2_overlap_union_t1to4_3T3L1_mm39.bed"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml"
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
        """


rule compute_matrix_CTCF_signal_on_GPS2_peaks:
    input:
        signal_files = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        reference_peaks = "Silencing/GPS2/results/filtered/sictl_bf.narrowPeak"
    output:
        matrix_file = "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_peaks.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8
    shell:
        '''
        mkdir -p $(dirname {output.matrix_file})
        computeMatrix reference-point \
            -S {input.signal_files} \
            -R {input.reference_peaks} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -out {output.matrix_file}
        '''

rule plot_CTCF_signal_profile_on_GPS2_peaks:
    input:
        matrix_file = rules.compute_matrix_CTCF_signal_on_GPS2_peaks.output.matrix_file
    output:
        plot_file = "CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_peaks_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot_file})
        plotProfile -m {input.matrix_file} \
            -out {output.plot_file} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal around GPS2 Binding Sites (siCtrl condition)"
        """

rule compute_matrix_CTCF_signal_on_CommonGeneFiltered_GPS2_peaks:
    input:
        signal_files = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        reference_peaks = "Silencing/GPS2/results/annotation/sictl_filtered_by_GPS2_CTCF_common.bed"
    output:
        matrix_file = "CTCF_3T3L1/results/matrix/CTCF_signal_over_CommonGeneFiltered_GPS2_peaks.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8 
    shell:
        '''
        mkdir -p $(dirname {output.matrix_file})
        computeMatrix reference-point \
            -S {input.signal_files} \
            -R {input.reference_peaks} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -out {output.matrix_file}
        '''

rule plot_CTCF_signal_profile_on_CommonGeneFiltered_GPS2_peaks:
    input:
        matrix_file = rules.compute_matrix_CTCF_signal_on_CommonGeneFiltered_GPS2_peaks.output.matrix_file
    output:
        plot_file = "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_peaks_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot_file})
        plotProfile -m {input.matrix_file} \
            --perGroup \
            --samplesLabel "CTCF_t1" "CTCF_t2" "CTCF_t3" "CTCF_t4" \
            --colors red blue green orange \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal around commongenes filtered GPS2 Peak (siCtrl condition)" \
            -out {output.plot_file}
        """


rule compute_matrix_CTCF_signal_on_GPS2_promoters:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        bed = "Silencing/GPS2/results/annotation/sictl_promoter_peaks.bed"
    output:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_promoters.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8 
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
            --binSize 20\
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix}
        """

rule plot_CTCF_signal_on_GPS2_promoters:
    input:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_promoters.gz"
    output:
        plot = "CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_promoters_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over GPS2 Promoter Peaks"
        """

rule compute_matrix_CTCF_signal_on_CommonGeneFiltered_GPS2_promoter_only_peaks:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        bed = "Silencing/GPS2/results/annotation/sictl_filtered_by_GPS2_CTCF_common_promoter_only.bed"
    output:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_CommonGeneFiltered_GPS2_promoter_only_peaks.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8 
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix}
        """

rule plot_CTCF_signal_on_CommonGeneFiltered_GPS2_promoter_only_peaks:
    input:
        matrix = rules.compute_matrix_CTCF_signal_on_CommonGeneFiltered_GPS2_promoter_only_peaks.output.matrix
    output:
        plot = "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_promoter_only_peaks_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over filtered GPS2 Promoter Peaks"
        """

rule compute_matrix_CTCF_signal_on_CommonGeneFiltered_GPS2_non_promoter_peaks:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        bed = "Silencing/GPS2/results/annotation/sictl_filtered_by_GPS2_CTCF_common_non_promoter.bed"
    output:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_CommonGeneFiltered_GPS2_non_promoter_peaks.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8 
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix}
        """

rule plot_CTCF_signal_on_CommonGeneFiltered_GPS2_non_promoter_peaks:
    input:
        matrix = rules.compute_matrix_CTCF_signal_on_CommonGeneFiltered_GPS2_non_promoter_peaks.output.matrix
    output:
        plot = "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_non_promoter_peaks_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over filtered GPS2 nonPromoter Peaks"
        """

rule compute_matrix_GPS2_day0_on_CTCF_promoters:
    input:
        bw = "Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw",
        bed = "CTCF_3T3L1/results/annotation/CTCF_t1_promoter_only.bed"
    output:
        matrix = "Adipocyte_differentiation/GPS2/results/matrix/gps2_day0_on_CTCF_t1_promoters.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8 
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix}
        """

rule plot_GPS2_day0_on_CTCF_promoters:
    input:
        matrix = "Adipocyte_differentiation/GPS2/results/matrix/gps2_day0_on_CTCF_t1_promoters.gz"
    output:
        plot = "Adipocyte_differentiation/GPS2/results/plots/gps2_day0_on_CTCF_t1_promoters_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "GPS2 Day 0 Signal on CTCF t1 Promoter Peaks"
        """

rule compute_matrix_CTCF_signal_on_d6_GPS2_peaks:
    input:
        signal_files = [
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        reference_peaks = "Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_bf.narrowPeak"
    output:
        matrix_file = "CTCF_3T3L1/results/matrix/CTCF_signal_on_d6_GPS2_peaks.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 16 
    shell:
        '''
        mkdir -p $(dirname {output.matrix_file})
        computeMatrix reference-point \
            -S {input.signal_files} \
            -R {input.reference_peaks} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -out {output.matrix_file}
        '''

rule plot_CTCF_signal_on_d6_GPS2_peaks:
    input:
        matrix_file = rules.compute_matrix_CTCF_signal_on_d6_GPS2_peaks.output.matrix_file
    output:
        plot_file = "CTCF_3T3L1/results/plots/CTCF_signal_on_d6_GPS2_peaks_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot_file})
        plotProfile -m {input.matrix_file} \
            -out {output.plot_file} \
            --legendLocation upper-right \
            --plotTitle "CTCF_signal_on_d6_GPS2_peaks"
        """

rule plot_heatmap_CTCF_signal_on_d6_GPS2_peaks:
    input:
        matrix_file = rules.compute_matrix_CTCF_signal_on_d6_GPS2_peaks.output.matrix_file
    output:
        heatmap_file = "CTCF_3T3L1/results/plots/CTCF_signal_on_d6_GPS2_peaks_heatmap.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.heatmap_file})
        plotHeatmap -m {input.matrix_file} \
            -out {output.heatmap_file} \
            --colorMap YlGnBu \
            --plotTitle "CTCF_signal_on_d6_GPS2_peaks"
        """

rule compute_matrix_GPS2_CTCF_d6_on_common_genes:
    input:
        signal_files = [
            "Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        reference_regions = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_CTCF_common_genes_filtered.bed"
    output:
        matrix_file = "CTCF_3T3L1/results/matrix/GPS2_CTCF_signal_on_common_genes_d6.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 16 
    shell:
        '''
        mkdir -p $(dirname {output.matrix_file})
        computeMatrix reference-point \
            -S {input.signal_files} \
            -R {input.reference_regions} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -out {output.matrix_file}
        '''

rule plot_heatmap_GPS2_CTCF_d6_on_common_genes:
    input:
        matrix_file = rules.compute_matrix_GPS2_CTCF_d6_on_common_genes.output.matrix_file
    output:
        heatmap_file = "CTCF_3T3L1/results/plots/GPS2_CTCF_signal_on_common_genes_d6_heatmap.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.heatmap_file})
        plotHeatmap -m {input.matrix_file} \
            -out {output.heatmap_file} \
            --colorMap YlGnBu \
            --plotTitle "GPS2 & CTCF on common genes (d6)"
        """

rule plot_profile_GPS2_CTCF_d6_on_common_genes:
    input:
        matrix_file = rules.compute_matrix_GPS2_CTCF_d6_on_common_genes.output.matrix_file
    output:
        profile_file = "CTCF_3T3L1/results/plots/GPS2_CTCF_signal_on_common_genes_d6_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.profile_file})
        plotProfile -m {input.matrix_file} \
            --samplesLabel "GPS2_d6" "CTCF_t3" "CTCF_t4" \
            --perGroup \
            --legendLocation upper-right \
            --plotTitle "GPS2 & CTCF on common genes (d6)" \
            -out {output.profile_file}
        """

rule compute_matrix_CTCF_on_filtered_GPS2d6_promoter:
    input:
        signal_files = [
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        reference_regions = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_filtered_by_CTCF_GPS2_d6_common_genes_promoter_only.bed"
    output:
        matrix_file = "CTCF_3T3L1/results/matrix/CTCF_signal_on_filtered_GPS2d6_promoter.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 16
    shell:
        """
        mkdir -p $(dirname {output.matrix_file})
        computeMatrix reference-point \
            -S {input.signal_files} \
            -R {input.reference_regions} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -out {output.matrix_file}
        """

rule plot_heatmap_CTCF_on_filtered_GPS2d6_promoter:
    input:
        matrix_file = rules.compute_matrix_CTCF_on_filtered_GPS2d6_promoter.output.matrix_file
    output:
        heatmap_file = "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_GPS2d6_promoter_heatmap.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.heatmap_file})
        plotHeatmap -m {input.matrix_file} \
            --colorMap YlGnBu \
            --plotTitle "CTCF on GPS2d6 common promoter" \
            -out {output.heatmap_file}
        """

rule plot_profile_CTCF_on_filtered_GPS2d6_promoter:
    input:
        matrix_file = rules.compute_matrix_CTCF_on_filtered_GPS2d6_promoter.output.matrix_file
    output:
        profile_file = "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_GPS2d6_promoter_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.profile_file})
        plotProfile -m {input.matrix_file} \
            --perGroup \
            --legendLocation upper-right \
            --plotTitle "CTCF on GPS2d6 common promoter" \
            -out {output.profile_file}
        """

#----------------------------------------------------------------------------
rule generate_genome_file:
    input:
        fasta="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adapters_and_Annotations/GRCm39_annotation.fa"
    output:
        genome="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adapters_and_Annotations/mm39.genome"
    threads: 8
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/samtools_env.yml'
    shell:
        """
        samtools faidx {input.fasta}
        cut -f1,2 {input.fasta}.fai > {output.genome}
        """

#------------------------------------------------------------------------------------------
#use siCTL for overlap with t1/t2 and d6 for t3/t4 (similar growth/differentiation stages)


rule merge_CTCF_preadipocyte:
    input:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_t1_bf.narrowPeak",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_t2_bf.narrowPeak"
    output:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_preadipocyte_merged.bed"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml'
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
        """

rule merge_CTCF_adipocyte:
    input:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_t3_bf.narrowPeak",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_t4_bf.narrowPeak"
    output:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_adipocyte_merged.bed"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml'
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
        """

rule overlap_preadipocyte:
    input:
        gps2 = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/filtered/sictl_bf.narrowPeak",
        ctcf = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_preadipocyte_merged.bed"
    output:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/GPS2_siCTL_vs_CTCF_preadipocyte.bed"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml'
    shell:
        """
        bedtools intersect -a {input.gps2} -b {input.ctcf} -wa -wb > {output}
        """

rule overlap_adipocyte:
    input:
        gps2 = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/filtered/sictl_bf.narrowPeak",
        ctcf = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_adipocyte_merged.bed"
    output:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/GPS2_d6_vs_CTCF_adipocyte.bed"
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml'
    shell:
        """
        bedtools intersect -a {input.gps2} -b {input.ctcf} -wa -wb > {output}
        """

#---------------------------------------------------------------------------------------------------------------------------

#ATF4 Dataset
rule compute_matrix_CTCF_signal_on_ATF4:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        bed = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/ATF4_d6_basal_mm39.narrowpeak"
    output:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_ATF4.gz"
    threads: 8 
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix}
        """

rule plot_CTCF_signal_on_ATF4:
    input:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_ATF4.gz"
    output:
        plot = "CTCF_3T3L1/results/plots/CTCF_signal_on_ATF4_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over ATF4 Peaks"
        """



rule compute_matrix_CTCF_signal_on_ATF4_promoters:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        bed = "ATF4_3T3L1/results/annotation/ATF4_d6_basal_promoter.bed"
    output:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_ATF4_promoters.gz"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8 
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix}
        """

rule plot_CTCF_signal_on_ATF4_promoters:
    input:
        matrix = "CTCF_3T3L1/results/matrix/CTCF_signal_over_ATF4_promoters.gz"
    output:
        plot = "CTCF_3T3L1/results/plots/CTCF_signal_on_ATF4_promoters_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over ATF4 Promoter Peaks"
        """

rule intersect_ctcf_atf4:
    input:
        ctcf_adipocyte = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_adipocyte_merged.bed",
        atf4 = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/ATF4_d6_basal_mm39.narrowpeak"
    output:
        overlap = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/intersect/ctcf_atf4_overlap.bed"
    conda: 
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml"
    shell:
        '''
        mkdir -p results/intersect
        bedtools intersect -a <(cut -f1-3 {input.ctcf_adipocyte}) \
                           -b <(cut -f1-3 {input.atf4}) -u \
                           > {output.overlap}
        '''
rule compute_matrix_all_CTCF_signal_on_filtered_ATF4:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        filtered_bed = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_common_genes.bed",
        promoter_only = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_common_genes_promoter_only.bed",
        non_promoter = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_common_genes_non_promoter.bed"
    output:
        matrix_filtered = "CTCF_3T3L1/results/matrix/CTCF_signal_over_filtered_ATF4.gz",
        matrix_promoter_only = "CTCF_3T3L1/results/matrix/CTCF_signal_over_filtered_promoter_only_ATF4.gz",
        matrix_non_promoter = "CTCF_3T3L1/results/matrix/CTCF_signal_over_filtered_non_promoter_ATF4.gz"
    threads: 16 
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.matrix_filtered})

        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.filtered_bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix_filtered}

        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.promoter_only} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix_promoter_only}

        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.non_promoter} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix_non_promoter}
        """


rule plot_all_CTCF_signal_on_filtered_ATF4:
    input:
        matrix_filtered = "CTCF_3T3L1/results/matrix/CTCF_signal_over_filtered_ATF4.gz",
        matrix_promoter_only = "CTCF_3T3L1/results/matrix/CTCF_signal_over_filtered_promoter_only_ATF4.gz",
        matrix_non_promoter = "CTCF_3T3L1/results/matrix/CTCF_signal_over_filtered_non_promoter_ATF4.gz"
    output:
        plot_filtered = "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_ATF4_profile.png",
        plot_promoter_only = "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_promoter_only_ATF4_profile.png",
        plot_non_promoter = "CTCF_3T3L1/results/plots/CTCF_signal_on_filtered_non_promoter_ATF4_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.plot_filtered})

        plotProfile -m {input.matrix_filtered} \
            -out {output.plot_filtered} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over All Filtered ATF4 Peaks"

        plotProfile -m {input.matrix_promoter_only} \
            -out {output.plot_promoter_only} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over Promoter-only ATF4 Peaks"

        plotProfile -m {input.matrix_non_promoter} \
            -out {output.plot_non_promoter} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over Non-promoter ATF4 Peaks"
        """

#3 means filtered by CTCF/ATF4/GPS2 common genes
rule compute_matrix_all_CTCF_signal_on_3filtered_ATF4:
    input:
        bw = [
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        filtered_bed = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_GPS2_common_genes.bed",
        promoter_only = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_GPS2_common_genes_promoter_only.bed",
        non_promoter = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_GPS2_common_genes_non_promoter.bed"
    output:
        matrix_filtered = "CTCF_3T3L1/results/matrix/CTCF_signal_over_3filtered_ATF4.gz",
        matrix_promoter_only = "CTCF_3T3L1/results/matrix/CTCF_signal_over_3filtered_promoter_only_ATF4.gz",
        matrix_non_promoter = "CTCF_3T3L1/results/matrix/CTCF_signal_over_3filtered_non_promoter_ATF4.gz"
    threads: 16 
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.matrix_filtered})

        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.filtered_bed} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix_filtered}

        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.promoter_only} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix_promoter_only}

        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.non_promoter} \
            --referencePoint center -b 3000 -a 3000 \
            -out {output.matrix_non_promoter}
        """


rule plot_all_CTCF_signal_on_3filtered_ATF4:
    input:
        matrix_filtered = "CTCF_3T3L1/results/matrix/CTCF_signal_over_3filtered_ATF4.gz",
        matrix_promoter_only = "CTCF_3T3L1/results/matrix/CTCF_signal_over_3filtered_promoter_only_ATF4.gz",
        matrix_non_promoter = "CTCF_3T3L1/results/matrix/CTCF_signal_over_3filtered_non_promoter_ATF4.gz"
    output:
        plot_filtered = "CTCF_3T3L1/results/plots/CTCF_signal_on_3filtered_ATF4_profile.png",
        plot_promoter_only = "CTCF_3T3L1/results/plots/CTCF_signal_on_3filtered_promoter_only_ATF4_profile.png",
        plot_non_promoter = "CTCF_3T3L1/results/plots/CTCF_signal_on_3filtered_non_promoter_ATF4_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.plot_filtered})

        plotProfile -m {input.matrix_filtered} \
            -out {output.plot_filtered} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over All 3Filtered ATF4 Peaks"

        plotProfile -m {input.matrix_promoter_only} \
            -out {output.plot_promoter_only} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over 3Filtered Promoter-only ATF4 Peaks"

        plotProfile -m {input.matrix_non_promoter} \
            -out {output.plot_non_promoter} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal over 3Filtered Non-promoter ATF4 Peaks"
        """

#compute matrix_ gps2/ctcf t3/t4 on atf4 signal
rule compute_matrix_gps2_ctcf_signal_on_ATF4:
    input:
        bw = [
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw",
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/bigwig/CTCF_t4.bw"
        ],
        bed = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_d6_basal_promoter.bed"
    output:
        matrix = "ATF4_3T3L1/results/matrix/GPS2d6_CTCF_on_ATF4_promoter_matrix.gz"
    threads: 8
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p ATF4_3T3L1/results/matrix
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -bs 50 \
            --skipZeros \
            --missingDataAsZero \
            -o {output.matrix}
        """

rule plot_gps2_ctcf_signal_on_ATF4_heatmap:
    input:
        matrix = "ATF4_3T3L1/results/matrix/GPS2d6_CTCF_on_ATF4_promoter_matrix.gz"
    output:
        png = "ATF4_3T3L1/results/plots/GPS2d6_CTCF_on_ATF4_promoter_heatmap.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p ATF4_3T3L1/results/plots
        plotHeatmap -m {input.matrix} \
            -out {output.png} \
            --colorMap RdBu_r \
            --refPointLabel "center" \
            --heatmapHeight 10 --heatmapWidth 6
        """

rule plot_gps2_ctcf_signal_on_ATF4_profile:
    input:
        matrix = "ATF4_3T3L1/results/matrix/GPS2d6_CTCF_on_ATF4_promoter_matrix.gz"
    output:
        png = "ATF4_3T3L1/results/plots/GPS2d6_CTCF_on_ATF4_promoter_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p ATF4_3T3L1/results/plots
        plotProfile -m {input.matrix} \
            -out {output.png} \
            --perGroup \
            --refPointLabel "center" \
            --plotTitle "GPS2/CTCF signal around ATF4 promoter"
        """
#ctcf_gps2d6 promoter_filtered heatmap
