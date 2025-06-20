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


rule all:
    input:
        # Download sample files
        expand("{orig_dir}/raw_samples/{orig_file}.fastq.gz",
               zip,
               orig_dir=sample_csv['Target Directory'],
               orig_file=sample_csv['File name']),

        # Files for downstream analysis
        expand("{dir}/samples/{file}.fastq.gz",
               zip,
               dir=updated_samples_sheet['Directory'],
               file=updated_samples_sheet['Name']),

        # Analysis outputs
        expand("{dir}/results/fastqc/{file}_fastqc.html", 
               zip, 
               dir=updated_samples_sheet['Directory'], 
               file=updated_samples_sheet['Name']),

        expand("{dir}/results/trimqc/{file}.trim_fastqc.html",
               zip, 
               dir=updated_samples_sheet['Directory'], 
               file=updated_samples_sheet['Name']),

        expand("{dir}/results/flagstat/{file}.txt",
               zip, 
               dir=updated_samples_sheet['Directory'], 
               file=updated_samples_sheet['Name']),

        expand("{directory}/multiqc/multiqc_report.html", directory=DIRECTORY),

        expand("{directory}/results/filtered/{peak_file}_bf.narrowPeak",
                zip, 
                directory=all_peak_calls['Directory'], 
                peak_file=all_peak_calls['Peak_file']),

        expand("{directory}/results/motifs/{peak_file}",
                zip, 
                directory=all_peak_calls['Directory'], 
                peak_file=all_peak_calls['Peak_file']),
        
        expand("{directory}/results/promoter_motifs/{peak_file}",
                zip, 
                directory=all_peak_calls['Directory'], 
                peak_file=all_peak_calls['Peak_file']),        

        expand("{dir}/results/bigwig/{file}.bw",
               zip,
               dir=updated_samples_sheet['Directory'],
               file=updated_samples_sheet['Name']),
        
        expand("CTCF_3T3L1/results/intersect/CTCF_{tp}_GPS2_overlap.bed", tp=CTCF_TIMEPOINTS),
        
       

        "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_promoters.gz",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_promoters_profile.png",
        #ComputeMatrix of CTCF signal on GPS2 Peak
        "CTCF_3T3L1/results/matrix/CTCF_signal_over_GPS2_peaks.gz",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_peaks_profile.png",
        "Silencing/GPS2/results/plots/GPS2_signal_on_CTCF_t1_promoters_profile.png",

        "Adipocyte_differentiation/GPS2/results/matrix/gps2_day0_on_CTCF_t1_promoters.gz",
        "Adipocyte_differentiation/GPS2/results/plots/gps2_day0_on_CTCF_t1_promoters_profile.png",

        expand(
            "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_extended_{size}.png",
            size=["2kb", "5kb", "10kb"]
        ),

        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_preadipocyte_merged.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/filtered/CTCF_adipocyte_merged.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/GPS2_siCTL_vs_CTCF_preadipocyte.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/intersect/GPS2_d6_vs_CTCF_adipocyte.bed",

        "CTCF_3T3L1/results/intersect/CTCF_GPS2_overlap_union_t1to4_3T3L1_mm39.bed",
        #computematrix ctcf over filtered gps2 peak
        "CTCF_3T3L1/results/matrix/CTCF_signal_over_CommonGeneFiltered_GPS2_peaks.gz",
        "CTCF_3T3L1/results/plots/CTCF_signal_on_CommonGeneFiltered_GPS2_peaks_profile.png"

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
            -out {output.plot_file} \
            --legendLocation upper-right \
            --plotTitle "CTCF Signal around commongenes filtered GPS2 Peak(siCtrl condition)"
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
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point \
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

rule compute_matrix_GPS2_day0_on_CTCF_promoters:
    input:
        bw = "Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw",
        bed = "CTCF_3T3L1/results/annotation/CTCF_t1_promoter_only.bed"
    output:
        matrix = "Adipocyte_differentiation/GPS2/results/matrix/gps2_day0_on_CTCF_t1_promoters.gz"
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

#----------------------------------------------------------------------------
#

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


rule extend_GPS2_peaks:
    input:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks.bed",
        genome="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adapters_and_Annotations/mm39.genome"
    output:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks_extended_2kb.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks_extended_5kb.bed",
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks_extended_10kb.bed"
    threads: 8
    conda:
        '/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/bedtools_env.yml'
    shell:
        """
        bedtools slop -i {input[0]} -g {input.genome} -b 2000 > {output[0]}
        bedtools slop -i {input[0]} -g {input.genome} -b 5000 > {output[1]}
        bedtools slop -i {input[0]} -g {input.genome} -b 10000 > {output[2]}
        """

rule compute_matrix_CTCF_on_GPS2_extended:
    input:
        bw=expand("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/bigwig/CTCF_t{i}.bw", i=[1, 2, 3, 4]),
        bed="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks_extended_{size}.bed"
    output:
        matrix="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/matrix/CTCF_signal_on_GPS2_extended_{size}.gz"
    threads: 16
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    params:
        flank=lambda wildcards: {"2kb": 2000, "5kb": 5000, "10kb": 10000}[wildcards.size]
    shell:
        """
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center \
            --beforeRegionStartLength {params.flank} \
            --afterRegionStartLength {params.flank} \
            --binSize 50 \
            --skipZeros \
            -o {output.matrix}
        """

rule plot_profile_CTCF_on_GPS2_extended:
    input:
        matrix="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/matrix/CTCF_signal_on_GPS2_extended_{size}.gz"
    output:
        png="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/plots/CTCF_signal_on_GPS2_extended_{size}.png"
    threads: 8
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        plotProfile -m {input.matrix} -out {output.png} \
            --perGroup \
            --plotTitle "CTCF signal ±{wildcards.size} from GPS2 peak center" \
            --refPointLabel "GPS2 peak"
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
