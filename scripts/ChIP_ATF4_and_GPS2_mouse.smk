import pandas as pd
import os

sample_csv = pd.read_csv('Sample_Data_for_GPS2_and_ATF4.csv')
updated_samples_sheet = pd.read_csv('Updated_Sample_Sheet.csv')
peak_calls_with_control = pd.read_csv("Peak_Calls_with_Control.csv")
peak_calls_without_control = pd.read_csv("Peak_Calls_without_Control.csv")

# Add a 'Control' column to the without_control DataFrame (fill with 'None')
peak_calls_without_control['Control'] = 'None'

# Combine both DataFrames
all_peak_calls = pd.concat([peak_calls_with_control, peak_calls_without_control], ignore_index=True)

REPEATS = ['rep1', 'rep2', 'rep3']
EXTENSIONS = [1, 2, 3, 4, 'rev.1', 'rev.2']
DIRECTORY = ['Silencing', 'Mito_stress_treatment', 'Adipocyte_differentiation']





def get_input_files(wildcards, column, dataframe):
    files = dataframe.loc[dataframe['Peak_file'] == wildcards.peak_file, column].values[0]
    return [f"{wildcards.directory}/results/sorted/{file.strip()}.sorted.bam" for file in files.split(';')]

rule all:
    input:
        # Download sample files
        expand("{orig_dir}/samples/{orig_file}.fastq.gz", 
               zip, orig_dir=sample_csv['Target Directory'], orig_file=sample_csv['File name']),
        # Files for downstream analysis (including concatenated)
        expand("{dir}/samples/{file}.fastq.gz", 
               zip, dir=updated_samples_sheet['Directory'], file=updated_samples_sheet['Name']),
        # Analysis outputs
        expand("{dir}/results/fastqc/{file}_fastqc.html", 
               zip, dir=updated_samples_sheet['Directory'], file=updated_samples_sheet['Name']),
        expand("{dir}/results/trimqc/{file}.trim_fastqc.html",
               zip, dir=updated_samples_sheet['Directory'], file=updated_samples_sheet['Name']),
        expand("{dir}/results/flagstat/{file}.txt",
               zip, dir=updated_samples_sheet['Directory'], file=updated_samples_sheet['Name']),
        expand("{directory}/multiqc/multiqc_report.html", directory=DIRECTORY),
        expand("{directory}/results/filtered/{peak_file}_bf.narrowPeak",
                zip, directory=all_peak_calls['Directory'], peak_file=all_peak_calls['Peak_file']),
        expand("{directory}/results/motifs/{peak_file}",
                zip, directory=all_peak_calls['Directory'], peak_file=all_peak_calls['Peak_file']),
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Differentiation_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Treatment_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Silencing_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Differentiation_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Treatment_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Silencing_plotprofile.png',
        #RNA_seq_comparisons
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Differentiation_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Treatment_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Silencing_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Differentiation_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Treatment_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Silencing_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Treatment_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Treatment_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Silencing_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Treatment_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Treatment_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Silencing_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Differentiation_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Treatment_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Silencing_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Differentiation_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Treatment_plotprofile.png',
        'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Silencing_plotprofile.png',
        #Plotprofile_over_TSS
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Differentiation1_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Silencing_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Treatment_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Differentiation_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Silencing_plotprofile.png',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Treatment_plotprofile.png',
        #Overlaps_with_AAD
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_aad.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_aad.bed',
        #Overlaps_with_WT_TMC
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_wt_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_wt_tmc.bed',
        #Overlaps_withKO_TMC
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_ko_tmc.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_ko_tmc.bed',
        #Overlaps_with_WT_vs_KO
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_wt_vs_ko.bed',
        'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_wt_vs_ko.bed'
 
rule wget_files:
    output:
        "{orig_dir}/samples/{orig_file}.fastq.gz"
    params:
        link = lambda wildcards: sample_csv.loc[
            (sample_csv['Target Directory'] == wildcards.orig_dir) &
            (sample_csv['File name'] == wildcards.orig_file), 'FTP Link'
        ].iloc[0]
    threads: 8
    shell:
        '''
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
        ATF4_WT = expand('Mito_stress_treatment/ER_stress/ATF4/samples/ATF4_WT_{rep}.fastq.gz', rep = REPEATS),
        ATF4_KO = expand('Mito_stress_treatment/ER_stress/ATF4/samples/ATF4_KO_{rep}.fastq.gz', rep = REPEATS)
    output:
        ATF4_WT = 'Mito_stress_treatment/ER_stress/ATF4/samples/ATF4_WT_Tmc.fastq.gz',
        ATF4_KO = 'Mito_stress_treatment/ER_stress/ATF4/samples/ATF4_KO_Tmc.fastq.gz'
    shell:
        '''
        cat {input.ATF4_WT} > {output.ATF4_WT}
        cat {input.ATF4_KO} > {output.ATF4_KO}
        '''

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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/fastqc_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/trimmomatic_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bowtie2_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/fastqc_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bowtie2_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/samtools_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/samtools_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/samtools_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/samtools_env.yml'
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
        "/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/multiqc_env.yml"
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
        blacklist = 'samples/mm39.excluderanges.bed'
    output:
        bigwig = '{dir}/results/bigwig/{file}.bw'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.bigwig)
    threads: 8
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/macs3_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/macs3_env.yml'
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
        blacklist = 'samples/mm39.excluderanges.bed'
    output:
        filtered = '{directory}/results/filtered/{peak_file}_bf.narrowPeak'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.filtered)
    threads: 4
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bedtools_env.yml'
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
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/homer_env.yml'
    threads: 8
    shell:
        '''
        findMotifsGenome.pl {input.filtered} {input.fasta} {output.motifs} -p {threads}
        '''
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

rule computeMatrix_Differentiation:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/filtered/ATF4_AAD_bf.narrowPeak',
        ATF4_ER = 'Mito_stress_treatment/ER_stress/ATF4/results/filtered/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_bf.narrowPeak'
    output:
        ATF4_AAD = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Differentiation_matrix.gz',
        ATF4_ER = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Differentiation_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_AAD} --outFileName {output.ATF4_AAD} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_ER} --outFileName {output.ATF4_ER} -p {threads}
        '''

rule computeMatrix_Treatment:
    input:
        gps2_with_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep1.bw',
        gps2_with_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep2.bw',
        gps2_no_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep1.bw',
        gps2_no_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/filtered/ATF4_AAD_bf.narrowPeak',
        ATF4_ER = 'Mito_stress_treatment/ER_stress/ATF4/results/filtered/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_bf.narrowPeak'
    output:
        ATF4_AAD = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Treatment_matrix.gz',
        ATF4_ER = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Treatment_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_AAD} --outFileName {output.ATF4_AAD} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_ER} --outFileName {output.ATF4_ER} -p {threads}
        '''

rule computeMatrix_Silencing:
    input:
        sictl = 'Silencing/GPS2/results/bigwig/chip_sictl.bw',
        sigps2 = 'Silencing/GPS2/results/bigwig/chip_sigps2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/filtered/ATF4_AAD_bf.narrowPeak',
        ATF4_ER = 'Mito_stress_treatment/ER_stress/ATF4/results/filtered/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_bf.narrowPeak'
    output:
        ATF4_AAD = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Silencing_matrix.gz',
        ATF4_ER = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Silencing_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.sictl} {input.sigps2} -R {input.ATF4_AAD} --outFileName {output.ATF4_AAD} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.sictl} {input.sigps2} -R {input.ATF4_ER} --outFileName {output.ATF4_ER} -p {threads}
        '''

rule Plotprofile:
    input:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Differentiation_matrix.gz',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Treatment_matrix.gz',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Silencing_matrix.gz',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Differentiation_matrix.gz',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Treatment_matrix.gz',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Silencing_matrix.gz'
    output:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Differentiation_plotprofile.png',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Treatment_plotprofile.png',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/ATF4_AAD_Silencing_plotprofile.png',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Differentiation_plotprofile.png',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Treatment_plotprofile.png',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/ATF4_ER_Silencing_plotprofile.png'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    shell:
        '''
        plotProfile -m {input.diff_aad} -o {output.diff_aad} --perGroup --plotTitle "ATF4_AAD with GPS2 Differentiation"
        plotProfile -m {input.treat_aad} -o {output.treat_aad} --perGroup --plotTitle "ATF4_AAD with GPS2 Treatment"
        plotProfile -m {input.sil_aad} -o {output.sil_aad} --perGroup --plotTitle "ATF4_AAD with GPS2 Silencing"
        plotProfile -m {input.diff_er} -o {output.diff_er} --perGroup --plotTitle "ATF4_ER with GPS2 Differentiation"
        plotProfile -m {input.treat_er} -o {output.treat_er} --perGroup --plotTitle "ATF4_ER with GPS2 Treatment"
        plotProfile -m {input.sil_er} -o {output.sil_er} --perGroup --plotTitle "ATF4_ER with GPS2 Silencing"
        '''

rule computeMatrix_Treatment_over_RNAseq_GSE76771:
    input:
        gps2_with_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep1.bw',
        gps2_with_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep2.bw',
        gps2_no_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep1.bw',
        gps2_no_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_WTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_KO.bed',
        ATF4_AAD_WTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_KO-TT.bed',
        ATF4_AAD_WTvsWT_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_WT-TT.bed',
        ATF4_AAD_WT_TTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT-TT_vs_KO.bed',
        ATF4_AAD_WT_TTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT-TT_vs_KO-TT.bed',
        ATF4_AAD_KOvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_KO_vs_KO-TT.bed',
        ATF4_ER_WTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_KO.bed',
        ATF4_ER_WTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_KO-TT.bed',
        ATF4_ER_WTvsWT_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_WT-TT.bed',
        ATF4_ER_WT_TTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT-TT_vs_KO.bed',
        ATF4_ER_WT_TTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT-TT_vs_KO-TT.bed',
        ATF4_ER_KOvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_KO_vs_KO-TT.bed',

    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Treatment_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Treatment_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_AAD_WTvsKO} {input.ATF4_AAD_WTvsWT_TT} {input.ATF4_AAD_KOvsKO_TT} {input.ATF4_AAD_WTvsKO_TT} {input.ATF4_AAD_WT_TTvsKO} {input.ATF4_AAD_WT_TTvsKO_TT} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_ER_WTvsKO} {input.ATF4_ER_WTvsWT_TT} {input.ATF4_ER_KOvsKO_TT} {input.ATF4_ER_WTvsKO_TT} {input.ATF4_ER_WT_TTvsKO} {input.ATF4_ER_WT_TTvsKO_TT} --outFileName {output.er} -p {threads}
        '''

rule computeMatrix_Differentiation_over_RNAseq_GSE76771:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_WTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_KO.bed',
        ATF4_AAD_WTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_KO-TT.bed',
        ATF4_AAD_WTvsWT_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_WT-TT.bed',
        ATF4_AAD_WT_TTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT-TT_vs_KO.bed',
        ATF4_AAD_WT_TTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT-TT_vs_KO-TT.bed',
        ATF4_AAD_KOvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_KO_vs_KO-TT.bed',
        ATF4_ER_WTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_KO.bed',
        ATF4_ER_WTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_KO-TT.bed',
        ATF4_ER_WTvsWT_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_WT-TT.bed',
        ATF4_ER_WT_TTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT-TT_vs_KO.bed',
        ATF4_ER_WT_TTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT-TT_vs_KO-TT.bed',
        ATF4_ER_KOvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_KO_vs_KO-TT.bed',
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Differentiation_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Differentiation_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_AAD_WTvsKO} {input.ATF4_AAD_WTvsWT_TT} {input.ATF4_AAD_KOvsKO_TT} {input.ATF4_AAD_WTvsKO_TT} {input.ATF4_AAD_WT_TTvsKO} {input.ATF4_AAD_WT_TTvsKO_TT} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_ER_WTvsKO} {input.ATF4_ER_WTvsWT_TT} {input.ATF4_ER_KOvsKO_TT} {input.ATF4_ER_WTvsKO_TT} {input.ATF4_ER_WT_TTvsKO} {input.ATF4_ER_WT_TTvsKO_TT} --outFileName {output.er} -p {threads}
        '''

rule computeMatrix_Silencing_over_RNAseq_GSE76771:
    input:
        sictl = 'Silencing/GPS2/results/bigwig/chip_sictl.bw',
        sigps2 = 'Silencing/GPS2/results/bigwig/chip_sigps2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_WTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_KO.bed',
        ATF4_AAD_WTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_KO-TT.bed',
        ATF4_AAD_WTvsWT_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT_vs_WT-TT.bed',
        ATF4_AAD_WT_TTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT-TT_vs_KO.bed',
        ATF4_AAD_WT_TTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_WT-TT_vs_KO-TT.bed',
        ATF4_AAD_KOvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/AAD_overlaps/ATF4_AAD_peaks_overlap_with_KO_vs_KO-TT.bed',
        ATF4_ER_WTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_KO.bed',
        ATF4_ER_WTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_KO-TT.bed',
        ATF4_ER_WTvsWT_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT_vs_WT-TT.bed',
        ATF4_ER_WT_TTvsKO = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT-TT_vs_KO.bed',
        ATF4_ER_WT_TTvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_WT-TT_vs_KO-TT.bed',
        ATF4_ER_KOvsKO_TT = 'ATF4_RNAseq/WT_and_KO_Tmc_treatment_LIVER_26960794/ER_overlaps/ATF4_ER_peaks_overlap_with_KO_vs_KO-TT.bed',
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Silencing_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Silencing_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.sictl} {input.sigps2} -R {input.ATF4_AAD_WTvsKO} {input.ATF4_AAD_WTvsWT_TT} {input.ATF4_AAD_KOvsKO_TT} {input.ATF4_AAD_WTvsKO_TT} {input.ATF4_AAD_WT_TTvsKO} {input.ATF4_AAD_WT_TTvsKO_TT} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.sictl} {input.sigps2} -R {input.ATF4_ER_WTvsKO} {input.ATF4_ER_WTvsWT_TT} {input.ATF4_ER_KOvsKO_TT} {input.ATF4_ER_WTvsKO_TT} {input.ATF4_ER_WT_TTvsKO} {input.ATF4_ER_WT_TTvsKO_TT} --outFileName {output.er} -p {threads}
        '''

rule Plotprofile_RNAseq_GSE76771:
    input:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Differentiation_matrix.gz',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Treatment_matrix.gz',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Silencing_matrix.gz',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Differentiation_matrix.gz',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Treatment_matrix.gz',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Silencing_matrix.gz'
    output:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Differentiation_plotprofile.png',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Treatment_plotprofile.png',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4Tmc_Silencing_plotprofile.png',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Differentiation_plotprofile.png',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Treatment_plotprofile.png',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4Tmc_Silencing_plotprofile.png'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    shell:
        '''
        plotProfile -m {input.diff_aad} -o {output.diff_aad} --perGroup --plotTitle "ATF4_AAD with Tmc across GPS2 Differentiation"
        plotProfile -m {input.treat_aad} -o {output.treat_aad} --perGroup --plotTitle "ATF4_AAD with Tmc across GPS2 Treatment"
        plotProfile -m {input.sil_aad} -o {output.sil_aad} --perGroup --plotTitle "ATF4_AAD with Tmc across GPS2 Silencing"
        plotProfile -m {input.diff_er} -o {output.diff_er} --perGroup --plotTitle "ATF4_ER with Tmc across GPS2 Differentiation"
        plotProfile -m {input.treat_er} -o {output.treat_er} --perGroup --plotTitle "ATF4_ER with Tmc across GPS2 Treatment"
        plotProfile -m {input.sil_er} -o {output.sil_er} --perGroup --plotTitle "ATF4_ER with Tmc across GPS2 Silencing"
        '''

rule computeMatrix_Treatment_over_RNAseq_mTORC1_GSE158605:
    input:
        gps2_with_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep1.bw',
        gps2_with_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep2.bw',
        gps2_no_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep1.bw',
        gps2_no_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_mTORC1_IR = 'ATF4_RNAseq/mTORC1_33646118/AAD/ATF4_AAD_from_mTORC1_ATF4-IR.bed',
        ATF4_AAD_mTORC1_TT = 'ATF4_RNAseq/mTORC1_33646118/AAD/ATF4_AAD_from_mTORC1_ATF4-TT.bed',
        ATF4_ER_MTORC1_IR = 'ATF4_RNAseq/mTORC1_33646118/ER/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_from_mTORC1_ATF4-IR.bed',
        ATF4_ER_MTORC1_TT = 'ATF4_RNAseq/mTORC1_33646118/ER/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_from_mTORC1_ATF4-TT.bed'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Treatment_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Treatment_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_AAD_mTORC1_IR} {input.ATF4_AAD_mTORC1_TT} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_ER_MTORC1_IR} {input.ATF4_ER_MTORC1_TT} --outFileName {output.er} -p {threads}
        '''

rule computeMatrix_Differentiation_over_RNAseq_mTORC1_GSE158605:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_mTORC1_IR = 'ATF4_RNAseq/mTORC1_33646118/AAD/ATF4_AAD_from_mTORC1_ATF4-IR.bed',
        ATF4_AAD_mTORC1_TT = 'ATF4_RNAseq/mTORC1_33646118/AAD/ATF4_AAD_from_mTORC1_ATF4-TT.bed',
        ATF4_ER_MTORC1_IR = 'ATF4_RNAseq/mTORC1_33646118/ER/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_from_mTORC1_ATF4-IR.bed',
        ATF4_ER_MTORC1_TT = 'ATF4_RNAseq/mTORC1_33646118/ER/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_from_mTORC1_ATF4-TT.bed'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Differentiation_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Differentiation_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_AAD_mTORC1_IR} {input.ATF4_AAD_mTORC1_TT} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_ER_MTORC1_IR} {input.ATF4_ER_MTORC1_TT} --outFileName {output.er} -p {threads}
        '''

rule computeMatrix_Silencing_over_RNAseq_mTORC1_GSE158605:
    input:
        sictl = 'Silencing/GPS2/results/bigwig/chip_sictl.bw',
        sigps2 = 'Silencing/GPS2/results/bigwig/chip_sigps2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_mTORC1_IR = 'ATF4_RNAseq/mTORC1_33646118/AAD/ATF4_AAD_from_mTORC1_ATF4-IR.bed',
        ATF4_AAD_mTORC1_TT = 'ATF4_RNAseq/mTORC1_33646118/AAD/ATF4_AAD_from_mTORC1_ATF4-TT.bed',
        ATF4_ER_MTORC1_IR = 'ATF4_RNAseq/mTORC1_33646118/ER/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_from_mTORC1_ATF4-IR.bed',
        ATF4_ER_MTORC1_TT = 'ATF4_RNAseq/mTORC1_33646118/ER/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_from_mTORC1_ATF4-TT.bed'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Silencing_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Silencing_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.sictl} {input.sigps2} -R {input.ATF4_AAD_mTORC1_IR} {input.ATF4_AAD_mTORC1_TT} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.sictl} {input.sigps2} -R {input.ATF4_ER_MTORC1_IR} {input.ATF4_ER_MTORC1_TT} --outFileName {output.er} -p {threads}
        '''

rule Plotprofile_RNAseq_mTORC1_GSE15860:
    input:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Differentiation_matrix.gz',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Treatment_matrix.gz',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Silencing_matrix.gz',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Differentiation_matrix.gz',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Treatment_matrix.gz',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Silencing_matrix.gz'
    output:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Differentiation_plotprofile.png',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Treatment_plotprofile.png',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_across_ATF4mTORC1_Silencing_plotprofile.png',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Differentiation_plotprofile.png',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Treatment_plotprofile.png',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_across_ATF4mTORC1_Silencing_plotprofile.png'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    shell:
        '''
        plotProfile -m {input.diff_aad} -o {output.diff_aad} --perGroup --plotTitle "ATF4_AAD in mTORC1 across GPS2 Differentiation"
        plotProfile -m {input.treat_aad} -o {output.treat_aad} --perGroup --plotTitle "ATF4_AAD in mTORC1 across GPS2 Treatment"
        plotProfile -m {input.sil_aad} -o {output.sil_aad} --perGroup --plotTitle "ATF4_AAD in mTORC1 across GPS2 Silencing"
        plotProfile -m {input.diff_er} -o {output.diff_er} --perGroup --plotTitle "ATF4_ER in mTORC1 across GPS2 Differentiation"
        plotProfile -m {input.treat_er} -o {output.treat_er} --perGroup --plotTitle "ATF4_ER in mTORC1 across GPS2 Treatment"
        plotProfile -m {input.sil_er} -o {output.sil_er} --perGroup --plotTitle "ATF4_ER in mTORC1 across GPS2 Silencing"
        '''
##
rule computeMatrix_Treatment_over_RNAseq_multiomics_GSE84631:
    input:
        gps2_with_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep1.bw',
        gps2_with_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep2.bw',
        gps2_no_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep1.bw',
        gps2_no_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_DE = 'ATF4_RNAseq/Multiomics_28566324/ATF4_AAD_overlaps_with_DEgenes_Multiomics.bed',
        ATF4_ER_DE = 'ATF4_RNAseq/Multiomics_28566324/ATF4_ER_overlaps_with_DEgenes_Multiomics.bed'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Treatment_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Treatment_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_AAD_DE} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.ATF4_ER_DE} --outFileName {output.er} -p {threads}
        '''

rule computeMatrix_Differentiation_over_RNAseq_multiomics_GSE84631:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_DE = 'ATF4_RNAseq/Multiomics_28566324/ATF4_AAD_overlaps_with_DEgenes_Multiomics.bed',
        ATF4_ER_DE = 'ATF4_RNAseq/Multiomics_28566324/ATF4_ER_overlaps_with_DEgenes_Multiomics.bed'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Differentiation_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Differentiation_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_AAD_DE} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_day0} {input.gps2_day6} -R {input.ATF4_ER_DE} --outFileName {output.er} -p {threads}
        '''

rule computeMatrix_Silencing_over_RNAseq_multiomics_GSE84631:
    input:
        sictl = 'Silencing/GPS2/results/bigwig/chip_sictl.bw',
        sigps2 = 'Silencing/GPS2/results/bigwig/chip_sigps2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        ATF4_AAD_DE = 'ATF4_RNAseq/Multiomics_28566324/ATF4_AAD_overlaps_with_DEgenes_Multiomics.bed',
        ATF4_ER_DE = 'ATF4_RNAseq/Multiomics_28566324/ATF4_ER_overlaps_with_DEgenes_Multiomics.bed'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Silencing_matrix.gz',
        er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Silencing_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint center -S {input.ATF4_AAD_bw} {input.sictl} {input.sigps2} -R {input.ATF4_AAD_DE} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint center -S {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.sictl} {input.sigps2} -R {input.ATF4_ER_DE} --outFileName {output.er} -p {threads}
        '''

rule Plotprofile_RNAseq_multiomics_GSE84631:
    input:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Differentiation_matrix.gz',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Treatment_matrix.gz',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Silencing_matrix.gz',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Differentiation_matrix.gz',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Treatment_matrix.gz',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Silencing_matrix.gz'
    output:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Differentiation_plotprofile.png',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Treatment_plotprofile.png',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/ATF4_AAD_overlap_Multiomics_Silencing_plotprofile.png',
        diff_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Differentiation_plotprofile.png',
        treat_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Treatment_plotprofile.png',
        sil_er = 'Mito_stress_treatment/ER_stress/ATF4/results/matrix/RNA_seq/ATF4_ER_overlap_Multiomics_Silencing_plotprofile.png'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    shell:
        '''
        plotProfile -m {input.diff_aad} -o {output.diff_aad} --perGroup --plotTitle "ATF4_AAD from Multiomics across GPS2 Differentiation"
        plotProfile -m {input.treat_aad} -o {output.treat_aad} --perGroup --plotTitle "ATF4_AAD from Multiomics across GPS2 Treatment"
        plotProfile -m {input.sil_aad} -o {output.sil_aad} --perGroup --plotTitle "ATF4_AAD from Multiomics across GPS2 Silencing"
        plotProfile -m {input.diff_er} -o {output.diff_er} --perGroup --plotTitle "ATF4_ER from Multiomics across GPS2 Differentiation"
        plotProfile -m {input.treat_er} -o {output.treat_er} --perGroup --plotTitle "ATF4_ER from Multiomics across GPS2 Treatment"
        plotProfile -m {input.sil_er} -o {output.sil_er} --perGroup --plotTitle "ATF4_ER from Multiomics across GPS2 Silencing"
        '''

rule computeMatrix_Treatment_over_TSS_of_RNAseq_multiomics_GSE84631:
    input:
        gps2_with_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep1.bw',
        gps2_with_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_withFCCP_rep2.bw',
        gps2_no_rep1 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep1.bw',
        gps2_no_rep2 = 'Mito_stress_treatment/GPS2/results/bigwig/gps2_noFCCP_rep2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        All_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_All_DEgenes_Genomicregions.txt',
        RNAseq_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_RNAseq_DE_genes_genomicregions.txt',
        Prot_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_Proteomics_DE_genes_genomicregions.txt',
        Alz_regions = 'ATF4_RNAseq/Alzheimers_36867706/Alzheimers_genomic_regions.txt'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Treatment_matrix.gz',
        alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Treatment_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint TSS -S {input.ATF4_AAD_bw} {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.All_DE_regions} {input.RNAseq_DE_regions} {input.Prot_DE_regions} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint TSS -S {input.ATF4_AAD_bw} {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_with_rep1} {input.gps2_with_rep2} {input.gps2_no_rep1} {input.gps2_no_rep2} -R {input.Alz_regions} --outFileName {output.alz} -p {threads}
        '''

rule computeMatrix_Differentiation_over_TSS_of_RNAseq_multiomics_GSE84631:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day0.bw',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/bigwig/gps2_day6.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        All_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_All_DEgenes_Genomicregions.txt',
        RNAseq_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_RNAseq_DE_genes_genomicregions.txt',
        Prot_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_Proteomics_DE_genes_genomicregions.txt',
        Alz_regions = 'ATF4_RNAseq/Alzheimers_36867706/Alzheimers_genomic_regions.txt'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Differentiation1_matrix.gz',
        alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Differentiation_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint TSS -S {input.ATF4_AAD_bw} {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_day0} {input.gps2_day6} -R {input.All_DE_regions} {input.RNAseq_DE_regions} {input.Prot_DE_regions} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint TSS -S {input.ATF4_AAD_bw} {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.gps2_day0} {input.gps2_day6} -R {input.Alz_regions} --outFileName {output.alz} -p {threads}
        '''

rule computeMatrix_Silencing_over_TSS_of_RNAseq_multiomics_GSE84631:
    input:
        sictl = 'Silencing/GPS2/results/bigwig/chip_sictl.bw',
        sigps2 = 'Silencing/GPS2/results/bigwig/chip_sigps2.bw',
        ATF4_AAD_bw = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/bigwig/ATF4_AAD.bw',
        ATF4_WT_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_WT_Tmc.bw',
        ATF4_KO_ER_bw = 'Mito_stress_treatment/ER_stress/ATF4/results/bigwig/ATF4_KO_Tmc.bw',
        All_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_All_DEgenes_Genomicregions.txt',
        RNAseq_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_RNAseq_DE_genes_genomicregions.txt',
        Prot_DE_regions = 'ATF4_RNAseq/Multiomics_28566324/Multiomics_Proteomics_DE_genes_genomicregions.txt',
        Alz_regions = 'ATF4_RNAseq/Alzheimers_36867706/Alzheimers_genomic_regions.txt'
    output:
        aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Silencing_matrix.gz',
        alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Silencing_matrix.gz'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    threads: 16
    shell:
        '''
        computeMatrix reference-point --referencePoint TSS -S {input.ATF4_AAD_bw} {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.sictl} {input.sigps2} -R {input.All_DE_regions} {input.RNAseq_DE_regions} {input.Prot_DE_regions} --outFileName {output.aad} -p {threads}
        computeMatrix reference-point --referencePoint TSS -S {input.ATF4_AAD_bw} {input.ATF4_WT_ER_bw} {input.ATF4_KO_ER_bw} {input.sictl} {input.sigps2} -R {input.Alz_regions} --outFileName {output.alz} -p {threads}
        '''

rule Plotprofile_TSS_of_RNAseq_multiomics_GSE84631:
    input:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Differentiation1_matrix.gz',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Silencing_matrix.gz',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Treatment_matrix.gz',
        diff_alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Differentiation_matrix.gz',
        treat_alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Treatment_matrix.gz',
        sil_alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Silencing_matrix.gz'
    output:
        diff_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Differentiation1_plotprofile.png',
        sil_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Silencing_plotprofile.png',
        treat_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Multiomics_DE_regions_overTSS_Treatment_plotprofile.png',
        diff_alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Differentiation_plotprofile.png',
        treat_alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Treatment_plotprofile.png',
        sil_alz = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/matrix/RNA_seq/Alzheimers_DE_regions_overTSS_Silencing_plotprofile.png'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/deeptools_env.yml'
    shell:
        '''
        plotProfile -m {input.diff_aad} -o {output.diff_aad} --perGroup --plotTitle "Differentiation - Centered across gene TSS from Multiomics"
        plotProfile -m {input.sil_aad} -o {output.sil_aad} --perGroup --plotTitle "Silencing - Centered across gene TSS from Multiomics"
        plotProfile -m {input.treat_aad} -o {output.treat_aad} --perGroup --plotTitle "Treatment - Centered across gene TSS from Multiomics"
        plotProfile -m {input.diff_alz} -o {output.diff_alz} --perGroup --plotTitle "Differentiation - Centered across gene TSS from Alzheimers"
        plotProfile -m {input.treat_alz} -o {output.treat_alz} --perGroup --plotTitle "Treatment - Centered across gene TSS from Alzheimers"
        plotProfile -m {input.sil_alz} -o {output.sil_alz} --perGroup --plotTitle "Silencing - Centered across gene TSS from Alzheimers"
        '''

rule overlaps_AAD:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day0_bf.narrowPeak',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_bf.narrowPeak',
        day6_vs_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_vs_gps2_day0_bf.narrowPeak',
        withFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_noFCCP_bf.narrowPeak',
        noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_bf.narrowPeak',
        withFCCP_vs_noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_vsgps2_noFCCP_bf.narrowPeak',
        sictl = 'Silencing/GPS2/results/filtered/sictl_bf.narrowPeak',
        sigps2 = 'Silencing/GPS2/results/filtered/sigps2_bf.narrowPeak',
        sictl_vs_sigps2 = 'Silencing/GPS2/results/filtered/sictl_vs_sigps2_bf.narrowPeak',
        atf4_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/filtered/ATF4_AAD_bf.narrowPeak'
    output:
        day0_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_aad.bed',
        day6_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_aad.bed',
        day6_vs_day0_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_aad.bed',
        withFCCP_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_aad.bed',
        noFCCP_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_aad.bed',
        withFCCP_vs_noFCCP_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_aad.bed',
        sictl_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_aad.bed',
        sigps2_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_aad.bed',
        sictl_vs_sigps2_aad = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_aad.bed'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bedtools_env.yml'
    shell:
        '''
        bedtools intersect -a {input.gps2_day0} -b {input.atf4_aad} > {output.day0_aad}
        bedtools intersect -a {input.gps2_day6} -b {input.atf4_aad} > {output.day6_aad}
        bedtools intersect -a {input.day6_vs_day0} -b {input.atf4_aad} > {output.day6_vs_day0_aad}
        bedtools intersect -a {input.withFCCP} -b {input.atf4_aad} > {output.withFCCP_aad}
        bedtools intersect -a {input.noFCCP} -b {input.atf4_aad} > {output.noFCCP_aad}
        bedtools intersect -a {input.withFCCP_vs_noFCCP} -b {input.atf4_aad} > {output.withFCCP_vs_noFCCP_aad}
        bedtools intersect -a {input.sictl} -b {input.atf4_aad} > {output.sictl_aad}
        bedtools intersect -a {input.sigps2} -b {input.atf4_aad} > {output.sigps2_aad}
        bedtools intersect -a {input.sictl_vs_sigps2} -b {input.atf4_aad} > {output.sictl_vs_sigps2_aad}
        '''

rule overlaps_WT_TMC:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day0_bf.narrowPeak',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_bf.narrowPeak',
        day6_vs_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_vs_gps2_day0_bf.narrowPeak',
        withFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_noFCCP_bf.narrowPeak',
        noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_bf.narrowPeak',
        withFCCP_vs_noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_vsgps2_noFCCP_bf.narrowPeak',
        sictl = 'Silencing/GPS2/results/filtered/sictl_bf.narrowPeak',
        sigps2 = 'Silencing/GPS2/results/filtered/sigps2_bf.narrowPeak',
        sictl_vs_sigps2 = 'Silencing/GPS2/results/filtered/sictl_vs_sigps2_bf.narrowPeak',
        atf4_wt_tmc = 'Mito_stress_treatment/ER_stress/ATF4/results/filtered/ATF4_WT_Tmc_bf.narrowPeak'
    output:
        day0_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_wt_tmc.bed',
        day6_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_wt_tmc.bed',
        day6_vs_day0_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_wt_tmc.bed',
        withFCCP_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_wt_tmc.bed',
        noFCCP_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_wt_tmc.bed',
        withFCCP_vs_noFCCP_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_wt_tmc.bed',
        sictl_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_wt_tmc.bed',
        sigps2_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_wt_tmc.bed',
        sictl_vs_sigps2_wt_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_wt_tmc.bed'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bedtools_env.yml'
    shell:
        '''
        bedtools intersect -a {input.gps2_day0} -b {input.atf4_wt_tmc} > {output.day0_wt_tmc}
        bedtools intersect -a {input.gps2_day6} -b {input.atf4_wt_tmc} > {output.day6_wt_tmc}
        bedtools intersect -a {input.day6_vs_day0} -b {input.atf4_wt_tmc} > {output.day6_vs_day0_wt_tmc}
        bedtools intersect -a {input.withFCCP} -b {input.atf4_wt_tmc} > {output.withFCCP_wt_tmc}
        bedtools intersect -a {input.noFCCP} -b {input.atf4_wt_tmc} > {output.noFCCP_wt_tmc}
        bedtools intersect -a {input.withFCCP_vs_noFCCP} -b {input.atf4_wt_tmc} > {output.withFCCP_vs_noFCCP_wt_tmc}
        bedtools intersect -a {input.sictl} -b {input.atf4_wt_tmc} > {output.sictl_wt_tmc}
        bedtools intersect -a {input.sigps2} -b {input.atf4_wt_tmc} > {output.sigps2_wt_tmc}
        bedtools intersect -a {input.sictl_vs_sigps2} -b {input.atf4_wt_tmc} > {output.sictl_vs_sigps2_wt_tmc}
        '''

rule overlaps_KO_TMC:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day0_bf.narrowPeak',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_bf.narrowPeak',
        day6_vs_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_vs_gps2_day0_bf.narrowPeak',
        withFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_noFCCP_bf.narrowPeak',
        noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_bf.narrowPeak',
        withFCCP_vs_noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_vsgps2_noFCCP_bf.narrowPeak',
        sictl = 'Silencing/GPS2/results/filtered/sictl_bf.narrowPeak',
        sigps2 = 'Silencing/GPS2/results/filtered/sigps2_bf.narrowPeak',
        sictl_vs_sigps2 = 'Silencing/GPS2/results/filtered/sictl_vs_sigps2_bf.narrowPeak',
        atf4_ko_tmc = 'Mito_stress_treatment/ER_stress/ATF4/results/filtered/ATF4_KO_Tmc_bf.narrowPeak',
    output:
        day0_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_ko_tmc.bed',
        day6_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_ko_tmc.bed',
        day6_vs_day0_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_ko_tmc.bed',
        withFCCP_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_ko_tmc.bed',
        noFCCP_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_ko_tmc.bed',
        withFCCP_vs_noFCCP_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_ko_tmc.bed',
        sictl_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_ko_tmc.bed',
        sigps2_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_ko_tmc.bed',
        sictl_vs_sigps2_ko_tmc = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_ko_tmc.bed'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bedtools_env.yml'
    shell:
        '''
        bedtools intersect -a {input.gps2_day0} -b {input.atf4_ko_tmc} > {output.day0_ko_tmc}
        bedtools intersect -a {input.gps2_day6} -b {input.atf4_ko_tmc} > {output.day6_ko_tmc}
        bedtools intersect -a {input.day6_vs_day0} -b {input.atf4_ko_tmc} > {output.day6_vs_day0_ko_tmc}
        bedtools intersect -a {input.withFCCP} -b {input.atf4_ko_tmc} > {output.withFCCP_ko_tmc}
        bedtools intersect -a {input.noFCCP} -b {input.atf4_ko_tmc} > {output.noFCCP_ko_tmc}
        bedtools intersect -a {input.withFCCP_vs_noFCCP} -b {input.atf4_ko_tmc} > {output.withFCCP_vs_noFCCP_ko_tmc}
        bedtools intersect -a {input.sictl} -b {input.atf4_ko_tmc} > {output.sictl_ko_tmc}
        bedtools intersect -a {input.sigps2} -b {input.atf4_ko_tmc} > {output.sigps2_ko_tmc}
        bedtools intersect -a {input.sictl_vs_sigps2} -b {input.atf4_ko_tmc} > {output.sictl_vs_sigps2_ko_tmc}
        '''

rule overlaps_WT_vs_KO_TMC:
    input:
        gps2_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day0_bf.narrowPeak',
        gps2_day6 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_bf.narrowPeak',
        day6_vs_day0 = 'Adipocyte_differentiation/GPS2/results/filtered/gps2_day6_vs_gps2_day0_bf.narrowPeak',
        withFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_noFCCP_bf.narrowPeak',
        noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_bf.narrowPeak',
        withFCCP_vs_noFCCP = 'Mito_stress_treatment/GPS2/results/filtered/gps2_withFCCP_vsgps2_noFCCP_bf.narrowPeak',
        sictl = 'Silencing/GPS2/results/filtered/sictl_bf.narrowPeak',
        sigps2 = 'Silencing/GPS2/results/filtered/sigps2_bf.narrowPeak',
        sictl_vs_sigps2 = 'Silencing/GPS2/results/filtered/sictl_vs_sigps2_bf.narrowPeak',
        wt_vs_ko = 'Mito_stress_treatment/ER_stress/ATF4/results/filtered/ATF4_WT_Tmc_vs_ATF4_KO_Tmc_bf.narrowPeak'
    output:
        day0_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day0_over_atf4_wt_vs_ko.bed',
        day6_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_over_atf4_wt_vs_ko.bed',
        day6_vs_day0_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_day6_vs_gps2_day0_over_atf4_wt_vs_ko.bed',
        withFCCP_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_over_atf4_wt_vs_ko.bed',
        noFCCP_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_noFCCP_over_atf4_wt_vs_ko.bed',
        withFCCP_vs_noFCCP_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/gps2_withFCCP_vs_gps2_noFCCP_over_atf4_wt_vs_ko.bed',
        sictl_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_over_atf4_wt_vs_ko.bed',
        sigps2_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sigps2_over_atf4_wt_vs_ko.bed',
        sictl_vs_sigps2_wt_vs_ko = 'Mito_stress_treatment/Amino_acid_deprivation/ATF4/results/overlaps/sictl_vs_sigps2_over_atf4_wt_vs_ko.bed'
    conda:
        '/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/envs/bedtools_env.yml'
    shell:
        '''
        bedtools intersect -a {input.gps2_day0} -b {input.wt_vs_ko} > {output.day0_wt_vs_ko}
        bedtools intersect -a {input.gps2_day6} -b {input.wt_vs_ko} > {output.day6_wt_vs_ko}
        bedtools intersect -a {input.day6_vs_day0} -b {input.wt_vs_ko} > {output.day6_vs_day0_wt_vs_ko}
        bedtools intersect -a {input.withFCCP} -b {input.wt_vs_ko} > {output.withFCCP_wt_vs_ko}
        bedtools intersect -a {input.noFCCP} -b {input.wt_vs_ko} > {output.noFCCP_wt_vs_ko}
        bedtools intersect -a {input.withFCCP_vs_noFCCP} -b {input.wt_vs_ko} > {output.withFCCP_vs_noFCCP_wt_vs_ko}
        bedtools intersect -a {input.sictl} -b {input.wt_vs_ko} > {output.sictl_wt_vs_ko}
        bedtools intersect -a {input.sigps2} -b {input.wt_vs_ko} > {output.sigps2_wt_vs_ko}
        bedtools intersect -a {input.sictl_vs_sigps2} -b {input.wt_vs_ko} > {output.sictl_vs_sigps2_wt_vs_ko}
        '''
