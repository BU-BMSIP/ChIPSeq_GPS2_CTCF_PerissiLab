import pandas as pd
import os

rule all:
    input:
        expand("Silencing/GPS2/results/matrix/GPS2_signal_on_CTCF_{tp}_promoters.gz", tp=["t1", "t2", "t3", "t4"]),
        expand("Silencing/GPS2/results/plots/GPS2_signal_on_CTCF_{tp}_promoters_profile.png", tp=["t1", "t2", "t3", "t4"]),
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/plots/ctcf_gps2_signal_on_common_promoters.png"

rule compute_matrix_GPS2_signal_on_CTCF_promoters:
    input:
        bw = "Silencing/GPS2/results/bigwig/chip_sictl.bw",
        bed = lambda wildcards: f"CTCF_3T3L1/results/annotation/CTCF_{wildcards.tp}_promoter_only.bed"
    output:
        matrix = "Silencing/GPS2/results/matrix/GPS2_signal_on_CTCF_{tp}_promoters.gz"
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

rule plot_GPS2_signal_on_CTCF_promoters:
    input:
        matrix = "Silencing/GPS2/results/matrix/GPS2_signal_on_CTCF_{tp}_promoters.gz"
    output:
        plot = "Silencing/GPS2/results/plots/GPS2_signal_on_CTCF_{tp}_promoters_profile.png"
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} \
            -out {output.plot} \
            --legendLocation upper-right \
            --plotTitle "GPS2 signal at CTCF {wildcards.tp} promoter peaks"
        """

rule compute_matrix_signal_on_common_promoters:
    input:
        bed= "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_CTCF_common_gene_peaks_centered_±3kb.bed",
        bw = [ 
            "CTCF_3T3L1/results/bigwig/CTCF_t1.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t2.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t3.bw",
            "CTCF_3T3L1/results/bigwig/CTCF_t4.bw",
            "Silencing/GPS2/results/bigwig/chip_sictl.bw"
        ]
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    output:
        matrix = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/matrix/ctcf_gps2_signal_on_common_promoters.gz",
    threads: 8
    shell:
        """
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.bed} \
            --referencePoint center \
            --beforeRegionStartLength 3000 \
            --afterRegionStartLength 3000 \
            --binSize 50 \
            --skipZeros \
            -o {output.matrix}
        """

rule plot_signal_on_common_promoters:
    input:
        matrix = rules.compute_matrix_signal_on_common_promoters.output.matrix
    output:
        plot = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/plots/ctcf_gps2_signal_on_common_promoters.png"
    threads: 8
    conda:
        "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/envs/deeptools_env.yml"
    shell:
        """
        plotProfile \
            -m {input.matrix} \
            -out {output.plot} \
            --perGroup \
            --plotTitle "GPS2 vs CTCF signal on common promoters"
        """
