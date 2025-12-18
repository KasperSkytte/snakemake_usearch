rule concat_all:
    input:
        expand(
            os.path.join(
                config["tmp_dir"],
                "01-sample_prep",
                "{sample}",
                "{sample}_filtered_renamed.fastq",
            ),
            sample=sample_dirs,
        ),
    output:
        temp(
            os.path.join(
                config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed.fastq"
            )
        ),
    log:
        os.path.join(config["log_dir"], "02-denoise", "concat_all.log"),
    message:
        "Concatenating all samples before denoising (generate zOTUs/ASVs)"
    resources:
        mem_mb=1024,
        runtime=30,
        cpus_per_task=1,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
        
        cat {input} > {output}

        # error if output file is empty
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
      """


rule cutadapt:
    input:
        os.path.join(
            config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed.fastq"
        ),
    output:
        temp(
            os.path.join(
                config["tmp_dir"],
                "02-denoise",
                "all_samples_filtered_renamed_oriented_trimmed.fastq",
            )
        ),
    log:
        os.path.join(config["log_dir"], "02-denoise", "cutadapt.log"),
    params:
        cutadapt_args=config["cutadapt_args"],
    message:
        "Orienting and trimming reads using cutadapt"
    resources:
        mem_mb=2048,
        runtime=120,
        cpus_per_task=config["max_threads"],
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: config["max_threads"]
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail

        # this step both orients and trims at once
        cutadapt {params.cutadapt_args} -o {output} {input} -j {threads}

        # error if output file is empty
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
        """


# according to Edgar, QC filtering should be done here, not in sample_prep
# rule qc_filter:
#   shell:
#     """
#     usearch -fastq_filter trimmed.fq -fastq_maxee 1.0 -fastaout filtered.fa
#     """


rule derep:
    input:
        os.path.join(
            config["tmp_dir"],
            "02-denoise",
            "all_samples_filtered_renamed_oriented_trimmed.fastq",
        ),
    output:
        temp(
            os.path.join(
                config["tmp_dir"],
                "02-denoise",
                "all_samples_filtered_renamed_oriented_trimmed_derep.fa",
            )
        ),
    log:
        os.path.join(config["log_dir"], "02-denoise", "derep.log"),
    message:
        "Dereplicating all samples"
    resources:
        mem_mb=lambda wc, input: max(1.5 * input.size_mb, 4096),  # a 64GB file took 72GB mem
        runtime=600,
        cpus_per_task=4,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    params:
        derep_minsize=config["derep_minsize"],
    threads: 4  # this command spends most of the time just reading in the file, increasing isn't beneficial
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
        
        usearch -fastx_uniques \
          {input} \
          -fastaout {output} \
          -sizeout \
          -threads {threads} \
          -minuniquesize {params.derep_minsize} \
          -relabel Uniq

        # error if output file is empty
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
      """


rule unoise:
    input:
        os.path.join(
            config["tmp_dir"],
            "02-denoise",
            "all_samples_filtered_renamed_oriented_trimmed_derep.fa",
        ),
    output:
        os.path.join(config["output_dir"], "zOTUs.fa"),
    log:
        os.path.join(config["log_dir"], "02-denoise", "unoise.log"),
    message:
        "Denoising/generating zOTUs/ASVs"
    resources:
        mem_mb=lambda wc, input: max(3 * input.size_mb, 512),
        runtime=600,
        cpus_per_task=1,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    params:
        unoise_minsize=config["unoise_minsize"],
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
        
        usearch -unoise3 \
          {input} \
          -zotus {output} \
          -sizeout \
          -minsize {params.unoise_minsize}

        # error if output file is empty
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
    """
