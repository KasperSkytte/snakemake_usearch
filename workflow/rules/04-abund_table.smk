# MUUUCH faster to do this per sample and merge afterwards
# usearch -otutab does NOT scale linearly with more threads
# more than 16 has diminishing returns
rule abund_table:
    input:
        zotus=os.path.join(config["output_dir"], "zOTUs.fa"),
        allreads_unfiltered=os.path.join(
            config["tmp_dir"], "01-sample_prep", "{sample}", "{sample}_renamed.fastq"
        ),
    output:
        temp(os.path.join(config["tmp_dir"], "abund_table", "abund_table_{sample}.tsv")),
    log:
        os.path.join(config["log_dir"], "05-abund_table", "abund_table_{sample}.log"),
    message:
        "{wildcards.sample}: Estimating abundances of zOTUs/ASVs"
    resources:
        mem_mb=1024,  # this needs to be calculated dynamically based on input file size
        runtime=300,
        cpus_per_task=lambda wc, input: min(config["max_threads"], 16),
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: lambda wc, input: min(config["max_threads"], 16)
    params:
        sample_sep=config["sample_sep"],
    shell:
        """
            exec &> "{log}"
            set -euxo pipefail
        
            usearch -otutab \
                "{input.allreads_unfiltered}" \
                -zotus "{input.zotus}" \
                -otutabout "{output}" \
                -threads "{threads}" \
                -sample_delim "{params.sample_sep}"
        """


rule merge_abund_tables:
    input:
        expand(
            os.path.join(config["tmp_dir"], "abund_table", "abund_table_{sample}.tsv"),
            sample=sample_dirs,
        ),
    output:
        os.path.join(config["output_dir"], "abund_table.tsv"),
    log:
        os.path.join(config["log_dir"], "05-abund_table", "merge_abund_tables.log"),
    message:
        "Merging abundance tables"
    resources:
        mem_mb=2048,
        runtime=30,
        cpus_per_task=1,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 1
    params:
        input_csv=lambda wildcards, input: ",".join(input),
    shell:
        """
            exec &> "{log}"
            set -euxo pipefail

            # filter empty tables from the list to avoid usearch fails
            non_empty_tables=""
            for table in {input}; do
                if [ "$(wc -l < "$table")" -gt 1 ]; then
                    if [ -z "$non_empty_tables" ]; then
                        non_empty_tables="$table"
                    else
                        non_empty_tables="$non_empty_tables,$table"
                    fi
                fi
            done

            # check if we actually have anything left to merge
            if [ -z "$non_empty_tables" ]; then
                echo "Error: All abundance tables are empty"
                exit 1
            fi

            usearch -otutab_merge "$non_empty_tables" -output "{output}"
        """


rule rarefy_abund_table:
    input:
        os.path.join(config["output_dir"], "abund_table.tsv"),
    output:
        os.path.join(config["output_dir"], "abund_table_rarefied.tsv"),
    log:
        os.path.join(config["log_dir"], "05-abund_table", "rarefy_abund_tables.log"),
    message:
        "Rarefying abundance table"
    resources:
        mem_mb=2048,
        runtime=120,
        cpus_per_task=1,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 1
    params:
        rarefy_sample_size=config["rarefy_sample_size"],
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail

        usearch -otutab_rare {input} -sample_size {params.rarefy_sample_size} -output {output}
        """
