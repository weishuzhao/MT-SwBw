"""
 * @Date: 2022-06-21 11:16:40
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-06 22:45:08
 * @FilePath: /2021_09-MT10kSW/workflow/remote/drep.smk
 * @Description:
"""

rule drep_create:
    input:
        bin_dirs = [file_path.dastool_bins(f"{site}") for site in site_layer_dict],
        checkms  = [file_path.checkm_of(file_path.dastool_bins(f"{site}")) for site in site_layer_dict]
    output:
        drep_out = directory(file_path.all_bins("drep")),
    threads: THREADS-1
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        mkdir -p {output.drep_out}

        for i in {input.bin_dirs}
        do
            ls $i/*.fa >> {output.drep_out}/all_bins.txt
        done
        echo "genome,completeness,contamination,strain_heterogeneity" > {output.drep_out}/checkm.tsv
        for i in {input.checkms}
        do
            awk -v FS='\t' -v OFS=',' '! /^Bin/ {{print $1".fa",$12,$13,$14}}' \
                $i \
            >> {output.drep_out}/checkm.tsv
        done

        dRep dereplicate \
            {output.drep_out} \
            -p {threads} \
            -comp 50 -con 10 \
            -pa 0.9 -sa 0.95 \
            -N50W 2 \
            -g {output.drep_out}/all_bins.txt \
            --genomeInfo {output.drep_out}/checkm.tsv \
            --debug
        """


rule gtdbtk_create:
    input:
        drep_out = file_path.all_bins("drep"),
    output:
        taxon_out = directory(file_path.all_bins("gtdbtk")),
    threads: THREADS-1
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate gtdbtk
        GTDBTK_DATA_PATH=~/Data/Database2/GTDB/release202/

        mkdir scratch_dir

        gtdbtk classify_wf \
            --cpus {threads} \
            --genome_dir {input.drep_out}/dereplicated_genomes -x fa \
            --out_dir {output.taxon_out} \
            --scratch_dir scratch_dir
        """
