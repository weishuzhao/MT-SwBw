"""
 * @Date: 2022-06-18 16:43:26
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-08-24 17:42:54
 * @FilePath: /2021_09-MT10kSW/workflow/remote/binning.smk
 * @Description:
"""
methods = [
    *(f"metabat2_{p}_{s}" for p in (60, 75, 90) for s in (60, 75, 90)),
    *(f"maxbin2_{m}" for m in (107, 40)),
    "concoct"
]
MIN_BIN_CONTIG_LEN=1500

contig = str(file_path.contig("{site}", "_cut"))
faa = file_path.gene("{site}", suffix=".faa")
jgi = str(file_path.contig("{site}", "_cut", "-jgi.depth"))


rule DASTool_create:
    input:
        contig  = contig,
        faa     = faa,
        ctg2mag = [file_path.ctg2mag("{site}", method) for method in methods]
    output:
        log     = file_path.dastool("{site}", ".log"),
        scf2bin = file_path.dastool("{site}", "_scaffolds2bin.tsv"),
        summary = file_path.dastool("{site}", "_summary.tsv"),
        bin_dir = directory(file_path.dastool_bins("{site}")),
    params:
        ctg2mag = ",".join([str(file_path.ctg2mag("{site}", method)) for method in methods]),
        methods = ",".join(methods),
    log:
        file_path.log("04_bin_dastool", "{site}"),
    threads: THREADS
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        DAS_Tool \
            -i {params.ctg2mag} -l {params.methods} \
            -c {input.contig} --proteins {input.faa} \
            -o DASTool \
            --write_bins 1 --search_engine diamond --score_threshold 0 \
            -t {threads} \
            --debug

        cp DASTool_DASTool.log {output.log}
        cp DASTool_DASTool_scaffolds2bin.txt {output.scf2bin}
        cp DASTool_DASTool_summary.txt {output.summary}
        mv DASTool_DASTool_bins {output.bin_dir}
        """


rule checkm:
    input:
        bin_dir = file_path.dastool_bins("{site}"),
    output:
        checkm = file_path.checkm_of(file_path.dastool_bins("{site}")),
    threads: THREADS-1
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        checkm lineage_wf \
            -x fa -t {threads} \
            {input.bin_dir} tmp \
            --tab_table \
            -f {output.checkm}

        """


rule metabat2:
    input:
        contig  = contig,
        jgi     = jgi,
    output:
        ctg2mag = file_path.ctg2mag("{site}", "metabat2_{p}_{s}")
    params:
        dout = "{site}-metabat2_{p}_{s}"
    log:
        file_path.log("04_bin_dastool/" + "metabat2_{p}_{s}", "{site}"),
    threads: 1
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python36

        # mkdir -p $dout  # DONOT create the directory

        metabat2 \
            -i {input.contig} \
            -a {input.jgi} \
            -o {params.dout}/{params.dout} \
            -t {threads} \
            -v \
            --minContig {MIN_BIN_CONTIG_LEN} \
            --minS {wildcards.p} --maxP {wildcards.s}

        ~/software/anaconda3/envs/python39/bin/Fasta_to_Scaffolds2Bin.sh \
            -i {params.dout} -e fa \
        > {output.ctg2mag}
        """


rule maxbin2:
    input:
        contig  = contig,
        jgi     = jgi,
    output:
        ctg2mag = file_path.ctg2mag("{site}", "maxbin2_{m}")
    params:
        dout = "{site}-maxbin2_{m}"
    log:
        file_path.log("04_bin_dastool/" + "maxbin2_{m}", "{site}"),
    threads: 1
    priority: 1
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python36

        dout=`basename {output.ctg2mag}`
        mkdir -p {params.dout}

        run_MaxBin.pl \
            -min_contig_length {MIN_BIN_CONTIG_LEN} \
            -contig {input.contig} \
            -abund {input.jgi} \
            -out {params.dout}/{params.dout} \
            -markerset {wildcards.m} -thread {threads}

        ~/software/anaconda3/envs/python39/bin/Fasta_to_Scaffolds2Bin.sh \
            -i {params.dout} -e fasta \
        > {output.ctg2mag}
        """


rule concoct_contig:
    input:
        contig  = contig,
    output:
        contig  = temp(f"{contig}_concoct_tmp"),
    message:
        "concoct cannot filter short contigs itself"
    run:
        from Bio import SeqIO
        SeqIO.write(
            (i for i in SeqIO.parse(input.contig, "fasta") if len(i.seq) > MIN_BIN_CONTIG_LEN),
            output.contig,
            format="fasta",
        )


rule concoct:
    input:
        contig  = f"{contig}_concoct_tmp",
        bams = lambda _: [file_path.bam(_.site, layer) for layer in site_layer_dict[_.site]],
    output:
        ctg2mag = file_path.ctg2mag("{site}", "concoct")
    params:
        dout = "{site}-concoct"
    log:
        file_path.log("04_bin_dastool/" + "concoct", "{site}"),
    threads: THREADS - 2  # for maxbin
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python36

        dout=`basename {output.ctg2mag}`
        mkdir -p {params.dout}

        cut_up_fasta.py \
            {input.contig} -c 10000 -o 0 \
            --merge_last -b {params.dout}.bed \
            > {params.dout}.fasta
        concoct_coverage_table.py \
            {params.dout}.bed {input.bams} \
            > {params.dout}.tsv
        concoct --coverage_file {params.dout}.tsv \
            --composition_file {params.dout}.fasta \
            -t {threads} -r 150 -b {params.dout}/ -s 599 --no_original_data
        merge_cutup_clustering.py \
            {params.dout}/clustering_gt1000.csv \
            > {params.dout}/clustering_gt1000_merge.csv

        awk -v FS="," -v OFS="\t" \
            '{{if (NR==1) {{next}}; {{print $1,"{params.dout}_"$2}}}}' \
            {params.dout}/clustering_gt1000_merge.csv \
        | sort | awk '{{print $2,$1}}' \
        | sort | awk -v OFS="\t" '{{print $2,$1}}' \
        > {output.ctg2mag}
        """


rule init_binny_config:
    input:
        contig  = contig,
        jgi     = jgi,
        bams    = lambda _: [file_path.bam(_.site, layer) for layer in site_layer_dict[_.site]],
    output:
        config  = temp(str(file_path.ctg2mag("{site}", "binny")).rsplit(".tsv", 1)[0] + ".yaml")
    params:
        dout = "{site}-concoct"  # unused
    run:
        from pathlib import Path
        with open(input.jgi) as fi:
            line = fi.readline()
        bam = Path(input.jgi).parent / line.split()[3]
        config_str = f"""
        mem:
          big_mem_avail: FALSE
          big_mem_per_core_gb: 40
          normal_mem_per_core_gb: 4
        tmp_dir: tmp
        raws:
          assembly:  "{input.contig}"
          metagenomics_alignment:  "{input.bams}"
          contig_depth: "{input.jgi}"
        sample: "{wildcards.site}"
        outputdir: "{output.config}-dir"
        db_path: "$HOME/software/binny/database"
        binning:
          binny:
            kmers: '2,3,4'
            cutoff: {MIN_BIN_CONTIG_LEN}
            cutoff_marker: 0
            max_n_contigs: 3.5e5
            distance_metric: 'manhattan'
            embedding:
              max_iterations: 50
              tsne_early_exag_iterations: 250
              tsne_main_iterations: 750
            clustering:
              hdbscan_epsilon: 0.25
              hdbscan_min_samples: 2
              include_depth_initial: 'False'
              include_depth_main: 'True'
            bin_quality:
              completeness: 80
              purity: 85
        """[1:]
        with open(output.config, 'w') as fout:
            for line in config_str.split("\n"):
                print(line[8:], file=fout)


rule binny:
    input:
        config  = str(file_path.ctg2mag("{site}", "binny")).rsplit(".tsv", 1)[0] + ".yaml"
    output:
        ctg2mag = file_path.ctg2mag("{site}", "binny")
    log:
        file_path.log("04_bin_dastool/" + "binny", "{site}"),
    shadow: "shallow"
    threads: THREADS - 2  # for maxbin
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate smk
            set -vx

        ~/software/binny/binny -l -t {threads} -n {wildcards.site} -r {input.config}

        ~/software/anaconda3/envs/python39/bin/Fasta_to_Scaffolds2Bin.sh \
            -i {input.config}-dir/bins -e fasta \
        > {output.ctg2mag}
        """
