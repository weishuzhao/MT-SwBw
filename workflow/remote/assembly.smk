contig = str(file_path.contig("{site}", "_cut"))
contig_raw = str(file_path.contig("{site}", "..megahit"))
jgi = str(file_path.contig("{site}", "_cut", "-jgi.depth"))

MIN_CONTIG_LEN = 500


rule assembly_megahit:
    input:
        r1s = lambda _: [file_path.trimmed_reads(_.site, layer, 1) for layer in site_layer_dict[_.site]],
        r2s = lambda _: [file_path.trimmed_reads(_.site, layer, 2) for layer in site_layer_dict[_.site]],
    output:
        contig = protected(contig_raw),
    params:
        outdir = file_path.contig("{site}", "..megahit", "-dir"),
    log:
        file_path.log("02_assem_megahit", "{site}"),
    threads: THREADS
    shell:
        """
        mkdir -p pipe/{wildcards.site}

            set +u
            source workflow/utils/.conda_init
            conda activate python36

        r1s=`printf "{input.r1s}" | tr "$IFS" ","`
        r2s=`printf "{input.r2s}" | tr "$IFS" ","`

        megahit \
            -1 $r1s -2 $r2s \
            -o {params.outdir} \
            -t {threads} \
        2>&1 |tee {log}
        cp {params.outdir}/final.contigs.fa {output.contig}
        """


rule assembly_rename_cut:
    input:
        contig = contig_raw
    output:
        contig = contig
    params:
        site = "{site}"
    run:
        from utils.contig_rename_cut import rename_cut
        with open(input.contig) as fi, open(output.contig, "w") as fo:
            rename_cut(fi, fo, MIN_CONTIG_LEN, params.site)


rule map_bbmap:
    input:
        r1 = str(file_path.trimmed_reads("{site}", "{layer}", 1)),
        r2 = str(file_path.trimmed_reads("{site}", "{layer}", 2)),
        contig = contig,
    output:
        bbdepth = file_path.bam("{site}", "{layer}", "-bbmap.depth"),
        bam     = protected(file_path.bam("{site}", "{layer}")),
        bai     = protected(file_path.bam("{site}", "{layer}", ".bam.bai")),
    threads: THREADS
    log:
        file_path.log("02_assem_map_bbmap", "{site}", "{layer}"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python36

        bbmap.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.contig} \
            nodisk k=13 minid=0.95 keepnames=t \
            covstats={output.bbdepth} \
            minaveragequality=5 \
            out={output.bam}.sam \
            trimreaddescriptions=t pairlen=350 rescuedist=650 \
        2>&1 |tee {log}

        samtools view \
            -@ {threads} -b {output.bam}.sam \
        | samtools sort \
            -@ {threads} -o {output.bam}

        samtools index {output.bam} -@ {threads}
        rm {output.bam}.sam

        """


rule depth_jgi:
    input:
        bams = lambda _: [file_path.bam(_.site, layer) for layer in site_layer_dict[_.site]],
    output:
        jgi = jgi
    log:
        file_path.log("02_assem_depth_jgi", "{site}"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python36

        jgi_summarize_bam_contig_depths \
            --outputDepth {output.jgi} \
            {input.bams} \
        2>&1 |tee {log}
        """


rule plass:
    input:
        r1 = "pipe/{sample}/01_trim..{sample}_1.fq.gz",
        r2 = "pipe/{sample}/01_trim..{sample}_2.fq.gz",
    output:
        out = "pipe/{sample}/02_assem..{sample}.prot.fa",
    threads: THREADS-1
    shell:
        """
        source ~/.conda_init
        conda activate python39

        plass assemble \
            {input.r1} {input.r2} \
            {output.out} \
            {output.out}-tmp \
            --threads {threads} \
            -v
        """
