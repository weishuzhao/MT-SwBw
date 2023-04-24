r1 = str(file_path.trimmed_reads("{site}", "{layer}", 1))
r2 = str(file_path.trimmed_reads("{site}", "{layer}", 2))
drop1 = str(file_path.trimmed_reads("{site}", "{layer}", "drop.1"))


rule trim_sickle:
    #input:
    output:
        r1 = r1,
        r2 = r2,
        drop1 = drop1,
    params:
        ini_raw = "data/{site}..{layer}"
    threads: THREADS
    log:
        file_path.log("01_trim_sickle", "{site}", "{layer}"),
    shell:
        """
        mkdir -p pipe/{wildcards.site}

            set +u
            source workflow/utils/.conda_init
            conda activate python36

        echo "#"$(date +%F%n%T)
        # length-threshold = 60% * 150bp
        sickle pe \
            -f {params.ini_raw}_1.fq.gz -r {params.ini_raw}_2.fq.gz \
            -o {output.r1} -p {output.r2} \
            -s {output.drop1} \
            --length-threshold 90 \
            -t sanger -g \
        2>&1 |tee {log}

        """


rule unzip_reads:
    input:
        r = file_path.trimmed_reads("{site}", "{layer}", "{any}"),
    output:
        r = temp(str(file_path.trimmed_reads("{site}", "{layer}", "{any}")).rsplit(".gz", 1)[0]),
    shell:
        """
        zcat {input.r} > {output.r}
        """


rule nonpareil:
    input:
        r1 = r1.rsplit(".gz", 1)[0],
    output:
        nonpareil = directory(file_path.stat_reads("{site}", "{layer}", "nonpareil", "/")),
        npo = protected(file_path.stat_reads("{site}", "{layer}", "nonpareil", ".npo")),
    threads: THREADS
    shell:
        """
        mkdir -p {output.nonpareil}

        nonpareil -s {input.r1} -T kmer -f fastq -b {output.nonpareil}/ -t {threads}
        mv {output.nonpareil}/.npo {output.npo}
        """

rule phyloflash:
    input:
        r1 = r1.rsplit(".gz", 1)[0],
        r2 = r2.rsplit(".gz", 1)[0],
    output:
        phyloFlash_gz     = file_path.stat_reads("{site}", "{layer}", "phyloFlash", ".tar.gz"),
        NTUfull_abundance = file_path.stat_reads("{site}", "{layer}", "phyloFlash", ".csv"),
    threads: THREADS
    shadow: "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate py27

        declare PHYLOFLASH_DB=~/Data/Database2/phyloFlash/138.1

        #cd `dirname {output.NTUfull_abundance}`

        phyloFlash.pl \
            -dbhome $PHYLOFLASH_DB \
            -lib phyloFlash -almosteverything -log \
            -CPUs {threads} \
            -read1 {input.r1} -read2 {input.r2}

        tar -xzvf phyloFlash.phyloFlash.tar.gz \
                phyloFlash.phyloFlash.NTUfull_abundance.csv
        mv phyloFlash.phyloFlash.tar.gz \
            {output.phyloFlash_gz}
        mv phyloFlash.phyloFlash.NTUfull_abundance.csv \
            {output.NTUfull_abundance}

        #cd -

        """
