<!--
 * @Date: 2022-09-11 16:17:53
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-11 17:32:13
 * @FilePath: /2021_09-MT10kSW/docs/Results.1.md
 * @Description:
-->

===
## Results

### Description of samples and environmental parameters
-
    Three dives of the lander “Tianya”, four dives of the lander “Yuanwei” and three dives of the lander “Wanquan” collected water samples at ten stations from the Challenger Deep in the Mariana Trench during “TS-09” cruise on aboard the R/V “Tansuo-01”, varying with water depth ranging from 7-11 km b.s.l under hydrostatic pressure from 75 to 115 MPa and covering the bottom-axis (seven stations) and southern slope (three stations) of the Challenger Deep (Fig. 1, Methods).
    The temperature of these sampling sites ranged from 1.8 to 2.6 °C and varies with the water depth (Table 1),
        while the salinity was around 34.7 psu (practical salinity unit), ~3.5% (w/v).
    We further measured the concentration of nitrogen species in seawater samples,
        and compared with the sediment samples from the similar stations reported in previous studies
        [Ref: 2022, NC, Microbiomes in the Challenger Deep slope and bottom-axis sediments; Ref: Microbial community and geochemical analyses of trans-trench sediments for understanding the roles of hadal environments. ISME J 14, 740–756 (2020); Ref: Organic matter diagenesis in hadal setting: Insights from the pore-water geochemistry of the Mariana Trench sediments. Deep Sea Res I 147, 22-31 (2019)].
    Levels of the concentration of ammonium and nitrite in seawater samples
        was comparable to those in sediment samples at the similar water depths in the bottom-axis of the Challenger Deep,
        where ammonium was measured as 1.9-3.9 μM in seawater samples and reported as 0.8-6.3 μM in the top layer of sediment samples (Additional file: Table S1).
    However, the levels of nitrate concentration in seawater samples
        were two orders of magnitude lower than those in sediment samples at the similar stations,
        which was measured as 0.12-0.14 μM in seawater but 14 -33 μM of the top surface of sediments (Additional file: Table S1),
        which indicated a potential large difference of nitrogen cycling capabilities
        in the microbiome living in the seawater and sediments in the Challenger Deep.

### [Metagenomic profile](../results/Rmd/check_results.Rmd)
-
    The shot-gun sequencing of ten near-benthic seawater samples
        yielded 59,679,900 clean raw read pairs after being trimmed and filtered, composing 98.36% of the raw reads
        (Methods).
    In addition, we collected
        23 metagenomes at seven stations of sediment samples from the bottom-axis (four stations) and southern slope (three stations) of the Challenger Deep as comparisons
        (Fig. 1, Additional file: Table S1) [Ref: 2022, NC, Microbiomes in the Challenger Deep slope and bottom-axis sediments].
    Assembly of the above 33 metagenomes generated 11,162,542 scaffolds (average length 1127 bp).
    Binning of these contigs resulted in 792 draft genomes with qualities above the middle-quality threshold (≥ 50% completeness and ≤ 10% contamination),
        among which 289 were above the high-quality threshold (≥ 90% completeness and ≤ 5% contamination) according to a widely accepted standard of MAGs
        [Ref: Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea].

### Comparison of prokaryotic community between seawater and sediments at bottom and slope
-
    A total of 831,935 reads were mapped to 5,556 Molecular Operational Taxonomic Unit (mOTU) sequences at the species level or higher based on the MiTAG method using phyloFlash (v3.4)
        (Methods, Additional file 1: Fig. S1, Additional file 2: Table S1).
    92.21% mOTUs can be annotated to class level in all samples
        (Additional file 1: Fig. S1),
        therefore, we further use class-level taxa to describe the community compositions among four groups
        which are slope seawater (Sw), slope sediments (Ss), bottom seawater (Bs) and bottom sediments (Bs).
    Considering the existence of all taxa at class level among above four groups, the PCoA show a separated clustering between the seawater and sediments samples (Binary Jaccard distance, Adonis test, P value < 0.001)
        (Fig. 1A).
    However, when considering the abundance of each class-level taxa in the community among samples of above four groups, the cluster of slope seawater is much closer to the slope sediments, but significantly separated from the cluster of bottom seawater (Bray-Curtis distance, Adonis test, P value < 0.001)
        (Fig. 1B).

-
    The community composition of prokaryotes of the bottom seawater is the most distinct from the other three groups.
        (Fig. 1B-C, Additional file 1: Fig. S2)
    Gammaproteobacteria, Alphaproteobacteria, Actinomycetia, and Bacteroidota are dominated accounting for more than 90% of the total in bottom seawater, where Gammaproteobacteria account for more than 50% in all seven samples.
    Interestingly, the propotion of Nitrososphaeria (ammonia oxidizing archaea, formerly named as Thaumarchaeota) [Ref] is quiet low among all samples in the bottom seawater which is less than 1%.
    In contrast, the composition of prokaryotic community indicated that relative abundance of main taxa at class level of the slope seawater is similar to the slope sediments
        (Fig. 1C).
    The percentage of Gammaproteobacteria in the slope seawater is significantly lower than those in bottom seawater, while significant increasing of Nitrososphaeria and Marinisomatota were observed in slope seawater than the bottom seawater
        (Additional file 1: Fig. S2).
    The different percentage of the ammonia oxidizing archaea, Nitrososphaeria, suggested a potentially different capacity of nitrogen cycle in the seawater at the slope and the bottom.

-
    The similar trend that more common taxa across the sediments and seawater samples at the slope than the bottom were also observed at metagenomic assembled genomes (MAGs).
    We constructed a phylogenetic tree by using 176 representative MAGs (dereplicated at 95% average nucleotide identity (ANI)) at slope and 323 representative MAGs at bottom together with 161 reference genomes from the GTDB
        [13] (Methods, Fig. 3).
    These MAGs were distributed into 25 phyla at slope while 28 phyla at bottom.
        (Additional file 1: Fig. S2)

### Metabolic versatile and distinct abundancy between seawater and sediments in slope and bottom
