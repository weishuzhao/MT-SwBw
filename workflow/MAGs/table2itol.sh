#!/bin/bash
:<<!EOF!
 * @Date: 2022-06-27 21:34:40
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-08-27 15:24:05
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/table2itol.sh
 * @Description:
!EOF!
set -e && echo "$0 $*" >&2 && source ~/.miniconda_init
conda activate R4.1

set -vx
Rscript ~/software/table2itol.R \
    -a  `#Abort if a requested column cannot be found` \
    -d  `#Create bar charts, not gradients, from double` \
    -D 04_bin/04_tree/table2itol  `#output files in this directory` \
    \
    -b Color  `#Column to define the background colours` \
    -i Genome  `#Mandatory identifier column` \
    -l Label  `#Column to define the tip labels displayed in the picture` \
    -t "%s" \
    -s "," \
    04_bin/04_tree/table2itol.csv


# LEGEND_COLORS	#ffffff	#984ea3	#d7191c	#2b83ba
# LEGEND_LABELS		cross	sediment	water

# #20854E Completeness
# #E18727 Contamination
# #0072B5 GenomeSize
# #BC3C29 GC

grep \
    "LEGEND_" \
    04_bin/04_tree/table2itol/iTOL_treecolors-Color.txt \
|grep -v "_SHAPES" |grep -v "_TITLE" \
|awk \
    'BEGIN{FS=" "} {for(i=1;i<=NF;i++)a[NR,i]=$i}
     END{for(j=1;j<=NF;j++)
         {for(k=1;k<=NR;k++) {printf a[k,j]"\t"}
         printf "\n"}}' \
> data/taxon_color.tsv
