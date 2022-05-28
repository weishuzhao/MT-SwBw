<!--
 * @Date: 2021-09-15 20:57:09
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-09-15 21:14:11
 * @FilePath: /2021_09-MT10kSW/README.md
 * @Description:
-->

===

---
## log
```bash
mkdir -p Analyze/alphadiv


declare KAIJU_DB=~/Data/Database/Kaiju/nr_euk_20200525/

declare samples=(`cat 00_data/sample_meta.tsv |awk '! /^#/ {print $1}'`)
kaiju_outs=""
for i in ${samples[@]}
do
    kaiju_outs+=Pipe/${i}-megahit/01_alpha_div/kaiju-${i}-megahit.out" "
done

conda activate python36
kaiju2table \
    -t ${KAIJU_DB}/nodes.dmp \
    -n ${KAIJU_DB}/names.dmp \
    -o Analyze/alphadiv/kaiju.order.tsv -r order \
    -l superkingdom,phylum,class,order \
    ${kaiju_outs} \
    -v 2>&1 |tee Oerr/01.3_order-kaiju2table.log

```


# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)
