<!--
 * @Date: 2021-09-15 20:57:09
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-02 16:32:35
 * @FilePath: /2021_09-MT10kSW/README.md
 * @Description:
-->
Microbial communities reveal niche partitioning across the slope and bottom zones of the challenger deep
===

- Analysis for [Microbial communities reveal niche partitioning across the slope and bottom zones of the challenger deep](https://doi.org/10.1111/1758-2229.13314)

---
## Description
- This reposity contain files generate from metagenomic data and tables for figures and other conclusions in the manuscript.

## Run
```bash
# conda and snakemake required
python -m snakemake all --use-conda
```

- files are kept in `data_table.zip` and can be unzipped before use it
    with suffixes ([xaa](data_table.xaa) and [xab](data_table.xab)).

```bash
cat data_table.xa* > data_table.zip
unzip data_table.zip
```

- information of software and version are shown in manuscript.


# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)
