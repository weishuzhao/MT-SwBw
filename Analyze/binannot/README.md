<!--
 * @Date: 2021-11-16 13:34:23
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-11-22 11:26:32
 * @FilePath: /2021_09-MT10kSW/Analyze/binannot/README.md
 * @Description:
-->
Map general gene annotation to bins
===

---
## Map marker to Bins
- output format:
    - for each sample: generate a table:
        - BinId, KO/marker, RPb
- This is generated from:
    - BinId from drep
    - BinID -- contigId from DASTool
    - contig -- KO -- gene -- RPb from RPb
    - marker -- gene from checkm marker

## Map entry to pathway
```python
import pandas as pd
sgf = pd.read_csv("Analyze/binannot/sgf.tsv", sep="\t")
entries = [entry
           for paths in sgf[sgf.columns[1:]].values
           for path in paths if not pd.isna(path)
           for entry in path.split(", ")]

module_name = pd.read_csv("Analyze/pathway/module_name.tsv", sep="\t")
entries
names = [module_name[module_name["entry"] == entry
                     ]["name"].values[0]
         for entry in entries]
PATHs = list({
    PATH
    for name in names
    for PATH in (
        name[name.index("[PATH:") + 6: -1
             ] if "[PATH:" in name else ""
    ).split(" ") if PATH})
entry_path = pd.DataFrame(
    [[int(PATH in name) for PATH in PATHs] for name in names],
    index = entries, columns = PATHs)
entry_path.to_csv("tmp.tsv", sep="\t")
```


# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)
