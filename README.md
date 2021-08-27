# UK_barley_popgen
# This is a tutorial for UK collegues to run Hierfstats and plot the Fst across the chromosome position

---
Note: Since R can not handle the full dataset, so we have to split each chromesome into two parts and run each part seperately.Then when you concatenate the Fst file when you calculta the Fst.

### Prepare the input file for Hierstats:

- Run vcftools to convert the vcf file into 12 coded files:

```bash
vcftools --vcf yourvcf.vcf --012 --out yourvcf
```

Then we can get several files:

yourvcf.012

yourvcf.012.indv

yourvcf.012.pos

-    Run [VCF012_to_hierfstat.py](https://github.com/lilei1/UK_barley_popgen/blob/main/scripts/VCF012_to_hierfstat.py) to convert those files into the input file for HierFstats:

The population are assigned according to if they are wild, landraces, or elite lines. The details of the population assignment can be avaiable [here](https://github.com/lilei1/UK_barley_popgen/blob/main/data/pop.txt)

```bash
./VCF012_to_hierfstat.py -012 yourvcf.012 -indv yourvcf.012.indv -pos yourvcf.012.pos -pop pop.txt
```

-    Run HierFstats to get the Fst values [HierFstats_Fst.R](https://github.com/lilei1/UK_barley_popgen/blob/main/scripts/HierFstats_Fst.R):

-   Then plot with [Fst_plot.R](https://github.com/lilei1/UK_barley_popgen/blob/main/scripts/Fst_plot.R)