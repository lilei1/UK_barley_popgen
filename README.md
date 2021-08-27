# UK_barley_popgen
# This is a tutorial for UK collegues to run Hierfstats and plot the Fst across the chromosome position

---

## Dataset1: Ahmand's GBS data

### Prepare the input file for Hierstats:

- Run vcftools to convert the vcf file into 12 coded files:

```bash
vcftools --vcf yourvcf.vcf --012 --out yourvcf
```

Then we can get several files:

yourvcf.012

yourvcf.012.indv

yourvcf.012.pos

- Run [VCF012_to_hierfstat.py](https://github.com/lilei1/UK_barley_popgen/blob/main/scripts/VCF012_to_hierfstat.py) to convert those files into the input file for HierFstats:

The population are assigned according to if they are wild, landraces, or elite lines. The details of the population assignment can be avaiable [here]()

```bash
./VCF012_to_hierfstat.py -012 yourvcf.012 -indv yourvcf.012.indv -pos yourvcf.012.pos -pop Ahmand_302_pop.txt
```

---

-    Run HierFstats to get the Fst values:

```R

```

### Plot the Fst:

```R

```

---

### Assign population based on Structure results and run HierFstats:

### Run structure

- Prepare the input file for structure

    - Randomly choose 1000 SNPs from the vcf file

      Get rid of the header for the vcf file:

```bash
grep -v "#" WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf >no_head_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf
```

      Randomly pick 1000 SNPs from the vcf files:

```bash
 shuf -n 1000 no_head_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf >1000_no_head_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf
```

    Sort the vcf file:

```bash
sort -k1,1 -k2,2n 1000_no_head_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf >sorted_1000_no_head_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf
```

    Add the header to sorted vcf file

```bash
cat <(grep "#" WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf) sorted_1000_no_head_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf >1000_WBDC_July2016_production_PTP_filt_noIntro_SNPs.recode.vcf
```


    - Run PGD spider to convert the vcf file into the structure format:

    - Run Structure with 3 burning and 1000 replicates, K set as 1-10

    - assign the populations based on the structure results, and can be avaible [here](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/Files/Ahmad_306_popbyStructure.txt)

- Run [VCF012_to_hierfstat.py](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/VCF012_to_hierfstat.py) to convert those files into the input file for HierFstats

-    Run HierFstats to get the Fst values [HierFstats_Fst.R](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/HierFstats_Fst.R):

-   Then plot with [R_plot_Fst_str.R](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/R_plot_Fst_str.R) and [R_plot_Fst_gpspop.R](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/R_plot_Fst_gpspop.R)


## Dataset2: SNPs called from Exome resequencing data

### Prepare the input data:

The SNPs data we used are here:
`/panfs/roc/scratch/llei/GATK_filter_recalibrator/Inversion_SAMPLE/179_inversion/inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf`

Since R can not handle the full dataset, so we have to split each chromesome into two parts and run each part seperately.

-   We create bed files based on the length of each chromesome. In the dir `/panfs/roc/scratch/llei/GATK_filter_recalibrator/Inversion_SAMPLE/179_inversion/split_for_Fst`, here are 14 bed files:

chr1H_part1.bed

chr1H_part2.bed

chr2H_part1.bed

chr2H_part2.bed

chr3H_part1.bed

chr3H_part2.bed

chr4H_part1.bed

chr4H_part2.bed

chr5H_part1.bed

chr5H_part2.bed

chr6H_part1.bed

chr6H_part2.bed

chr7H_part1.bed

chr7H_part2.bed

-   Then use vcftools to cut the vcf files into 14 vcf files:

```bash
vcfintersect -b chr1H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr1H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf


vcfintersect -b chr2H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr2H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr3H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr3H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf


vcfintersect -b chr4H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr4H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr4H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr4H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr5H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr5H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr6H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr6H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr7H_part1.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf

vcfintersect -b chr7H_part2.bed .././inversion_179sample_onlySNP_Recal_50xfilter_heter_missing_filtering_mafgt0.002.recode.vcf >chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf
```

-   Use vcftools to creat the files for [Paul's script](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/VCF012_to_hierfstat.py):

```bash
#chr1H
vcftools --vcf chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing 

vcftools --vcf chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

#chr2H
vcftools --vcf chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

vcftools --vcf chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

#chr3H
vcftools --vcf chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

vcftools --vcf chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

#chr4H
vcftools --vcf chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

vcftools --vcf chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

#chr5H
vcftools --vcf chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

vcftools --vcf chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

#chr6H
vcftools --vcf chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

vcftools --vcf chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

#chr7H
vcftools --vcf chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

vcftools --vcf chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.recode.vcf --012 --out chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing

```

-   Then run [Paul's script](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/VCF012_to_hierfstat.py) to creat the hierfstats format to feed the hierfstats:

```bash
#chr1H
python /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt
#chr2H
python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

#chr3H
python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

#chr4H
python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

#chr5H
python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

#chr6H
python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

#chr7H
python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

python3 /Users/lilei/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Ahamd_GBS/VCF012_to_hierfstat.py -012 chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012 -indv chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.indv -pos chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.012.pos -pop GPS_pop.txt

```
-   Calculate the Fst with R script [HierFstats_Fst.R](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/HierFstats_Fst.R):

```bash
./HieterFstats.R chr1H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr1H_win2_strpop

./HieterFstats.R chr1H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr1H_win1_strpop


./HieterFstats.R chr2H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr2H_win2_strpop

./HieterFstats.R chr2H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr2H_win1_strpop


./HieterFstats.R chr3H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr3H_win2_strpop

./HieterFstats.R chr3H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr3H_win1_strpop


./HieterFstats.R chr4H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr4H_win2_strpop

./HieterFstats.R chr4H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr4H_win1_strpop


./HieterFstats.R chr5H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr5H_win2_strpop

./HieterFstats.R chr5H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr5H_win1_strpop


./HieterFstats.R chr6H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr6H_win2_strpop

./HieterFstats.R chr6H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr6H_win1_strpop


./HieterFstats.R chr7H_win2_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr7H_win2_strpop

./HieterFstats.R chr7H_win1_inversion_179sample_onlySNP_Recal_50xfilter_heter_missing.hierfstat chr7H_win1_strpop
```

-   Then we need to concatnate the two files into single file for each chromesome and sorted them based on the physical position.

-   Then plot 1H, 3H, 6H, the Fst with [Fst_plot.R](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/Fst_plot.R)

```shell
./Fst_plot.R chr1H_strpop_forR.Fst 1H 560 ~/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Exome chr1H_strpop_Fst_plot

./Fst_plot.R chr3H_strpop_forR.Fst 3H 700 ~/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Exome chr3H_strpop_Fst_plot

./Fst_plot.R chr4H_strpop_forR.Fst 4H 650 ~/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Exome chr4H_strpop_Fst_plot

./Fst_plot.R chr7H_strpop_forR.Fst 7H 650 ~/Backup_20170601/Li_befre_mac_clean/Barley_project/inversion_project/HierFsstats/Exome chr7H_strpop_Fst_plot
```
- since 2H and 5H is specific and need to highlight the putative inverted region and centeralmere region, I wrote additional [hard coded script](https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/Hierfstats/script/plot_Fst_2H_5H.R) to do that: