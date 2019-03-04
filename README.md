# asaMap

Article
https://www.biorxiv.org/content/10.1101/014001v3

The implementation is still a work in progress. More features might be added at some point:
-- Poisson regression
-- larger number of populations (crazy high number of parameters - so maybe not the best idea)

# Install
git clone https://github.com/angsd/asaMap.git;
cd asaMap;
make

# Options
```
Options:

   -p <filename>       plink prefix filename
   -o <filename>       output filename
   -c <filename>       covariates filename
   -y <filename>       phenotypes
   -a <filename>       admixproportions (for source pop1)
   -Q <filename>       .Q file from ADMIXTURE
   -f <filename>       allele frequencies (.P file)
   -m <INT>            model 0=add 1=rec (default: 0)
   -l <INT>            regression 0=linear regression, 1=logistic regression (default: 0)
   -b <filename>       file containing the start
   -i <INT>            maximum iterations (default: 40)
   -t <FLOAT>          tolerance for breaking EM (default: 0.0001)
   -r <INT>            seed for rand
   -P <INT>            number of threads
   -e <INT>            estimate standard error of coefficients (0: no (default), 1: yes)
   -w <INT>            run M0/R0 model that models effect of other allele (0: no, 1: yes (default))


All files must be specified: -p -c -y -a/-Q -f -o
```
All files must be specified: -p -c -y -a/-Q -f -o

# Run example using ADMIXTURE

```
PH=pheno.files
COV=cov.file
PF=plinkFile #without postfix e.g. no .bim / .fam / .bed
#run admixture
admixture $PF.bed 2

#run asaMap with .Qfile
./asaMap -p $PF -o out -c $COV -y $PH -Q $PF.2.Q -f $PF.2.P

#get admixture proportions for population 1
cut -f1 -d" " $PF.2.Q > Q

#run asaMap with admix proportions
./asaMap -p $PF -o out -c $COV -y $PH -a Q -f $PF.2.P
```

### p-values
Easy to optain in R

Using R/getPvalues.R (will create file: out.res.Pvalues)
```
Rscript R/getPvalues.R out.res
```

```
res <- read.table("out.res",head=T)
pval <- function(a,b,df)
1-pchisq(-2*(a-b),df=df)
#pvalues for M1 vs. M5
pval(res$llh.M1.,res$llh.M5.,df=2)
#pvalues for M1 vs. M4
pval(res$llh.M1.,res$llh.M4.,df=1)
```

### models

| model | parameters | notes | #effect Parameters |
| --- | --- | --- | --- |
| M0 | (beta_1, beta_2, delta_1) in R^3 | effect of non-assumed effect allele | 1 |
| M1 | (beta_1, beta_2) in R^2  | population specific effects | 2 |
| M2 | beta_1=0, beta_2 in R | no effect in population 1  | 1 |
| M3 | beta_1 in R, beta_2=0 | no effect in population 2 | 1 |
| M4 | beta_1=beta_2 in R | same effect in both populations | 1 |
| M5 | beta_1=beta_2=0 | no effect in any population | 0 |

| model | parameters | notes | #effect Parameters |
| --- | --- | --- | --- |
| R0 | (beta_1, beta_m, beta_2, delta_1, delta_2) in R^5 | recessive effect of non-assumed effect alleles | 2 |
| R1 | (beta_1, beta_m, beta_2) in R^3 | population specific effects | 3 |
| R2 | beta_1 in R, beta_m=beta_2 in R | same effect when one or both variant alleles are from pop 2 | 2 |
| R3 | beta_1=beta_m in R, beta_2 in R | same effect when one or both variant alleles are from pop 1 | 2 |
| R4 | beta_1 in R, beta_m=beta_2=0 | only an effect when both variant alleles are from pop 1 | 1 |
| R5 | beta_1=beta_m=0, beta_2 in R | only an effect when both variant alleles are from pop 2 | 1 |
| R6 | beta_1=beta_m=beta_2 in R | same effect regardless of ancestry | 1 |
| R7 | beta_1=beta_m=beta_2=0 | no effect in any population | 0 |


# Input files
### Genotypes
plink binary (.bim .bam .fam)

### Phenotypes (response)
A file with each individuals phenotypes on each line. e.g., has to have same number of rows as .fam file.
```
>head pheno
-0.712027291121767
-0.158413122435864
-1.77167888612947
-0.800940619551485
0.3016297021294
0.596892506547882
-0.661786423692485
-0.405728519330873
-1.04224674183241
0.0881848860116932
```
### extra covariates (in addition to the intercept and genotypes)
A file where each column is a covariate and each row is an individual - should NOT have columns of 1s for intercept (intercept will be included automatically).
This file has to have same number of rows as phenotype file and .fam file.
```
>head cov
0.0127096117618385 -0.0181281029917176 -0.0616739439849275 -0.0304606694443973
0.0109944672768584 -0.0205785925514037 -0.0547523583405743 -0.0208813157640705
0.0128395346453956 -0.0142116856067135 -0.0471689997039534 -0.0266186436009881
0.00816783754598649 -0.0189271733933446 -0.0302259313905976 -0.0222247658768436
0.00695928218989132 -0.0089960963981644 -0.0384886176827146 -0.0126490197701687
0.00908359304129912 -0.019562503526549 -0.0276058491506046 -0.0202388414332682
0.0193657006952317 -0.0219605099975189 -0.050537627417191 -0.0236411635865132
0.015862252334236 -0.0134969241244036 -0.0336244748700029 -0.0222294652006281
0.0194100156955457 -0.0371103372950621 0.00813012568415838 -0.015311879434991
0.0190516629849255 -0.0194012185542486 -0.0413589828106922 -0.0292318169458017
```
### admixture proportions


Right now the implementation is just for two populations. Either using a .Qfile from ADMIXTURE - has to have same number of rows as .fam file.

```
>head $PF.2.Q
0.398834 0.601166
0.491314 0.508686
0.399352 0.600648
0.634947 0.365053
0.000010 0.999990
0.467470 0.532530
0.497885 0.502115
0.800272 0.199728
0.999990 0.000010
0.032108 0.967892
```

Or using a column of the admixture proportions of the first population - has to have same number of rows as .fam file.
```
>head Q
0.398834
0.491314
0.399352
0.634947
0.000010
0.467470
0.497885
0.800272
0.999990
0.032108
```

### ancestral allele frequencies
A files where each column is a population and each row is a SNP. This file has to have same number of rows as .bim file.
```
>head P
0.907268 0.903723
0.498913 0.543486
0.865254 0.819332
0.737799 0.738876
0.812931 0.823395
0.783896 0.785741
0.711677 0.865518
0.592533 0.587921
0.834116 0.850966
0.837011 0.852711
```
