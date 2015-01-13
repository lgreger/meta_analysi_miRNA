# Meta analysis on AHUS and UCAM miRNA microarrays datasets 
## Liliana Greger

## Data description
The  dataset  consists of :
* The AHUS datasets are Agilent_miRNA_v2_4470B consisting of  4 benign, 15 invasive, 18 normal, no DICS samples
* Agilent_miRNA_v3_4470C consists of  23 benign, 55 invasive, 70 normal and  8 DICS
* UCAM dataset consists of  8 benign, 1283 invasive, 116 normal and 10 DICS

AHUS outputs comes from Agilents Feature Extraction Software (AFE). UCAM data  has been received  as a preprocessed matrix as described in the Nature paper (Nature. 2013 May 16;497(7449):378-82). 

```r
library(AgiMicroRna)
```

```
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
## 
## Loading required package: affy
## Loading required package: preprocessCore
## Loading required package: affycoretools
## Loading required package: GO.db
## Loading required package: AnnotationDbi
## Loading required package: DBI
## 
## Loading required package: KEGG.db
## 
## KEGG.db contains mappings based on older data because the original
##   resource was removed from the the public domain before the most
##   recent update was produced. This package should now be
##   considered deprecated and future versions of Bioconductor may
##   not have it available.  Users who want more current data are
##   encouraged to look at the KEGGREST or reactome.db packages
```

```
## Warning in namespaceImportMethods(ns, loadNamespace(j <- imp[[1L]],
## c(lib.loc, : No generic function found corresponding to requested imported
## methods for "summary" from package "AnnotationDbi" (malformed exports?)
```

```
## Warning in namespaceImportMethods(ns, loadNamespace(j <- imp[[1L]],
## c(lib.loc, : No generic function found corresponding to requested imported
## methods for "summary" from package "Category" (malformed exports?)
```

```
## Warning: found methods to import for function 'summary' but not the
## generic itself
```

```
## Warning in namespaceImportMethods(ns, loadNamespace(j <- imp[[1L]],
## c(lib.loc, : No generic function found corresponding to requested imported
## methods for "summary" from package "GOstats" (malformed exports?)
```

```
## Warning: found methods to import for function 'summary' but not the
## generic itself
```

```
## 
## Attaching package: 'AgiMicroRna'
## 
## The following object is masked from 'package:limma':
## 
##     readTargets
```

```r
library(Biobase)
library(MAMA)
```

```
## Loading required package: genefilter
## Loading required package: metaMA
## Loading required package: SMVar
## 
## Attaching package: 'metaMA'
## 
## The following object is masked from 'package:genefilter':
## 
##     rowVars
## 
## Loading required package: xtable
## Loading required package: multtest
## Loading required package: gtools
## Loading required package: grid
## Loading required package: GeneMeta
## 
## Attaching package: 'MAMA'
## 
## The following objects are masked from 'package:GeneMeta':
## 
##     multExpFDR, zScoreFDR, zScorePermuted, zScores
```

```r
library(plyr)
load("AHUS_agiMicroRNAdata.rdata")
```

## Preprocessing -  summarization
First we load  and preprocess the AHUS data. The use only one of the AHUS dataset which has DCIS samples. We choose to use  RMA summarisation of  the signal without background correction as described in BMC Research Notes 2010, 3:18.We also read the processed and already normalized UCAM data and the annotation file. 
We further annotate the samples using the annotation file and check if the data is paired. Only 5 samples from the AHUS Agilent
v2 and 31 samples Agilent v3 is paired.  No information on matched samples from UCAM.  Since we are interested in DCIS comparisons ( where there is no paired data available) we will not use paired data design.

```r
## read the annotation file  and the MIMAT mappings conversion file
aliases <- read.delim("aliases.txt", header = FALSE, sep = "") 
annotation <- read.delim("sample_annotationv3", header = FALSE, sep = "")
sampleannotation <- read.delim("sampleannotation_validatingvolinia.txt", header = TRUE, sep = "")
annot <- merge(annotation, sampleannotation, by.x = "V7", by.y = "fulldatapath")
## summarize ahus data with RMA
ahus_v3_4470C  = rmaMicroRna(AHUS_agiMicroRNAdata$Agilent_miRNA_v3_4470C, normalize = TRUE, background = FALSE)
```

```
## Calculating Expression
```

```r
## read the UCAM data
ucam = read.table("Agilent_ncRNA_60k_normalised_miRNA_expression_ORIGINAL.txt", header = TRUE, sep = "")
rowsd = apply(ucam[,-c(1,2)], MARGIN=1, FUN=sd)
ucam = ucam[order(rowsd, decreasing=TRUE), ]
ucam = ucam[!duplicated(ucam[,2]), ]
ucam_rownames  <- ucam$miRNA
rownames(ucam) <- ucam_rownames
ucam = as.matrix(ucam[, -c(1,2)])
colnames(ucam) = paste("UCAM_", colnames(ucam), sep="") 
ucam = ucam[, colnames(ucam) %in% sampleannotation$sample_id]
```
## Create the target files
We further annotate the samples using the annotation file.  No information on matched samples from UCAM.  

```r
## create ahus 3 target file 
ahus_v3_targets <- ahus_v3_4470C$targets
ahus_3 <- as.matrix(ahus_v3_targets$FileName)
status <- c()
 for ( i in 1: length(ahus_v3_targets$FileName)) {
status = c(status, as.vector(annot[which(annot$V7 %in% ahus_3[i]), 5]))}
var2 <- c()
for ( i in 1: length(ahus_v3_targets$FileName)) {
var2 = c(var2, as.vector(annot[which(annot$V7 %in% ahus_3[i]), 3]))}
ihc <- c()
for ( i in 1: length(ahus_v3_targets$FileName)) {
ihc = c(ihc, as.vector(annot[which(annot$V7 %in% ahus_3[i]), "IHC"]))}
GErep <- gsub("benign", "1", status)
GErep <- gsub("DCIS", "2", GErep)
GErep <- gsub("invasive", "3", GErep)
GErep <- gsub("normal", "4", GErep)
ahus_v3_targets <- lapply(ahus_v3_targets, function(x) gsub("microRNA/AHUS/data/Agilent_miRNA_v3_4470C/", "", x))
ahus_v3_targets  <- sapply(ahus_v3_targets , function(x) sub("_.*", "", x))
ahus_v3_targets <- as.data.frame(ahus_v3_targets)
targets_ahus_v3 <- cbind(ahus_v3_targets ,status,  var2,ihc,  GErep) 
targets_ahus_v3$ihc <- as.character.factor(targets_ahus_v3$ihc)
targets_ahus_v3$ihc[which(targets_ahus_v3$var2 == "DCIS") ] = "DCIS"
targets_ahus_v3$ihc <- as.factor(targets_ahus_v3$ihc)
##  create  ucam target file
status <- c()
for ( i in 1 : length(colnames(ucam))) {
status <- c(status, as.vector(annot[ annot$sample_id %in% grep( colnames(ucam)[i],annot$V9, value = "TRUE"), "V4"]))}
GErep <- gsub("benign", "1", status)
GErep <- gsub("DCIS", "2", GErep)
GErep <- gsub("invasive", "3", GErep)
GErep <- gsub("normal", "4", GErep)
var2<- c()
for ( i in 1 : length(colnames(ucam))) {
var2<- c(var2, as.vector(annot[ annot$sample_id %in% grep( colnames(ucam)[i],annot$V9, value = "TRUE"), "V2"]))}
ihc <- c()
for ( i in 1 : length(colnames(ucam))) {
ihc <- c(ihc, as.vector(annot[ annot$sample_id %in% grep( colnames(ucam)[i],annot$V9, value = "TRUE"), "IHC"]))}
targets_ucam = as.data.frame(cbind(FileName = colnames(ucam), status, var2,ihc,  GErep))
targets_ucam$ihc <- as.character.factor(targets_ucam$ihc)
targets_ucam$ihc[which(targets_ucam$var2 == "DCIS") ] = "DCIS"
targets_ucam$ihc <- as.factor(targets_ucam$ihc)
```
## Preprocessing - filtering
We then filter out the ahus data from control spots,  miRNA whose gMeanSignal is close to the expression of the negative controls(25\% limit) and  miRNA which are not expressed. We set up a limit of 75\% of the miRNA which must remain in at least one experimental condition. 288 features for ahus_v3 dataset are left  out of 961. 
UCAM data has been already filtered - 823 features.

```r
ahus_v3_filtered =filterMicroRna(ahus_v3_4470C , AHUS_agiMicroRNAdata$Agilent_miRNA_v3_4470C, control = TRUE, IsGeneDetected=TRUE,limIsGeneDetected=75,limNEG=25,targets = as.data.frame(targets_ahus_v3), verbose=TRUE,writeout=FALSE)
```

```
## FILTERING PROBES BY FLAGS 
## 
## 
## FILTERING BY ControlType 
## 
##    FEATURES BEFORE FILTERING:  961 
##  	FEATURES AFTER ControlType FILTERING:  939 
## ------------------------------------------------------ 
## FILTERING BY IsGeneDetected FLAG 
## 
## 	FLAG FILTERING OPTIONS - FLAG OK = 1 - limIsGeneDetected:  75 % 
## 	FEATURES AFTER IsGeneDetected FILTERING:  288 
## 	NON Gene Detected : 651 
## ------------------------------------------------------
```
## Creating the expression sets 
We then create the expression sets and later make sure the all datasets have common feature names in the exact order

```r
ahus_v3_eset =esetMicroRna(ahus_v3_filtered ,targets_ahus_v3 ,makePLOT=FALSE,verbose=TRUE)
```

```
## outPUT DATA: esetPROC 
## Features  Samples 
##      288      156 
## ------------------------------------------------------
```

```r
ucam_data <- targets_ucam
rownames(ucam_data) <- targets_ucam$FileName
pd <- new("AnnotatedDataFrame", data = ucam_data)
ucam_eset <- new("ExpressionSet", exprs = as.matrix(ucam) , phenoData = pd)
```
##  Conversion of   miRNA names to miRbase accession numbers

```r
aliases <- read.delim("aliases.txt", header = FALSE, sep = "") ## file from miRbase
aliases2 <- as.data.frame(aliases)
colnames(aliases2) <- c("V1", "V2")
library(data.table)
df <- data.table(aliases2, key="V1")
df <- df[, list(V2 = unlist(strsplit(as.character(V2), ";"))), by=V1]
df <- as.matrix(df)
df <- as.data.frame(df)
df <- ddply(df, "V2", summarize, ID = paste(V1, collapse="_")) 
featureNames(ahus_v3_eset) <- gsub("_v.*", "", featureNames(ahus_v3_eset))## remove the version attached to tehnames
miRnabase_ID2 <- data.frame(miRNA = c(), mirbase_an = c())
for (i in 1: length( featureNames(ahus_v3_eset ))) { 
miRBase_feature <- df[ which(df[,1] == featureNames(ahus_v3_eset)[i]),2]
tem_ID <- data.frame(miRNA   = featureNames(ahus_v3_eset)[i],  mirbase_an =  as.character(miRBase_feature)   )
miRnabase_ID2  <- rbind( miRnabase_ID2, tem_ID)}
miRnabase_id<- miRnabase_ID2[!duplicated(miRnabase_ID2), ]
featureNames(ahus_v3_eset) <- miRnabase_id[,2]
ind <- grep("put.*", featureNames(ucam_eset))
ucam_eset<- ucam_eset[-ind, ]
ind <- grep("DQ.*", featureNames(ucam_eset))
ucam_eset<- ucam_eset[-ind, ]
ind <- grep("FANTOM.*", featureNames(ucam_eset))
ucam_eset<- ucam_eset[-ind, ]
miRnabase_ID3 <- data.frame(miRNA = c(), mirbase_an = c())
for (i in 1: length( featureNames(ucam_eset ))) { 
miRBase_feature <- df[ which(df[,1] == featureNames(ucam_eset)[i]),2 ]
tem_ID <- data.frame(miRNA   = featureNames(ucam_eset)[i],  mirbase_an =  as.character(miRBase_feature)   )
miRnabase_ID3  <- rbind( miRnabase_ID3, tem_ID)
}
miRnabase_id<- miRnabase_ID3[!duplicated(miRnabase_ID3), ]
rownames(miRnabase_id) <- 1:nrow(miRnabase_id)
d <- as.character(miRnabase_id$mirbase_an)
d[389] <- "MIMAT0003257*"
featureNames(ucam_eset) <- d
```
## Find common features for both datasets

```r
ahus_v3_features <- featureNames(ahus_v3_eset) ## 266
ucam_features <- featureNames(ucam_eset)
all_features <- intersect(ahus_v3_features, ucam_features) ## ahus_v2 does not have DCIS data
ahus_v3_common_features <- ahus_v3_eset[featureNames(ahus_v3_eset)  %in% all_features , ]
ahus_v3_common_features <- ahus_v3_common_features[order( featureNames(ahus_v3_common_features)),]
ucam_common_features <- ucam_eset[featureNames(ucam_eset)  %in% all_features , ]
ucam_common_features <- ucam_common_features[order( featureNames(ucam_common_features)),]
```
## Perform meta analysis
We use the combined p value meta analysis using limma implemented in the MAMA package 
## First, prepare the ahus sub-dataseta for meta analysis

```r
ahus_v3_dcis_lumA <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "LumA")]
ahus_v3_dcis_lumA <- ahus_v3_dcis_lumA[ , ahus_v3_dcis_lumA$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_lumA))$var2 <- factor( pData(phenoData(ahus_v3_dcis_lumA))$var2) 
pData(phenoData(ahus_v3_dcis_lumA))$FileName <- factor( pData(phenoData(ahus_v3_dcis_lumA))$FileName)
pData(phenoData(ahus_v3_dcis_lumA))$status <- factor( pData(phenoData(ahus_v3_dcis_lumA))$status) 
ahus_v3_dcis_lumB <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "LumB")]
ahus_v3_dcis_lumB <- ahus_v3_dcis_lumB[ , ahus_v3_dcis_lumB$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_lumB))$var2 <- factor( pData(phenoData(ahus_v3_dcis_lumB))$var2) 
pData(phenoData(ahus_v3_dcis_lumB))$FileName <- factor( pData(phenoData(ahus_v3_dcis_lumB))$FileName) 
pData(phenoData(ahus_v3_dcis_lumB))$status <- factor( pData(phenoData(ahus_v3_dcis_lumB))$status) 
ahus_v3_dcis_Basal <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in% c("DCIS", "Basal")]
ahus_v3_dcis_Basal<- ahus_v3_dcis_Basal[ , ahus_v3_dcis_Basal$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_Basal))$var2 <- factor(pData(phenoData(ahus_v3_dcis_Basal))$var2) 
pData(phenoData(ahus_v3_dcis_Basal))$FileName <- factor(pData(phenoData(ahus_v3_dcis_Basal))$FileName)
pData(phenoData(ahus_v3_dcis_Basal))$status <- factor(pData(phenoData(ahus_v3_dcis_Basal))$status)
ahus_v3_dcis_Normal <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "Normal")]
ahus_v3_dcis_Normal <- ahus_v3_dcis_Normal[ , ahus_v3_dcis_Normal$status %in% c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_Normal))$var2 <- factor(pData(phenoData(ahus_v3_dcis_Normal))$var2) 
pData(phenoData(ahus_v3_dcis_Normal))$FileName <- factor(pData(phenoData(ahus_v3_dcis_Normal))$FileName) 
pData(phenoData(ahus_v3_dcis_Normal))$status <- factor(pData(phenoData(ahus_v3_dcis_Normal))$status) 
ahus_v3_dcis_Her2 <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "Her2")]
ahus_v3_normal_Her2<- ahus_v3_dcis_Her2[ , ahus_v3_dcis_Her2$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_Her2))$var2 <-factor( pData(phenoData(ahus_v3_dcis_Her2))$var2)  
pData(phenoData(ahus_v3_dcis_Her2))$FileName <- factor( pData(phenoData(ahus_v3_dcis_Her2))$FileName) 
pData(phenoData(ahus_v3_dcis_Her2))$status <- factor( pData(phenoData(ahus_v3_dcis_Her2))$status) 
ahus_v3_normal_DCIS <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("normal", "DCIS")]
pData(phenoData(ahus_v3_normal_DCIS))$status <- factor( pData(phenoData(ahus_v3_normal_DCIS))$status)# drop previous factor levels
pData(phenoData(ahus_v3_normal_DCIS))$FileName <- factor( pData(phenoData(ahus_v3_normal_DCIS))$FileName)# 
ahus_v3_normaltobenign <-  ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("normal", "benign")]
pData(phenoData(ahus_v3_normaltobenign ))$status <- factor( pData(phenoData(ahus_v3_normaltobenign ))$status)# drop previous factor levels
pData(phenoData(ahus_v3_normaltobenign ))$FileName <- factor( pData(phenoData(ahus_v3_normaltobenign ))$FileName)# 
ahus_v3_invasive_DCIS <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_invasive_DCIS))$status <- factor( pData(phenoData(ahus_v3_invasive_DCIS))$status)# drop previous factor levels
pData(phenoData(ahus_v3_invasive_DCIS))$FileName <- factor( pData(phenoData(ahus_v3_invasive_DCIS))$FileName)# 
ahus3_dcis_HER2neg_ERneg_PGRneg  <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2neg_ERneg_PGRneg <- ahus3_dcis_HER2neg_ERneg_PGRneg[, ahus3_dcis_HER2neg_ERneg_PGRneg$ihc %in% c("DCIS", "HER2neg_ERneg_PGRneg")]
pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$ihc) 
pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$FileName) 
pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$status <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$status) 
ahus3_dcis_HER2neg_ERpos  <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2neg_ERpos <- ahus3_dcis_HER2neg_ERpos[, ahus3_dcis_HER2neg_ERpos$ihc %in% c("DCIS", "HER2neg_ERpos")]
pData(phenoData(ahus3_dcis_HER2neg_ERpos))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERpos))$ihc) 
pData(phenoData(ahus3_dcis_HER2neg_ERpos))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERpos))$FileName) 
pData(phenoData(ahus3_dcis_HER2neg_ERpos))$status <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERpos))$status) 
ahus3_dcis_HER2pos_ERneg   <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2pos_ERneg  <- ahus3_dcis_HER2pos_ERneg [, ahus3_dcis_HER2pos_ERneg $ihc %in% c("DCIS", "HER2pos_ERneg")]
pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$ihc) 
pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$FileName) 
pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$status <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$status) 
ahus3_dcis_HER2pos_ERpos   <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2pos_ERpos  <- ahus3_dcis_HER2pos_ERpos [, ahus3_dcis_HER2pos_ERpos $ihc %in% c("DCIS", "HER2pos_ERpos")]
pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$ihc) 
pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$FileName) 
pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$status <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$status)
```
## Prepare ucam sub-datasets for meta analysis


```r
ucam_normallDCIS<- ucam_common_features[, ucam_common_features$status %in% c("normal", "DCIS") ]   
pData(phenoData(ucam_normalDCIS))$status  <- factor(pData(phenoData(ucam_normalDCIS))$status )
```

```
## Error in pData(phenoData(ucam_normalDCIS)): error in evaluating the argument 'object' in selecting a method for function 'pData': Error in phenoData(ucam_normalDCIS) : 
##   error in evaluating the argument 'object' in selecting a method for function 'phenoData': Error: object 'ucam_normalDCIS' not found
```

```r
pData(phenoData(ucam_normalDCIS))$FileName  <- factor(pData(phenoData(ucam_normalDCIS))$FileName )
```

```
## Error in pData(phenoData(ucam_normalDCIS)): error in evaluating the argument 'object' in selecting a method for function 'pData': Error in phenoData(ucam_normalDCIS) : 
##   error in evaluating the argument 'object' in selecting a method for function 'phenoData': Error: object 'ucam_normalDCIS' not found
```

```r
ucam_invasiveDCIS<- ucam_common_features[, ucam_common_features$status %in% c("invasive", "DCIS") ]
pData(phenoData(ucam_invasiveDCIS))$status  <- factor(pData(phenoData(ucam_invasiveDCIS))$status )               
pData(phenoData(ucam_invasiveDCIS))$FileName  <- factor(pData(phenoData(ucam_invasiveDCIS))$FileName )  
ucam_normalTobenign<- ucam_common_features[, ucam_common_features$status %in% c("benign", "normal") ]
pData(phenoData(ucam_normalTobenign))$status  <- factor(pData(phenoData(ucam_normalTobenign))$status )               
pData(phenoData(ucam_normalTobenign))$FileName  <- factor(pData(phenoData(ucam_normalTobenign))$FileName )
ucam_dcis_lumA <- ucam_common_features[, ucam_common_features$var2 %in% c("DCIS", "LumA") ]
ucam_dcis_lumA <- ucam_dcis_lumA[ , ucam_dcis_lumA$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_lumA))$var2  <- factor(pData(phenoData(ucam_dcis_lumA))$var2 )               
pData(phenoData(ucam_dcis_lumA))$FileName  <- factor(pData(phenoData(ucam_dcis_lumA))$FileName )
pData(phenoData(ucam_dcis_lumA))$status <- factor(pData(phenoData(ucam_dcis_lumA))$status)
ucam_dcis_lumB <- ucam_common_features[, ucam_common_features$var2 %in% c("DCIS", "LumB") ]
ucam_dcis_lumB <- ucam_dcis_lumB[ , ucam_dcis_lumB$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_lumB))$var2  <- factor(pData(phenoData(ucam_dcis_lumB))$var2 )               
pData(phenoData(ucam_dcis_lumB))$FileName  <- factor(pData(phenoData(ucam_dcis_lumB))$FileName )
pData(phenoData(ucam_dcis_lumB))$status <- factor(pData(phenoData(ucam_dcis_lumB))$status)
ucam_dcis_Her2 <- ucam_common_features[, ucam_common_features$var2 %in% c("DCIS", "Her2") ]
ucam_dcis_Her2 <- ucam_dcis_Her2[ , ucam_dcis_Her2$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_Her2))$var2  <- factor(pData(phenoData(ucam_dcis_Her2))$var2 )               
pData(phenoData(ucam_dcis_Her2))$FileName  <- factor(pData(phenoData(ucam_dcis_Her2))$FileName )
pData(phenoData(ucam_dcis_Her2))$status  <- factor(pData(phenoData(ucam_dcis_Her2))$status )
ucam_dcis_Basal <- ucam_common_features[ , ucam_common_features$var2 %in% c("DCIS", "Basal")]
ucam_dcis_Basal<- ucam_dcis_Basal[ , ucam_dcis_Basal$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_Basal))$var2 <- factor(pData(phenoData(ucam_dcis_Basal))$var2) 
pData(phenoData(ucam_dcis_Basal))$FileName <- factor(pData(phenoData(ucam_dcis_Basal))$FileName)
pData(phenoData(ucam_dcis_Basal))$status <- factor(pData(phenoData(ucam_dcis_Basal))$status)
ucam_dcis_Normal <- ucam_common_features[ , ucam_common_features$var2 %in%  c("DCIS", "Normal")]
ucam_dcis_Normal <- ucam_dcis_Normal[ , ucam_dcis_Normal$status %in% c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_Normal))$var2 <- factor(pData(phenoData(ucam_dcis_Normal))$var2) 
pData(phenoData(ucam_dcis_Normal))$FileName <- factor(pData(phenoData(ucam_dcis_Normal))$FileName) 
pData(phenoData(ucam_dcis_Normal))$status <- factor(pData(phenoData(ucam_dcis_Normal))$status) 
ucam_dcis_HER2neg_ERneg_PGRneg  <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2neg_ERneg_PGRneg <- ucam_dcis_HER2neg_ERneg_PGRneg[, ucam_dcis_HER2neg_ERneg_PGRneg$ihc %in% c("DCIS", "HER2neg_ERneg_PGRneg")]
pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$ihc <- factor( pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$ihc) 
pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$FileName <- factor( pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$FileName) 
pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$status <- factor( pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$status) 
ucam_dcis_HER2neg_ERpos  <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2neg_ERpos <- ucam_dcis_HER2neg_ERpos[, ucam_dcis_HER2neg_ERpos$ihc %in% c("DCIS", "HER2neg_ERpos")]
pData(phenoData(ucam_dcis_HER2neg_ERpos))$ihc <- factor( pData(phenoData(ucam_dcis_HER2neg_ERpos))$ihc) 
pData(phenoData(ucam_dcis_HER2neg_ERpos))$FileName <- factor( pData(phenoData(ucam_dcis_HER2neg_ERpos))$FileName) 
pData(phenoData(ucam_dcis_HER2neg_ERpos))$status <- factor( pData(phenoData(ucam_dcis_HER2neg_ERpos))$status) 
ucam_dcis_HER2pos_ERneg   <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2pos_ERneg  <- ucam_dcis_HER2pos_ERneg [, ucam_dcis_HER2pos_ERneg $ihc %in% c("DCIS", "HER2pos_ERneg")]
pData(phenoData(ucam_dcis_HER2pos_ERneg ))$ihc <- factor( pData(phenoData(ucam_dcis_HER2pos_ERneg ))$ihc) 
pData(phenoData(ucam_dcis_HER2pos_ERneg ))$FileName <- factor( pData(phenoData(ucam_dcis_HER2pos_ERneg ))$FileName) 
pData(phenoData(ucam_dcis_HER2pos_ERneg ))$status <- factor( pData(phenoData(ucam_dcis_HER2pos_ERneg ))$status) 
ucam_dcis_HER2pos_ERpos   <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2pos_ERpos  <- ucam_dcis_HER2pos_ERpos [, ucam_dcis_HER2pos_ERpos $ihc %in% c("DCIS", "HER2pos_ERpos")]
pData(phenoData(ucam_dcis_HER2pos_ERpos ))$ihc <- factor( pData(phenoData(ucam_dcis_HER2pos_ERpos ))$ihc) 
pData(phenoData(ucam_dcis_HER2pos_ERpos ))$FileName <- factor( pData(phenoData(ucam_dcis_HER2pos_ERpos ))$FileName) 
pData(phenoData(ucam_dcis_HER2pos_ERpos ))$status <- factor( pData(phenoData(ucam_dcis_HER2pos_ERpos ))$status)
```

## Prepare the comined sub-datatests


```r
all_normal_DCIS <- new("MetaArray", GEDM = list( exprs(ahus_v3_normal_DCIS), exprs(ucam_normalDCIS)), clinical = list(pData(phenoData(ahus_v3_normal_DCIS)), pData(phenoData(ucam_normalDCIS) )) , datanames = c( "ahus_v3_4470C_normal_DCIS", "ucam_normal_DCIS"))
```

```
## Error in exprs(ucam_normalDCIS): error in evaluating the argument 'object' in selecting a method for function 'exprs': Error: object 'ucam_normalDCIS' not found
```

```r
all_invasive_DCIS <- new("MetaArray", GEDM = list( exprs(ahus_v3_invasive_DCIS), exprs(ucam_invasiveDCIS )), clinical = list(pData(phenoData(ahus_v3_invasive_DCIS)), pData(phenoData(ucam_invasiveDCIS ) )) , datanames = c( "ahus_v3_4470C_invasive_DCIS", "ucam_invasive_DCIS"))
all_normnal_to_benign <- new("MetaArray", GEDM = list( exprs(ahus_v3_normaltobenign), exprs(ucam_normalTobenign )), clinical = list(pData(phenoData(ahus_v3_normaltobenign)), pData(phenoData(ucam_normalTobenign ) )) , datanames = c( "ahus_v3_4470C_invasive_DCIS", "ucam_normalTobenign"))
all_dcis_LumA  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_lumA), exprs(ucam_dcis_lumA )), clinical = list(pData(phenoData(ahus_v3_dcis_lumA)), pData(phenoData(ucam_dcis_lumA ) )) , datanames = c( "ahus_v3_dcis_lumA", "ucam_dcis_lumA"))
all_dcis_LumB  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_lumB), exprs(ucam_normal_lumB )), clinical = list(pData(phenoData(ahus_v3_dcis_lumB)), pData(phenoData(ucam_dcis_lumB ) )) , datanames = c( "ahus_v3_dcis_lumB", "ucam_dcis_lumB"))
```

```
## Error in exprs(ucam_normal_lumB): error in evaluating the argument 'object' in selecting a method for function 'exprs': Error: object 'ucam_normal_lumB' not found
```

```r
all_dcis_Her2  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_Her2), exprs(ucam_dcis_Her2 )), clinical = list(pData(phenoData(ahus_v3_dcis_Her2)), pData(phenoData(ucam_dcis_Her2 ) )) , datanames = c( "ahus_v3_dcis_lumB", "ucam_dcis_Her2"))
all_dcis_Normal  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_Normal),exprs(ucam_dcis_Normal )), clinical =list(pData(phenoData(ahus_v3_dcis_Normal)), pData(phenoData(ucam_dcis_Normal ))) , datanames = c( "ahus_v3_dcis_Normal", "ucam_dcis_Normal"))
all_dcis_Basal  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_Basal),exprs(ucam_dcis_Basal )), clinical = list(pData(phenoData(ahus_v3_dcis_Basal)),pData(phenoData(ucam_dcis_Basal) )) , datanames = c( "ahus_v3_dcis_Basal", "ucam_dcis_Basal"))
all_dcis_HER2neg_ERneg_PGRneg <- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2neg_ERneg_PGRneg),exprs(ucam_dcis_HER2neg_ERneg_PGRneg )), clinical =list(pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg)), pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg ))) , datanames = c( "ahus3_dcis_HER2neg_ERneg_PGRneg", "ucam_dcis_HER2neg_ERneg_PGRneg"))
all_dcis_HER2pos_ERpos<- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2pos_ERpos),exprs(ucam_dcis_HER2pos_ERpos )), clinical = list(pData(phenoData(ahus3_dcis_HER2pos_ERpos)), pData(phenoData(ucam_dcis_HER2pos_ERpos))) , datanames = c( "ahus3_dcis_HER2pos_ERpos", "ucam_dcis_HER2pos_ERpos"))
all_dcis_HER2pos_ERneg<- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2pos_ERneg),exprs(ucam_dcis_HER2pos_ERneg )), clinical = list(pData(phenoData(ahus3_dcis_HER2pos_ERneg)), pData(phenoData(ucam_dcis_HER2pos_ERneg))) , datanames = c( "ahus3_dcis_HER2pos_ERneg", "ucam_dcis_HER2pos_ERneg"))
all_dcis_HER2neg_ERpos<- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2neg_ERpos),exprs(ucam_dcis_HER2neg_ERpos )), clinical = list(pData(phenoData(ahus3_dcis_HER2neg_ERpos)), pData(phenoData(ucam_dcis_HER2neg_ERpos))) , datanames = c( "ahus3_dcis_HER2neg_ERpos", "ucam_dcis_HER2neg_ERpos"))
```
## Meta analysis using limma


```r
results_pval1 <- metaMA(all_normal_DCIS, varname = "status", moderated = "limma",   which = "pval") 
```

```
## Error in GEDM(data): object 'all_normal_DCIS' not found
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
```

```
## Error in eval(expr, envir, enclos): object 'results_pval1' not found
```

```r
write.table(test, "DCIS_to_normal_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_invasive_DCIS, varname = "status", moderated = "limma",   which = "pval") 
```

```
##   DE  IDD Loss  IDR  IRR 
##    0    0   28  NaN  100
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "DCIS_to_invasive_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_normnal_to_benign, varname = "status", moderated = "limma",   which = "pval")
```

```
##     DE    IDD   Loss    IDR    IRR 
## 169.00  10.00  27.00   5.92  14.52
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "benign_to_normal_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_dcis_LumA, varname = "var2", moderated = "limma",   which = "pval") 
```

```
##    DE   IDD  Loss   IDR   IRR 
## 10.00  2.00 16.00 20.00 66.67
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_lumA_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_dcis_LumB, varname = "var2", moderated = "limma",   which = "pval") 
```

```
## Error in GEDM(data): object 'all_dcis_LumB' not found
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_lumB_limma.txt", sep = "\t", quote = FALSE)
write.table(test, "dcis_to_lumB_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_dcis_Her2, varname = "var2", moderated = "limma",   which = "pval") 
```

```
##   DE  IDD Loss  IDR  IRR 
##   11    0   11    0   50
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_her2_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_dcis_Normal, varname = "var2", moderated = "limma",   which = "pval") 
```

```
##    DE   IDD  Loss   IDR   IRR 
##  7.00  3.00  1.00 42.86 20.00
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_normal(invasive)_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(all_dcis_Basal, varname = "var2", moderated = "limma",   which = "pval") 
```

```
##    DE   IDD  Loss   IDR   IRR 
##  7.00  0.00  2.00  0.00 22.22
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "basal_to_dcis(invasive)_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(  all_dcis_HER2neg_ERneg_PGRneg, varname = "ihc", moderated = "limma",   which = "pval") 
```

```
##    DE   IDD  Loss   IDR   IRR 
##  2.00  0.00  1.00  0.00 33.33
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_HER2neg_ERneg_PGRneg_limma.txt", sep = "\t", quote = FALSE)
ults_pval1 <- metaMA(  all_dcis_HER2neg_ERpos, varname = "ihc", moderated = "limma",   which = "pval") 
```

```
##    DE   IDD  Loss   IDR   IRR 
##  1.00  0.00 11.00  0.00 91.67
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_HER2neg_ERpos_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(  all_dcis_HER2pos_ERneg, varname = "ihc", moderated = "limma",   which = "pval") 
```

```
##    DE   IDD  Loss   IDR   IRR 
## 15.00  5.00  2.00 33.33 16.67
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_HER2pos_ERneg_limma.txt", sep = "\t", quote = FALSE)
results_pval1 <- metaMA(  all_dcis_HER2pos_ERpos, varname = "ihc", moderated = "limma",   which = "pval") 
```

```
##   DE  IDD Loss  IDR  IRR 
##    3    0    3    0   50
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  df[ which(df$ID ==  DE) , 1]
test_st <-  results_pval1$TestStatistic[i]
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), FC = test_st)
test <- rbind(test, test2)
}
write.table(test, "dcis_to_HER2pos_ERpos_limma.txt", sep = "\t", quote = FALSE)
```









