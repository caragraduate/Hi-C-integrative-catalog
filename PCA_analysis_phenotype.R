library(PCAForQTL)
library(preprocessCore)
#source("http://www.bioconductor.org/biocLite.R")
#devtools::install_github(repo = 'zrmacc/RNOmni')
#biocLite("preprocessCore")
#install.packages("RNOmni")
library(RNOmni)
library(RColorBrewer)
library(janitor)
require(dplyr)

dataGeneExpressionFP<-read.table("phenotype/input/directory/impute_geom_scale_bin_middle_noNA12329noGM19204.bed", header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char = "", check.names = FALSE)
dataGeneExpressionFP_t <- t(dataGeneExpressionFP)

### Quantile Normalization
expr1<- as.matrix(dataGeneExpressionFP[,-(1:4)])
expr.qt <- preprocessCore::normalize.quantiles(as.matrix(expr1))

### inverse normal transform
expr.int.qt = t(apply(expr.qt, 1, RNOmni::RankNorm)) #when MARGIN=1, it applies over rows, whereas with MARGIN=2, it works over columns.
colnames(expr.int.qt) <- c('NA18534','NA18939','NA19036','NA19240','NA19650','NA19833','NA20509','NA20847','HG00096',
                           'HG00171','HG00514','HG00733','HG00864','HG01114','HG01505','HG01573','HG01596','HG02011',
                           'HG02018','HG02492','HG02587','HG03009','HG03065','HG03371','HG03683','HG03732','HG00512', 
                           'HG00513','HG00731','HG00732','NA19238','NA19239','NA18507','NA18505','NA18508','NA18486',
                           'NA19099','NA19141','NA18516','NA18522')
data_title = dataGeneExpressionFP[,(1:4)]
full_expr = cbind(data_title, expr.int.qt)
write.table(full_expr, file='phenotype/input/directory/impute_geom_scale_bin_qt_int_noNA12329noGM19204.bed', sep="\t",row.names = FALSE, quote = FALSE)

### visualize the quantile normalization 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) # 74 colors in RcolorBrewer
## visualized the first 10 TADs
#before quantile transformation
plot(density(na.omit(unlist(expr1[1, ]))), col=col_vec[1])
for(i in 2:10) lines(density(na.omit(unlist(expr1[i,]))), col=col_vec[i])
#for(i in 2:27) lines(density(na.omit(expr1[,i])), col=col_vec[i])
#after quantile transformation
plot(density(na.omit(expr.qt[1,])), col = col_vec[1])
for(i in 2:10) lines(density(na.omit(unlist(expr.qt[i, ]))), col=col_vec[i])

## visualized the all 27 samples
#before quantile transformation
plot(density(na.omit(unlist(expr1[,1]))), col=col_vec[1])
for(i in 2:41) lines(density(na.omit(unlist(expr1[ ,i]))), col=col_vec[i])
#for(i in 2:27) lines(density(na.omit(expr1[,i])), col=col_vec[i])
#after quantile transformation
plot(density(na.omit(expr.qt[ ,1])), col = col_vec[1])
for(i in 2:41) lines(density(na.omit(unlist(expr.qt[ ,i]))), col=col_vec[i])

## visualized the first 10 TADs
# before inverse normal transform
plot(density(na.omit(unlist(expr.qt[1, ]))), col = col_vec[1])
for(i in 2:10) lines(density(na.omit(unlist(expr.qt[i,]))), col=col_vec[i])
#plot(density(na.omit(unlist(expr1[1, ]))), col=col_vec[1])
#for(i in 2:10) lines(density(na.omit(unlist(expr1[i,]))), col=col_vec[i])
# after inverse normal transform
plot(density(na.omit(expr.int.qt[1,])), col = col_vec[1])
for(i in 2:10) lines(density(na.omit(unlist(expr.int.qt[i, ]))), col=col_vec[i])

## visualized the all 41 samples
#before inverse normal transform
#for(i in 2:10) lines(density(na.omit(unlist(expr1[ ,i]))), col=col_vec[i])
plot(density(na.omit(unlist(expr1[,1]))), col=col_vec[1])
for(i in 2:40) lines(density(na.omit(unlist(expr1[ ,i]))), col=col_vec[i])
#for(i in 2:27) lines(density(na.omit(expr1[,i])), col=col_vec[i])
#after inverse normal trnasform
plot(density(na.omit(expr.int.qt[ ,1])), col = col_vec[1])
for(i in 2:40) lines(density(na.omit(unlist(expr.int.qt[ ,i]))), col=col_vec[i])

expr.qt.int_new <- full_expr[,-(1:3)]
expr2 <- t(expr.qt.int_new)

expr_noconstant = remove_constant(expr2, na.rm= TRUE)
new_expr_noconstant <- expr_noconstant[ , colSums(is.na(expr_noconstant))==0] # colSums(is.na(df)) counts the number of NAs per column
apply(new_expr_noconstant, 2, var) == 0
final <- t(new_expr_noconstant)

## subset the original phenotype file by the common TAD_id column
data_title = dataGeneExpressionFP[,(1:4)]
df = merge(x = data_title, y = final, by = "bin_id")

df_1 <- df %>% relocate(bin_id, .after = IS_END_samp)

### generate the phenotype input file without any zero/identical value
expr_pca <- t(full_expr[,-(1:4)]) 
dim(expr_pca)
prcompResult<-prcomp(expr_pca,center=TRUE,scale.=TRUE) #This should take less than a minute.
PCs<-prcompResult$x 
dim(PCs)
write.table(PCs, file='phenotype/input/directory/bin_covariates_phenotypePCA.txt', sep="\t", quote = FALSE, col.names=NA)

importanceTable<-summary(prcompResult)$importance
PVEs<-importanceTable[2,]
sum(PVEs) #Theoretically, this should be 1.
plot(PVEs,xlab="PC index",ylab="PVE")

### choose the number of PCs
# automatic elbow detection method
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow)

# BE algorithm
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
resultRunBE<-PCAForQTL::runBE(expr_pca,B=20,alpha=0.05)
print(resultRunBE$numOfPCsChosen)

K_elbow<-resultRunElbow #2.
K_BE<-resultRunBE$numOfPCsChosen #3.
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE"),values=c(K_elbow,K_BE),
                         titleText="IS-TAD-qtl")

### using PCs in the QTL analysis
dataCovariates<-read.table("phenotype/input/directory/bin_covariates_final_genopca3_noNA12329noGM19204_top2.txt", sep="\t", header = TRUE, stringsAsFactors=FALSE)
dataCovariates_t = t(dataCovariates)
names(dataCovariates_t) <- as.character(dataCovariates_t[1,])
dataCovariates_t <- dataCovariates_t[-1,]
knownCovariates<-dataCovariates_t[,4:9]

identical(rownames(knownCovariates),rownames(expr_pca))
colnames(knownCovariates)[2:6]<-paste0("genotypePC",1:5) #This is to avoid potential confusion.
colnames(knownCovariates)[1]<-'sex'

PCsTop<-PCs[,1:K_BE]
knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(knownCovariates,PCsTop,unadjustedR2_cutoff=0.9)

#PCsTop<-scale(PCsTop) #Optional. Could be helpful for avoiding numerical inaccuracies.
covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)
