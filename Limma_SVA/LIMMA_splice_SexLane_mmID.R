#################################
###To Run copy this script onto the scc4 cluster
###From the folder containing this script
##Copy DEXseq_Rrun_Comb_NoCage_Lane.sh to the folder
########NOTE: EDIT "#$ -M reeder@bu.edu" to include your e-mail (line 8 of .sh file)

########RUN COMMAND###########
##qsub -P rufy1 LIMMA_splice_SexLane_mmID.sh

###################################

require(limma)
require(edgeR)

#Read in count files
inDir<-"/restricted/projectnb/rufy1/Eric/tophat_4DEXseq/DEXseq_txt_output"
countFiles = list.files(inDir, pattern=".txt$", full.names=TRUE)

flattenedFile = list.files("/restricted/projectnb/rufy1/reference/ENSEMBL/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/", pattern="gff$", full.names=TRUE)

#Read in table of samples
sampleTable = data.frame(
  row.names = c( "CB_10_L7", "CB_10_L8", "CB_11_L7", "CB_11_L8", "CB_12_L7", "CB_12_L8", "CB_13_L8", "CB_14_L7", "CB_14_L8","CB_15_L7", "CB_15_L8","CB_16_L7", "CB_16_L8","CB_1_L5", "CB_1_L6","CB_2_L5","CB_2_L6","CB_3_L5", "CB_3_L6","CB_4_L5", "CB_4_L6","CB_5_L5", "CB_5_L6","CB_6_L5", "CB_6_L6","CB_7_L5", "CB_7_L6","CB_8_L5", "CB_8_L6","CB_9_L7", "CB_9_L8"),
  condition = c( "congenic","congenic","B6","B6","congenic","congenic","B6","congenic","congenic","B6","B6","congenic","congenic","B6","B6","congenic", "congenic","B6","B6","congenic","congenic","B6","B6","congenic","congenic","B6","B6","congenic","congenic","B6","B6"),
  Lane = c( "L7", "L8", "L7", "L8", "L7", "L8", "L8", "L7", "L8","L7", "L8","L7", "L8","L5", "L6","L5","L6","L5", "L6","L5", "L6","L5", "L6","L5", "L6","L5", "L6","L5", "L6","L7", "L8" ),
  Sex = c("M", "M", "F", "F", "F","F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "M", "M",  "M", "M", "F", "F", "F", "F", "M", "M"),
  Cage = c("C85", "C85", "C90", "C90", "C84", "C84", "C83", "C82", "C82", "C83", "C83", "C82", "C82", "C89", "C89", "C87", "C87","C81", "C81", "C80", "C80", "C89", "C89", "C87", "C87", "C86", "C86", "C90", "C90", "C90", "C90"),
  ID = c("CB10", "CB10", "CB11", "CB11", "CB12", "CB12", "CB13", "CB14", "CB14", "CB15", "CB15", "CB16", "CB16", "CB1","CB1", "CB2","CB2", "CB3","CB3", "CB4","CB4", "CB5","CB5", "CB6", "CB6", "CB7", "CB7", "CB8", "CB8", "CB9", "CB9"))

#Order sampleTable by file names
sampFiles<-list.files(inDir, pattern=".txt$", full.names=FALSE)
sampFiles<-sub("DEXseq_", "", sampFiles)
sampFiles<-sub(".txt", "", sampFiles)
sampFiles<-sub("00", "", sampFiles)
sampFiles<-sub("-", "_", sampFiles)
sampleTable<-sampleTable[order(match(row.names(sampleTable), sampFiles)),]

#Combine count files
countMat<-list()
for(i in 1:length(countFiles)){
  countMat[[i]]<-read.table(countFiles[i], row.names=1)}

CountDat<-do.call(cbind, countMat)
colnames(CountDat)<-rownames(sampleTable)

#Remove the last fixe rows, which specify non-mapped reads
CountDat<-CountDat[1:(nrow(CountDat)-5),]

#Get Sample Blocks
sampBlocks<-sampleTable$ID

#Create model matrix from 3 covariates, ie sex, lane, condition
design<-model.matrix(~ Sex + Lane + condition, data = sampleTable)

#Create exon annotation object
gff<-read.table(flattenedFile, sep="\t", stringsAsFactors = FALSE)
head(gff)
colnames(gff)<-c("chr", "py", "type", "start", "end", "NA1", "strand", "NA2", "info")

####sub for exons
gff<-gff[gff$type=="exonic_part",]
####Remove type
gff<-gff[,colnames(gff)!="type"] 

#Get gene name
infoName<-unlist(strsplit(gff$info, ";"))
infoName<-infoName[seq(3, length(infoName), 3)]

infoName<-unlist(strsplit(infoName, " "))
infoName<-infoName[seq(3, length(infoName), 3)] 
gff$groupID<-infoName

#Get exon number
exonName<-unlist(strsplit(gff$info, ";"))
exonName<-exonName[seq(2, length(exonName), 3)]

exonName<-unlist(strsplit(exonName, " "))
exonName<-exonName[seq(3, length(exonName), 3)]
gff$featureID<-exonName
gff$exonID<-paste(gff$groupID, gff$featureID, sep=":")

gff<-gff[,c("chr", "start", "end", "strand", "groupID", "featureID", "exonID")]

#Subset and sort gff file by countDat rownames
gff<-gff[gff$exonID%in%rownames(CountDat),]
gff<-gff[order(match(gff$exonID, rownames(CountDat))),]

geneAnnot<-gff

#Create dge object
dge<-DGEList(CountDat, genes=geneAnnot)

#Keep exons(bins) with row counts greater than 10
A <- rowSums(dge$counts)
dge <- dge[A>10,,keep.lib.sizes=FALSE]

#Normalize counts
dge <- calcNormFactors(dge)

#Transform to log2 counts per million
v <- voom(dge,design,plot=TRUE)

#Fit correlations
corfit <- duplicateCorrelation(v, block = sampBlocks)

#Create mixed model
fit <- lmFit(v,design, block = sampBlocks, correlation = corfit$consensus.correlation)

# Perform differential spicing analysis
ex <- diffSplice(fit, geneid="groupID", exonid = "start")

save(fit, ex, design, sampBlocks, geneAnnot, sampleTable, file="LIMMA_splice_SexLane_mmID.RData")
