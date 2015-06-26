fileName <- list.files(path = "/Users/neemoyaz/Google Drive/Dr. Camron Bryant/Dr. Evan Johnson/HTSeq & EdgeR output/", pattern="txt", full.names = T, recursive = T)

expr <- NULL
sampleName <- NULL
for (i in fileName[-7]){
  expr <- cbind(expr,read.table(i)[,2])
  sampleName <- c(sampleName, substr(i,4,nchar(i)-4)) 
}
 

rownames(expr) <- read.table(fileName[1])[,1]
colnames(expr) <- sampleName

expr1 <- expr[1:(nrow(expr)-7),]
expr2 <- expr1[,c(14:31,1:13)]

write.table(expr2,file="htseq_exprMatrix.txt", quote=F, sep="\t")