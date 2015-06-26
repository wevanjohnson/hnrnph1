

run4all<- function(trt){
#########################
# load genetic matrices #
#########################

load("idcfold.RData")

ibF2<- gmF2$AA/2
AAF2<- gmF2$AA
DDF2<- gmF2$DD
ADF2<- gmF2$AD
HHF2<- gmF2$HH
MHF2<- gmF2$MH

ibF8<- gmF8$AA/2
AAF8<- gmF8$AA
DDF8<- gmF8$DD
ADF8<- gmF8$AD
HHF8<- gmF8$HH
MHF8<- gmF8$MH

####################################
# load genotype and phenotype data #
# pedigree info and more           #
####################################

load("imFnold.RData")

idF2<- dimnames(prDatF2$pr)[[1]]
idF8<- dimnames(prDatF8$pr)[[1]]

pdatF2<- read.table(file="MAF2.txt",check.names=F,sep="\t", header=T, na.string="?")
   pdatF2<- pdatF2[!is.na(pdatF2[,trt]),]
ii<- match(idF2,pdatF2$"SUBJECT ID")
   sum(is.na(ii)) # 3
idF2<- idF2[!is.na(ii)]
ii<- ii[!is.na(ii)]
pdatF2<- pdatF2[ii,]

ii<- match(idF2,colnames(AAF2))
AAF2<- AAF2[ii,ii]
DDF2<- DDF2[ii,ii]
ADF2<- ADF2[ii,ii]
HHF2<- HHF2[ii,ii]
MHF2<- MHF2[ii,ii]
ibF2<- ibF2[ii]
ii<- match(idF2,dimnames(prDatF2$pr)[[1]])
prDatF2$pr<- prDatF2$pr[ii,,]

# hist(pdatF2[,trt],breaks=40,)

pdatF8<- read.table(file="MAF8.txt",check.names=F, sep="\t",header=T, na.string="?")
   pdatF8<- pdatF8[!is.na(pdatF8[,trt]),]
ii<- match(idF8,pdatF8$"SUBJECT ID")
   sum(is.na(ii)) # 3
idF8<- idF8[!is.na(ii)]
ii<- ii[!is.na(ii)]
pdatF8<- pdatF8[ii,]

ii<- match(idF8,colnames(AAF8))
AAF8<- AAF8[ii,ii]
DDF8<- DDF8[ii,ii]
ADF8<- ADF8[ii,ii]
HHF8<- HHF8[ii,ii]
MHF8<- MHF8[ii,ii]
ibF8<- ibF8[ii]
ii<- match(idF8,dimnames(prDatF8$pr)[[1]])
prDatF8$pr<- prDatF8$pr[ii,,]

##############
# condl prob #
##############

prDat<- array(0,dim=c(dim(prDatF2$pr)[1]+dim(prDatF8$pr)[1],3,dim(prDatF2$pr)[3]))
   prDat[,1,]<- rbind(prDatF2$pr[,1,],prDatF8$pr[,1,])
   prDat[,2,]<- rbind(prDatF2$pr[,2,],prDatF8$pr[,2,])
   prDat[,3,]<- rbind(prDatF2$pr[,3,],prDatF8$pr[,3,])
   prDat<- list(pr=prDat,chr=prDatF2$chr,dist=prDatF2$dist,snp=prDatF2$snp)
any(is.na(prDat$pr))
dim(prDat$pr)
c(length(prDat$chr),length(prDat$dist))

########################
# trait and covariates #
########################

y<- c(scale(pdatF2[,trt]),scale(pdatF8[,trt])) # **************
   y<- as.matrix(y,ncol=1)
x<- c(pdatF2[,c("Sex")],pdatF8[,c("Sex")])
   x<- as.matrix(x)
intcv<- NULL

isF8<- c(rep(F,length(idF2)),rep(T,length(idF8)))

#################
# model fitting #
#################

ovF2<- list(AA=AAF2,DD=DDF2,AD=NULL,HH=NULL,MH=NULL,EE=diag(sum(!isF8)))
cat(date(),"\n")
  ooF2<- estVC(y[!isF8,],x[!isF8,],v=ovF2)
#cat(date(),"\n")

ovF8<- list(AA=AAF8,DD=DDF8,AD=NULL,HH=NULL,MH=NULL,EE=diag(sum(isF8)))
#cat(date(),"\n")
  ooF8<- estVC(y[isF8,],x[isF8,],v=ovF8)
cat(date(),"\n")

gcvF2<- matrix(0,nrow=nrow(ooF2$v$EE),ncol=nrow(ooF2$v$EE))
for(i in 1:ooF2$nv)
   if(ooF2$nnl[i]) gcvF2<- gcvF2 + ooF2$v[[i]]*ooF2$par[ncol(ooF2$x)+ooF2$nn[i]]
gcvF8<- matrix(0,nrow=nrow(ooF8$v$EE),ncol=nrow(ooF8$v$EE))
for(i in 1:ooF8$nv)
   if(ooF8$nnl[i]) gcvF8<- gcvF8 + ooF8$v[[i]]*ooF8$par[ncol(ooF8$x)+ooF8$nn[i]]
rm(i)

gcv<- matrix(0,nrow=nrow(y),nrow(y))
gcv[1:sum(!isF8),1:sum(!isF8)]<- gcvF2
   gcv[1:sum(isF8)+sum(!isF8),1:sum(isF8)+sum(!isF8)]<- gcvF8

cat(date(),"\n")
  llk<- scanOne(y=y,x=cbind(x,isF8),vc=gcv,intcovar=intcv,prdat=prDat)
#cat(date(),"\n")

#cat(date(),"\n")
  llkF2<- scanOne(y=y[!isF8,],x=x[!isF8,],vc=gcvF2,intcovar=intcv[!isF8],prdat=prDatF2)
#cat(date(),"\n")

#cat(date(),"\n")
  llkF8<- scanOne(y=y[isF8,],x=x[isF8,],vc=gcvF8,intcovar=intcv[isF8],prdat=prDatF8)
cat(date(),"\n")

save(list=union(ls(all=TRUE,pos=1),ls(all=TRUE)),file=paste("R Workspace Images/ma",trt,".RData",sep="")) # **********


# q("no")

}

if(F){
################
# thresholds #
##############

nn<- length(isF8)
ntimes<- 1000
cvMtr<- NULL
for(n in 1:ntimes){
  cat(n,"/",ntimes,"\r"); flush.console() # track process
  idx<- sample(1:nn,replace=FALSE)
  prdTmp<- prDat
  prdTmp$pr<- prdTmp$pr[idx,,]
  tmp<- scanOne(y=y,x=cbind(x,isF8),vc=gcv,intcovar=intcv,prdat=prdTmp)
  cvMtr<- rbind(cvMtr,tmp$p)
}

save(list=union(ls(all=TRUE,pos=1),ls(all=TRUE)),file=paste("R Workspace Images/ma",trt,".RData",sep="")) # **********

################
# mapping plot #
################

cv<- NULL
llkB<- rbind(data.frame(p=llk$p,chr=llk$chr,dist=llk$dist),
               data.frame(p=llkF2$p,chr=llkF2$chr,dist=llkF2$dist),
               data.frame(p=llkF8$p,chr=llkF8$chr,dist=llkF8$dist))
   llkB$group<- c(rep("Integrated",length(llk$p)),rep("F2",length(llkF2$p)),rep("F8",length(llkF8$p)))
   llkB$y<- llkB$p/(2*log(10))
   llkB$chr<- reorder(llkB$chr)
   llkB$group<- reorder(as.factor(llkB$group))
   levels(llkB$group)
#postscript(paste("ma",trt,"_map.ps",sep=""),horizontal=F)
pdf(paste("ma",trt,"_map.pdf",sep=""))
   tmp<- plotit(llkB,cv=cv,type="l",lty=c(2,4,5),col=c(2,3,4),xlab="Genetic Position (cM)",ylab="LOD",main=trt,bychr=T)
   update(tmp,key=list(lines=list(lty=c(2,4,5),col=c(2,3,4)),
                       text=list(c(expression(F[2]),"AIL","Integrated")),
                       columns=3))
   

   plotit(data.frame(y=llkF2$p/(2*log(10)),chr=llkF2$chr,dist=llkF2$dist),cex=0.5,main="F2",xlab="Genetic Position",ylab="LOD",ylim=c(0,15))
   plotit(data.frame(y=llkF8$p/(2*log(10)),chr=llkF8$chr,dist=llkF8$dist),cex=0.5,main="F8",xlab="Genetic Position",ylab="LOD",ylim=c(0,15))
   plotit(data.frame(y=llk$p/(2*log(10)),chr=llk$chr,dist=llk$dist),cex=0.5,main="Integrated",xlab="Genetic Position",ylab="LOD",ylim=c(0,15))
dev.off()

llkAll<- data.frame(snp=prDat$snp,
                  chr=llk$chr,
                  genPos=llk$dist,
                  lr=llk$p,
                  lrF2=llkF2$p,
                  lrF8=llkF8$p)
est<- matrix(unlist(llk$parameters),nrow=length(llk$p),byrow=T)
#   rownames(est)<- names(llk$parameters)
   colnames(est)<- names(llk$parameters[[1]])

llkEst<- cbind(llkAll,est)

write.csv(llkEst,file=paste("llkma",trt,".csv",sep=""),row.names=F)

}

#########
# start #
#########

setwd("C:/Users/Palmer Lab/Documents/Clarissa Documents/B6 x D2 Data/Mapping Script")

library(QTLRel)

library(QTLRel)
library(latticedl)

allTraits<- c("D1TOTDIST5","D1TOTDIST10","D1TOTDIST15","D1TOTDIST20","D1TOTDIST25","D1TOTDIST30",
"D2TOTDIST5","D2TOTDIST10","D2TOTDIST15","D2TOTDIST20","D2TOTDIST25","D2TOTDIST30",
"D3TOTDIST5","D3TOTDIST10","D3TOTDIST15","D3TOTDIST20","D3TOTDIST25","D3TOTDIST30")

for(n in 1:length(allTraits)){
   trtTmp<- allTraits[n]
   print(trtTmp)
   run4all(trtTmp)
}


########################################################################
# the end #
###########

setwd("C:/Users/Palmer Lab/Documents/Clarissa Documents/B6 x D2 Data/Mapping Script")
library(QTLRel)

ntimesTmp<- 1000
allTraitsTmp<- c("D1TOTDIST5","D1TOTDIST10","D1TOTDIST15","D1TOTDIST20","D1TOTDIST25","D1TOTDIST30",
"D2TOTDIST5","D2TOTDIST10","D2TOTDIST15","D2TOTDIST20","D2TOTDIST25","D2TOTDIST30",
"D3TOTDIST5","D3TOTDIST10","D3TOTDIST15","D3TOTDIST20","D3TOTDIST25","D3TOTDIST30")
cv4All<- NULL
objTmp<- ls()
for(n in 1:length(allTraitsTmp)){
   cvMtr<- NULL
   trtTmp<- allTraitsTmp[n]
   print(trtTmp)
   load(paste("R Workspace Images/ma",trtTmp,".RData",sep=""))
   nn<- length(isF8)
   for(n in 1:ntimesTmp){
     idx<- sample(1:nn,replace=FALSE)
     prdTmp<- prDat
     prdTmp$pr<- prdTmp$pr[idx,,]
     tmp<- scanOne(y=y,x=cbind(x,isF8),vc=gcv,intcovar=intcv,prdat=prdTmp)
     cvMtr<- rbind(cvMtr,tmp$p)
     cat(n,"/",ntimesTmp,"\r") # track process
   }
   cv4All<- cbind(cv4All,apply(cvMtr,1,max))
   rm(list=setdiff(ls(),c(objTmp,"objTmp")))
}
colnames(cv4All)<- allTraitsTmp

# genome-wide (LOD) thresholds at 0.1, 0.05 and 0.01
 dim(cv4All)
 apply(cv4All,2,quantile,c(0.9,0.95,0.99))/(2*log(10))








##########################
# output: combine output #
##########################

# create output table
llkOutput<- function(llk,y,gdat,gmap){
# llk: from llr()
# gat: genotype data used in mapping (snp in column)
# gmap: genetic map
   fst<- function(y,gdat,grp){
      y<- as.matrix(y)
      gd<- gdat==grp
      out<- data.frame(n=rep(NA,ncol(gdat)),mean=rep(NA,ncol(gdat)),sd=rep(NA,ncol(gdat)))
      for(j in 1:ncol(gdat)){
         idx<- gd[,j]
         if(any(idx)){
            y0<- y[idx,1]
            out$n[j]<- sum(idx)
            out$mean[j]<- mean(y0)
            if(sum(idx)>1) out$sd[j]<- sd(y0)
         }
      }
      out
   }

   snp<- llk$snp
   idx<- match(snp,gmap$snp)
   chr<- gmap$ch[idx]
   genPos<- gmap$dist[idx]
   idx<- match(snp,colnames(gdat))
   gdt<- as.matrix(gdat)[,idx]
      gdt[is.na(gdt)]<- 0
   tmp<- fst(y,gdt,"1")
   nAA<- tmp$n
   meanAA<- tmp$mean
   sdAA<- tmp$sd
   tmp<- fst(y,gdt,"2")
   nAB<- tmp$n
   meanAB<- tmp$mean
   sdAB<- tmp$sd
   tmp<- fst(y,gdt,"3")
   nBB<- tmp$n
   meanBB<- tmp$mean
   sdBB<- tmp$sd

   out<- data.frame(snp=llk$snp,
                    chr=chr,
                    genPos=genPos,
                    nAA=nAA,
                    meanAA=meanAA,
                    sdAA=sdAA,
                    nAB=nAB,
                    meanAB=meanAB,
                    sdAB=sdAB,
                    nBB=nBB,
                    meanBB=meanBB,
                    sdBB=sdBB,
                    lr=llk$p)
   out
}

gmap<- read.table("gmapFn.txt",check.names=F,header=T)
gdatF8<- read.table("gdatF8.txt",check.names=F,header=T)
   gdatF8<- gdatF8[match(rownames(prDatF8$pr),rownames(gdatF8)),]
llktb<- llkOutput(llk=llk,y=pdatF8[,trt],gdat=gdatF8,gmap=gmap)

tmp<- cbind(llktb,
            lrF2=llkF2$p,
            lrF8=llkF8$p)
   tmp$snp<- prDatF8$snp
   tmp$chr<- prDatF8$chr
   tmp$genPos<- prDatF8$dist
c(nrow(tmp),sum(!is.na(tmp$pF2)))

# write.csv(tmp,file="llkma.csv",row.names=F)

#######################################################################################
# the end # the end # the end #
###############################

