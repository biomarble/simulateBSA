args=commandArgs(T)

# args: targetChr,targetSize,outfile

library(tidyverse)

#  simulate SNP-index result for BSA

lenfile='len.txt'               # chrid  length
simfile="delta.simulation.txt"  # simulation of L95 H95 L99 H99 for each depth
tChr=args[1]                    # target QTL Chr
tSize=as.numeric(args[2])       # target QTL Size
outfile=args[3]                 # output file name

minDP=10
maxDP=100
meanDP=40

lensize=read.delim(lenfile,header=F)[,1:2]
simdat=read.delim(simfile)
colnames(lensize)=c('chr','size')
nchr=nrow(lensize)

data=NULL
for (chr in unique(lensize$chr)) {
  chrsize=lensize[lensize$chr==chr,'size']
  snpNumbers=floor(runif(1,min=0.00005,max=0.0001)*chrsize)   #simulate SNP number
  pos=unique(sort(floor(runif(snpNumbers,min=1,max=chrsize)))) #simulate SNP position uniformly  across the chromosome
  snpNumbers=length(pos)                                
  DP=floor(rnorm(snpNumbers,mean=meanDP,sd=20))                 #simulate seq depth for each SNP 
  DP[DP>maxDP]=maxDP
  DP[DP<minDP]=minDP
  thrdat=simdat[match(DP,simdat$DEPTH),]%>%select(-DEPTH)    
  rownames(thrdat)=pos
  
  SNPindex1=rnorm(length(pos),mean = 0.5,sd = 0.15)             # use normal distribution to generate  SNP index for mutation pool
  SNPindex1[SNPindex1>1]=1
  SNPindex1[SNPindex1<0]=0
  names(SNPindex1)=pos
  SNPindex2=rnorm(length(pos),mean = 0.5,sd = 0.15)             # use normal distribution to generate SNP index for  wildtype pool
  SNPindex2[SNPindex2>1]=1
  SNPindex2[SNPindex2<0]=0
  names(SNPindex2)=pos
  
  if(chr == tChr){                                              
    tPos=floor(runif(1,min = 0,max=chrsize))                    # generate  random QTL position 
    start=tPos-tSize/2                                          # define linkage disequilibrium start 
    end=tPos+tSize/2                                            # define linkage disequilibrium end
    for (i in pos[pos>start & pos<end]) {                   
      targetIndex=1-abs(i-tPos)/tSize                           # more colser to target, more larger SNPindex1
      targetIndex1=rnorm(1,mean=targetIndex,sd=0.3)             # use normal distribution to generate final value around the desired SNPindex1
      targetIndex1=min(targetIndex1,1)
      targetIndex1=max(targetIndex1,0)
      SNPindex1[names(SNPindex1)==i]=targetIndex1
    }
  }

  df=data.frame(chr,pos,SNPindex1,SNPindex2)%>%mutate(deltaSNPindex=SNPindex1-SNPindex2)
  df=bind_cols(df,thrdat)
  data=bind_rows(data,df)
}

write.table(data,file=outfile,row.names = F,col.names = F,quote=F,sep="\t")
