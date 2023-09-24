
##############################################################
Arg<-commandArgs(T)  #
TAB_file<-Arg[1]
window_size<-Arg[2]
step<-Arg[3]
howmany_snp_number<-Arg[4]
out_put<-Arg[5]
##############################################################

options(scipen=100)
DATA_1<-read.table(TAB_file, header=F, as.is=TRUE, quote="", comment.char="", sep="\t")
window_size<-as.numeric(window_size)+0
minChr=window_size
step<-as.numeric(step)+0
howmany_snp_number<-as.numeric(howmany_snp_number)+0
window_size_print<-window_size/1000000

sortlist <- order(DATA_1$V1, pmax(DATA_1$V1, DATA_1$V2)) 
DATA_1 <-DATA_1[order(DATA_1$V1, DATA_1$V2),] 
chromosome<-unique(DATA_1$V1)

on_off<-0

for(chr_num in chromosome){
	DATA_1_M<-DATA_1[DATA_1$V1==chr_num,]
	max_position_of_chr<-max(DATA_1_M$V2)
	if(max_position_of_chr<minChr){
	  cat('sliding ',chr_num," skipped\r")
	  next
	}else{
	  cat ('sliding ',chr_num,'...')

  	max_value<-(max_position_of_chr-window_size/2)
  	window_start<-seq(window_size/2*-1,max_value, by=step)
  	
  	for(key_a in window_start){
  		DATA_1_M_AV<-DATA_1_M[DATA_1_M$V2>=key_a & DATA_1_M$V2<=key_a+window_size,]
  		position<-c(key_a+window_size/2)
  		nrow_number<-nrow(DATA_1_M_AV)
  		if(nrow_number > howmany_snp_number){
  			mean_SNP_index_a<-mean(DATA_1_M_AV$V3)
  			mean_SNP_index_b<-mean(DATA_1_M_AV$V4)
  			mean_SNP_index_D<-mean(DATA_1_M_AV$V5)
  			mean_95_l<-mean(DATA_1_M_AV$V6)
  			mean_95_u<-mean(DATA_1_M_AV$V7)
  			mean_99_l<-mean(DATA_1_M_AV$V8)
  			mean_99_u<-mean(DATA_1_M_AV$V9)
  			print_vector<-c(chr_num, position, mean_SNP_index_a, mean_SNP_index_b, mean_SNP_index_D,mean_95_l,mean_95_u,mean_99_l,mean_99_u)		
  			if(on_off==0){
  				write(t(print_vector), out_put, append=F,sep="\t", ncolumns=length(print_vector))
  				on_off=1
  			}else{
  				write(t(print_vector), out_put, append=T,sep="\t", ncolumns=length(print_vector))
  			}
  		}
  	}
    cat(" done\n")
	}
}

cat("\n")



