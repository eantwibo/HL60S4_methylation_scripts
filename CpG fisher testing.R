library(GetoptLong)

GetoptLong(c("filename=s", "file name"))

temp_dataset <-read.table(filename, header=FALSE)

inputseq <-temp_dataset[complete.cases(temp_dataset),]
chrm_name <- unique(inputseq[,1])

##############################
#alpha = 0.05
#########using methylation value and coverage for fischer exact test for significance
RA_control_data<-cbind(inputseq[,7:8],inputseq[,10:11])
RA_control_dataframe <-as.data.frame(RA_control_data)
RA_control_pvalues<-matrix(0,length(RA_control_dataframe[,1]),1)

TPA_control_data<-cbind(inputseq[,7],inputseq[,9],inputseq[,10],inputseq[,12])
TPA_control_dataframe <-as.data.frame(TPA_control_data)

RA_TPA_data <- cbind(inputseq[,8:9],inputseq[,11:12])
RA_TPA_dataframe <-as.data.frame(RA_TPA_data)


all_pvalues<-matrix(0,length(inputseq[,1]),3)
for(i in 1:length(inputseq[,1])){
RA <- fisher.test(matrix(as.numeric(RA_control_dataframe[i,]),2,2))$p.value
TPA <- fisher.test(matrix(as.numeric(TPA_control_dataframe[i,]),2,2))$p.value
RA_TPA <- fisher.test(matrix(as.numeric(RA_TPA_dataframe[i,]),2,2))$p.value
#ifelse(i==1, RA_control_pvalues <-g,RA_control_pvalues <-c(RA_control_pvalues,g))
all_values <- c(RA,TPA,RA_TPA)
all_pvalues[i,] <- all_values
}
#ifelse(i==1, TPA_control_pvalues <-g,TPA_control_pvalues<-c(TPA_control_pvalues,g))

all_data <- cbind(inputseq[,],all_pvalues)

write.table(all_data, paste("Fischer_tested_",chrm_name,"_raw_CpG.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
