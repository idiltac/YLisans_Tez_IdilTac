setwd("/mnt/NEOGENE4/projects/medical_2020/ancient/project_3501/mathieson2015_freqs/freqs/Ancient/analysis1")
library(tidyverse)
library(ggpubr)
library(data.table)
library(dplyr)
library(Hmisc)
library(glue)

df = read.csv("filename", header = T, sep = ",") #example dataset 
df <- subset (df, select = c(-X))
head(df)

df$T <- as.numeric(df$T) #total reads
df$R <- as.numeric(df$R) #risk allele containing reads

#check data
class(df$T);class(df$R);class(df$X);class(df)

df_subset <-  split(df,list(df$Group,df$SNP_ID,df$Phenotype))
df_subset_v1 = df_subset[sapply(df_subset, function(x) dim(x)[1]) > 0]
df_subset_renamed <-setNames(df_subset_v1, as.vector(1:length(df_subset_v1)))


### compute frequency

#r=number of risk allele copies
#t=number of total allele copies
#p=allele frequency -> we would like to find this
#error rate epsilon -> 0.01 

freq_loop <-  function(dataframe) {
  function_for_freq <- function( r  ,t ,p) { 
    (((p^2*dbinom(r,t,1-0.01,log=T)) + ( 2*p*(1-p)*dbinom(r,t,0.5,log=T) ) + ( (1-p)^2*dbinom(r,t,0.01,log=T))))
  } 
  
ps=seq(0,1,by=0.01)
logLike=matrix(NA,length(ps),1)
rownames(logLike)=ps
for (p in ps) { 
    x=0
    for (i in 1:nrow(dataframe)) {
      x=x+(function_for_freq(dataframe$R[i],dataframe$T[i],p))
      logLike[rownames(logLike)==p,]=x} 
}
period_name <- as.character(unique(dataframe$Group))
snp_name <- as.character(unique(dataframe$SNP_ID))
phenotype_name <- as.character(unique(dataframe$Phenotype))
pop_size <- as.numeric(nrow(dataframe))
dlogLike <- logLike - max(logLike)
pHat <- ps[logLike == max(logLike)] ; pHat #Maximum likelihood estimate
pHat_max <- as.numeric(pHat)

#Confidence interval
#An approximate likelihood-based 95% confidence interval is obtained as follows.
x= pHat_max*pop_size/pop_size;  alpha <- .05
CI= x + c(-qnorm(1-alpha/2), qnorm(1-alpha/2))*sqrt((1/100)*x*(1-x)); CI = as.data.frame(CI)
#put the result into a data frame
likeResults <- data.frame(frequency = ps, diffLogLike = dlogLike, Loglikelihood = logLike,  Period= period_name, SNPid= snp_name, Phenotype=phenotype_name, POPsize= pop_size ,  pHat= pHat_max, CI_Lower= CI[1,] ,   CI_Upper = CI[2,])
likeResults$Loglikelihood <- as.numeric(likeResults$Loglikelihood)
maxlikelihood_result <- likeResults[which.max(likeResults$Loglikelihood),] 
print(maxlikelihood_result)
}

df_empty = data.frame()

for (i in 1:length(df_subset_renamed)) { 
  output = freq_loop(df_subset_renamed[[i]])
  df_empty = rbind(df_empty, output)
  maxll_results <- df_empty
}

maxll_results

write.csv(maxll_results, "filename.csv")

