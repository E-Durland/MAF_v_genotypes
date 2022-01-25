# MAF_v_genotypes
A re-analysis of genotype data from Plough &amp; Hedgecock (2011) in the framework of MAF changes during larval oyster development.
"1/25/2022"

This script takes genotype data supplied in the supplementary information for [Plough & Hedgecock (2011)](https://doi.org/10.1534/genetics.111.131854) and re-analyzes it within a framework of minor allele frequency (MAF) changes over larval oyster development.

*note:  The markers used in this analysis are a subset of the larger marker set used in that study. Filtering criteria was:
1) Bi-allelic markers
2) >=3 time points during larval development (0-45 days post fertilization)

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r}
#load the data 
df<- read.table("/PATH/TO/FILE/Plough_2011_bi_allele.txt",
                header=TRUE)
head(df)

#now, calculate the frequency of each genotype:
df$pAA <- apply(df,1,function(x){
  pAA <- as.numeric(x[5])
  pAB <- as.numeric(x[6])
  pBB <- as.numeric(x[7])
  pAA <- pAA/sum(pAA,pAB,pBB)
})
df$pAB <- apply(df,1,function(x){
  pAA <- as.numeric(x[5])
  pAB <- as.numeric(x[6])
  pBB <- as.numeric(x[7])
  pAA <- pAB/sum(pAA,pAB,pBB)
})
df$pBB <- apply(df,1,function(x){
  pAA <- as.numeric(x[5])
  pAB <- as.numeric(x[6])
  pBB <- as.numeric(x[7])
  pAA <- pBB/sum(pAA,pAB,pBB)
})

#calculate frequencies of each allele:
df$pA <- apply(df, 1, function(x){
  pAA <- as.numeric(x[5])
  pAB <- as.numeric(x[6])
  pBB <- as.numeric(x[7])
  pA <- (pAA*2+pAB)/(sum(pAA,pAB,pBB)*2)
  #pB <- (pAB + pBB*2)/(sum(pAA,pAB,pBB)*2)
})

df$pB <- apply(df, 1, function(x){
  pAA <- as.numeric(x[5])
  pAB <- as.numeric(x[6])
  pBB <- as.numeric(x[7])
  #pA <- (pAA*2+pAB)/(sum(pAA,pAB,pBB)*2)
  pB <- (pAB + pBB*2)/(sum(pAA,pAB,pBB)*2)
})

#reduce dataset to larval stages:
df <- df[df$day<45,]

#calculate minor allele frequency (MAF)
df$MAF <- apply(df,1,function(x){
  alleles <-c(as.numeric(x[11]),as.numeric(x[12]))
  MAF <- min(alleles/sum(alleles))
})

#calculate change (delta, or 'd' in this case) in MAF:
d_MAF<-c()
old_marker <- "x"
for(a in 1:nrow(df)){
  new_marker <- df$Marker[a]
  if(new_marker != old_marker){
    d_af <- NA
  }else{
    d_af <- df$MAF[a]-df$MAF[a-1]
  }
  old_marker <- new_marker
  d_MAF<- append(d_MAF,d_af)
}
df$dMAF <- d_MAF
head(df)

gtp <- pivot_longer(df[,c(1:10)],
                       cols=pAA:pBB,
                       names_to = "gtp",
                       values_to = "freq")%>%
  as.data.frame()
frq <- pivot_longer(df[,c(1:4,11,12)],
                       cols=pA:pB,
                       names_to = "allele",
                       values_to = "freq")%>%
  as.data.frame()

#genotype frequency
ggplot(gtp,aes(day,freq))+
  geom_line(aes(color=gtp),size=1.5)+
  facet_wrap(~Marker)+
  labs(title="Bi-allelic markers, genotype frequency")+
  theme_bw()
#frequency of both alleles
ggplot(frq,aes(day,freq))+
  geom_line(aes(color=allele),size=1.5)+
  facet_wrap(~Marker)+
  labs(title ="Bi-allelic markers, allele frequency")+
  theme_bw()

#Minor allele frequency acros time:
ggplot(df,aes(day,MAF))+
  geom_line(aes(color=Marker),size=1.5)+
  facet_wrap(~Marker)+
  labs(title="Minor allele frequency (Bi-allelic)")+
  guides(color=FALSE)+
  theme_bw()
#change in minor allele frequency:
ggplot(df,aes(day,dMAF))+
  geom_line(aes(color=Marker),size=1.5)+
  geom_hline(yintercept = 0, lty=2)+
  #scale_color_continuous(low="red",high="blue")+
  labs(title=expression(paste(Delta, "Minor allele frequency (Bi-allelic)")))+
  facet_wrap(~Marker)+
  guides(color=FALSE)+
  theme_bw()
```
