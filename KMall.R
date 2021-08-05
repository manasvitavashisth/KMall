rm(list=ls())
setwd("~/Box Sync/Lab/KM Plot/LiverCancer")
df.input1 = read.table("LIHC_clinicalMatrix", header = TRUE,sep="\t",stringsAsFactors = FALSE)
df.drugs = read.table("drugs.txt", header = TRUE,sep="\t",stringsAsFactors = FALSE)
df.map = read.table("LiverCancerAll", header = FALSE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
name=df.input[1]
namepro=df.map[1,2:dim(df.map)[2]]
Trk1 <- data.frame(matrix(ncol = 50, nrow = dim(df.map)[1]-1))
patientid=namepro[1]
y=1
k=1
for(i in 1:423)
{
  if(substr(namepro[i],14,15)=='01')
  {
    Trk1[,k]=as.double(df.map[2:dim(df.map)[1],i+1])
    patientid[k]=(namepro[1,i])
    k=k+1
  }
}

# If you want to match patients with their primary tumor and normal tissue respectively
y=1
for(i in 1:(dim(df.map)[2]-2))
{
  for(j in (i+1):dim(df.map)[2]-1)
  {
    if(substr(namepro[1,i],1,12)==substr(namepro[1,j],1,12))
    {
      if(substr(namepro[1,i],14,15)=='01' && substr(namepro[1,j],14,15)=='06')
      {
        patientid[y]=substr(namepro[1,i],1,12)
        Trk1[1:(dim(df.map)[1]-1),y]=as.double(as.double(df.map[2:dim(df.map)[1],i+1])-as.double(df.map[2:dim(df.map)[1],j+1]))
        
        y=y+1
      }else if(substr(namepro[1,j],14,15)=='01' && substr(namepro[1,i],14,15)=='06')
      {
        patientid[y]=substr(namepro[1,i],1,12)
        Trk1[1:(dim(df.map)[1]-1),y]=as.double(as.double(df.map[2:dim(df.map)[1],j+1])-as.double(df.map[2:dim(df.map)[1],i+1]))
      
        y=y+1
      }
    }
  }
}

#if you

Trk <- data.frame(matrix(ncol = 5, nrow = 371))
for(i in 1:dim(Trk)[1])
{
  for(j in 1:dim(df.input1)[1])
  {
    #if(name[j,1]==Trk[i,ncol(df.map)])
    if((name[j,1])==patientid[1,i] )
    {
      #Trk[i,26]=df.input$hist_hepato_carc_fact[j]
      #Trk[i,(ncol(df.map)+1)]=df.input$fibrosis_ishak_score[j]
      Trk[i,2]=df.input1$X_OS[j]/365
      Trk[i,3]=df.input1$X_OS_IND[j]
      #Trk2[i,nrow(df.map)+2]=df.input$hist_hepato_carc_fact[j]
      #Trk[i,(ncol(df.map)+2)]=df.input$X_OS_IND[j]
      #print(y)
      break
    }
  }
}
Trk2 <- data.frame(matrix(ncol = 4, nrow = 50))
for(i in 1:dim(Trk2)[1])
{
  for(j in 1:dim(df.drugs)[1])
  {
    #if(name[j,1]==Trk[i,ncol(df.map)])
    if(df.drugs$bcr_patient_barcode[j]==patientid[1,i] )
    {
      #Trk[i,26]=df.input$hist_hepato_carc_fact[j]
      #Trk[i,(ncol(df.map)+1)]=df.input$fibrosis_ishak_score[j]
      Trk2[i,2]=df.drugs$pharmaceutical_therapy_drug_name[j]
      Trk2[i,3]=df.drugs$pharmaceutical_therapy_type[j]
      #Trk2[i,4]=df.drugs$pharmaceutical_tx_started_days_to[j]
      #Trk2[i,nrow(df.map)+2]=df.input$hist_hepato_carc_fact[j]
      #Trk[i,(ncol(df.map)+2)]=df.input$X_OS_IND[j]
      #print(y)
      break
    }
  }
}
kmsig=(Trk1[1,2:4])
pro_name=df.map[1,1]
library(survminer)
library(survival)
k=1
for (i in 1:dim(Trk1)[1])
{
  if(sum(Trk1[i,]==0)<40)
  {
  Trk[,1]=as.double(t(Trk1[iden,]))
  Trk[,4]=NULL
  Trk[,4] <- Trk[,1]> median(Trk[,1])
  fit <- survfit(Surv(Trk[,2],Trk[,3]) ~ Trk[,4],
                 data = Trk)
  
    kmsig[k,1]=surv_pvalue(fit)$pval
    kmsig[k,2]=surv_median(fit)$median[1]
    kmsig[k,3]=surv_median(fit)$median[2]
    pro_name[k]=df.map[i+1,1]
    k=k+1
  }
}
kmsig1=kmsig
kmsig2=cbind(kmsig1,pro_name)
kmsig2[,5]=kmsig2[,3]/kmsig2[,2]
plot(-log2(as.double(kmsig2[,1])),log2(as.double(kmsig2[,5])),pch = 19,cex=0.3,col="grey",xlab='-log2(p-value)',ylab='-log2(median survival factor change)')

CNN1 <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
LMNB1_km=kmsig2[1,]
for (i in 1:dim(CNN1)[1])
{
  for (j in 1:length(pro_name))
  {
    if(pro_name[j]==trimws(CNN1[i,1]))
    {
      LMNB1_km[i,]=kmsig2[j,]
      #cou[i]=j
      break
    }
  }
}
#LMNB1_km=na.omit(LMNB1_km)
points(-log2(LMNB1_km[,1]),log2(LMNB1_km[,5]),pch = 19,cex=0.5,col="red")
text(-log2(LMNB1_km[,1]),log2(LMNB1_km[,5]), LMNB1_km[,4],
     cex=0.5, pos=3,col="red") 
ACTA <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
ACTA1=kmsig2[1,]
for (i in 1:dim(ACTA)[1])
{
  for (j in 1:length(pro_name))
  {
    if(pro_name[j]==trimws(ACTA[i,1]))
    {
      ACTA1[i,]=kmsig2[j,]
      #cou[i]=j
      break
    }
  }
}
#LMNB1_km=na.omit(LMNB1_km)
points(-log2(ACTA1[,1]),log2(ACTA1[,5]),pch = 19,cex=0.5,col="blue")
text(-log2(ACTA1[,1]),log2(ACTA1[,5]), ACTA1[,4],
     cex=0.5, pos=3,col="blue") 
Angio <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
Angiogenesis=kmsig2[1,]
for (i in 1:dim(Angio)[1])
{
  for (j in 1:length(pro_name))
  {
    if(pro_name[j]==trimws(Angio[i,1]))
    {
      Angiogenesis[i,]=kmsig2[j,]
      #cou[i]=j
      break
    }
  }
}
#LMNB1_km=na.omit(LMNB1_km)
points(-log2(Angiogenesis[,1]),log2(Angiogenesis[,5]),pch = 19,cex=0.5,col="black")
text(-log2(Angiogenesis[,1]),log2(Angiogenesis[,5]), Angiogenesis[,4],
     cex=0.5, pos=3,col="black") 
points(c(-log2(0.05),-log2(0.05),-log2(0.05)),-1:1,pch = 3,cex=1)
Angio <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
Angiogenesis=kmsig2[1,]
for (i in 1:dim(Angio)[1])
{
  for (j in 1:length(pro_name))
  {
    if(pro_name[j]==trimws(Angio[i,1]))
    {
      Angiogenesis[i,]=kmsig2[j,]
      #cou[i]=j
      break
    }
  }
}
#LMNB1_km=na.omit(LMNB1_km)
points(-log2(Angiogenesis[,1]),log2(Angiogenesis[,5]),pch = 19,cex=0.5,col="chartreuse3")
text(-log2(Angiogenesis[,1]),log2(Angiogenesis[,5]), Angiogenesis[,4],
     cex=0.5, pos=3,col="chartreuse3") 

sortedarray=kmsig2[1,]
k=1
for (j in 1:16733)
{
  if(kmsig2[j,5]>1)
  {
    sortedarray[k,]=kmsig2[j,]
    k=k+1
    #cou[i]=j
    #break
  }
}
sortedarray <- sortedarray[order(sortedarray[,1]),]
for (j in 1:dim(sortedarray)[1])
{
  if(sortedarray[j,1]>0.05)
  {
    #sortedarray[k,]=kmsig2[j,]
    iden=j
    #cou[i]=j
    break
  }
}
write.table(sortedarray, file = "FC_HIGH", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

CNN1 <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
ts_km=kmsig2[1,]
for (i in 1:dim(CNN1)[1])
{
  for (j in 1:15184)
  {
    if(pro_name[j]==trimws(CNN1[i,1]))
    {
      ts_km[i,]=kmsig2[j,]
      #cou[i]=j
      break
    }
  }
}
points(-log2(ts_km[,1]),log2(ts_km[,5]),pch = 19,cex=0.5,col="black")
text(-log2(ts_km[,1]),log2(ts_km[,5]), ts_km[,4],
     cex=0.5, pos=3,col="black") 
points(c(-log2(0.05),-log2(0.05),-log2(0.05)),-1:1,pch = 3,cex=1)

for (j in 1:16733)
{
  if(pro_name[j]=='CDK1')
  {
    iden=j
    #cou[i]=j
    break
  }
}
kmsig2[iden,]
a=1
for(i in 1:29)
{
for (j in 2:20531)
{
  if(df.map[j,1]==CNN1[i,1])
  {
    iden=j-1
    a[i]=median(as.double(Trk1[iden,]))
    #cou[i]=j
    break
  }
}
}
for (j in 2:20531)
{
  if(df.map[j,1]=='MRVI1')
  {
    iden=j-1
    #cou[i]=j
    break
  }
}
diff(range(Trk1[iden,]))
median(as.double(Trk1[iden,]))
Trk1[iden,]

CNN1 <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
compare2 <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
compare1 <- read.table(pipe("pbpaste"), sep="\t", header = FALSE)
a=compare1[1:106,1]
k=1
for (i in 1:dim(CNN1)[1])
{
  for (j in 1:dim(compare2)[1])
  {
    if(trimws(CNN1[i,1])==trimws(compare2[j,1]))
    {
      a[k]=compare1[j,1]
      k=k+1
    }
  }
}
a=t(a)
