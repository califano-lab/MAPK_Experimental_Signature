#### Tutorial for the Generation of the MAPK Consensus Gene Expression Signature from the PlateSeq Drugs Screening

### 1) Load the logTPM-normalized expression data of both the cell lines (ASPC1 and PANC1) included in the screening
ASPC1_logTPM_Expression<-readRDS("ASPC1_logTPM_Expression_Final.rds")
PANC1_logTPM_Expression<-readRDS("PANC1_logTPM_Expression_Final.rds")

### 2) Load the metadata files with all the information about the treatment conditions of each well in the 4 Plates
ULAASPC1PSA_Final<-readRDS("ULAASPC1PSA_Final.rds")
ULAASPC1PSB_Final<-readRDS("ULAASPC1PSB_Final.rds")
ULAPANC1PSA_Final<-readRDS("ULAPANC1PSA_Final.rds")
ULAPANC1PSB_Final<-readRDS("ULAPANC1PSB_Final.rds")

### 3) For each one of the two cell lines, retrieve the column names corresponding to the DMSO wells, to be used
###    as a reference for the computation of the Gene Expression Signatures
DMSO_ASPC1<-c(ULAASPC1PSA_Final$colNames[which(ULAASPC1PSA_Final$`Compound 1`=="DMSO")],
              ULAASPC1PSB_Final$colNames[which(ULAASPC1PSB_Final$`Compound 1`=="DMSO")])

DMSO_PANC1<-c(ULAPANC1PSA_Final$colNames[which(ULAPANC1PSA_Final$`Compound 1`=="DMSO")],
              ULAPANC1PSB_Final$colNames[which(ULAPANC1PSB_Final$`Compound 1`=="DMSO")])

### 4) For each one of the two cell lines, compute z-Scores with respect to DMSO
mean_DMSO_ASPC1<-apply(ASPC1_logTPM_Expression[,DMSO_ASPC1],1,mean)
mean_DMSO_PANC1<-apply(PANC1_logTPM_Expression[,DMSO_PANC1],1,mean)

sd_DMSO_ASPC1<-apply(ASPC1_logTPM_Expression[,DMSO_ASPC1],1,sd)
sd_DMSO_PANC1<-apply(PANC1_logTPM_Expression[,DMSO_PANC1],1,sd)

signature_ASPC1_DMSORef<-(ASPC1_logTPM_Expression - mean_DMSO_ASPC1)/sd_DMSO_ASPC1
signature_ASPC1_DMSORef<-signature_ASPC1_DMSORef[which(sd_DMSO_ASPC1>0),]

signature_PANC1_DMSORef<-(PANC1_logTPM_Expression - mean_DMSO_PANC1)/sd_DMSO_PANC1
signature_PANC1_DMSORef<-signature_PANC1_DMSORef[which(sd_DMSO_PANC1>0),]


### 5) From both Signatures, remove columns corresponding to DMSOs and integrate the rest of the treatment conditions
####   
Drugs_ASPC1<-c(ULAASPC1PSA_Final$colNames[which(ULAASPC1PSA_Final$`Compound 1`!="DMSO")],
               ULAASPC1PSB_Final$colNames[which(ULAASPC1PSB_Final$`Compound 1`!="DMSO")])

Drugs_PANC1<-c(ULAPANC1PSA_Final$colNames[which(ULAPANC1PSA_Final$`Compound 1`!="DMSO")],
               ULAPANC1PSB_Final$colNames[which(ULAPANC1PSB_Final$`Compound 1`!="DMSO")])

signature_ASPC1_DMSORef_Integrated<-rowSums(signature_ASPC1_DMSORef[,Drugs_ASPC1])/sqrt(length(Drugs_ASPC1))
signature_PANC1_DMSORef_Integrated<-rowSums(signature_PANC1_DMSORef[,Drugs_PANC1])/sqrt(length(Drugs_PANC1))

common_genes<-intersect(names(signature_ASPC1_DMSORef_Integrated), 
                        names(signature_PANC1_DMSORef_Integrated))

consesus_sign_MAPK_GEXP<- - (signature_ASPC1_DMSORef_Integrated[common_genes] + 
                             signature_PANC1_DMSORef_Integrated[common_genes])/2

