

source("pD_functions.R")
library(dplyr)
library(reshape2)
library(stringr)
library(data.table)
library(diann)

#CL7496 analysis
DA <- "/eJD1262_1267/Report.tsv"#path the report.tsv
batches_pagelab <- "batches_pagelab.tsv"

dat <- data.frame(fread(DA))
prot_info <- dat %>% dplyr::distinct(Protein.Group,.keep_all=T) %>% dplyr::select("Protein.Group","Genes")

#plex_MS1 <- dat[grepl("1262|1263|1264",dat$Run),] #select 7506 cell lines
plex_MS1 <- dat[grepl("1265|1266|1267",dat$Run),] #select 7496 cell lines
plex_MS1 <- plex_MS1[which(plex_MS1$Lib.PG.Q.Value<0.01),]  #1% protein FDR
plex_MS1 <- plex_MS1[plex_MS1$Translated.Q.Value<0.01,]
plex_MS1 <- plex_MS1 %>% dplyr::rename("Intensity" = "Ms1.Area")
plex_MS1 <- pD_channel(plex_MS1)
plex_MS1 <- pD_seqcharge(plex_MS1)
plex_MS1$precRun <- paste0(plex_MS1$seqcharge, plex_MS1$Run)
plex_MS1$chanRun <- paste0(plex_MS1$Run, "_",plex_MS1$channel_name)
plex_MS1$File.Name <- plex_MS1$chanRun

prots <- data.frame(diann_maxlfq(plex_MS1, id.header = "seqcharge", 
                                 quantity.header = "Intensity", group.header = "Protein.Group"))
 
#normalize each sample for loading amounts, then relative within a run
prots$PG <- row.names(prots)
prots_m <- melt(prots)
prots_m$Run <- sub('\\_.*', '', prots_m$variable)
prots_m$log2_value <- log2(prots_m$value)
prots_m <- prots_m  %>%
  dplyr::group_by(PG,Run) %>% dplyr::mutate(norm = log2_value-mean(log2_value,na.rm=T)) %>% ungroup() %>%
  dplyr::group_by(variable) %>% 
  dplyr::mutate(norm = norm-median(norm,na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::group_by(PG) %>% dplyr::mutate(norm = norm-mean(norm,na.rm=T)) %>% ungroup()

prots_m <- prots_m %>% left_join(prot_info, by =c("PG"="Protein.Group"))
bat <- read.delim(batches_pagelab)
bat$run_chan <- paste0(bat$Run, "_", bat$Label) 
bat <- bat %>% dplyr::select(-c("Run"))
prots_m <- prots_m %>% inner_join(bat, by=c("variable" = "run_chan"))
prots_m$Cond_lab <- paste0(prots_m$Condition,"_",prots_m$Label)
#plot cors
d_plexMS1 <- dcast(prots_m, PG~Cond_lab, value.var = "norm")
row.names(d_plexMS1) <- d_plexMS1$PG
#prot_infos <- row.names(d_plexMS1)
d_plexMS1 <- d_plexMS1[,-c(1)]

#### batch correction:
na_matrix <- is.na(d_plexMS1)

d_plexMS1 <- data.frame(hknn(as.matrix(d_plexMS1),k=3))

bat$Cond_lab <- paste0(bat$Condition, "_",meta$Label)

bat_groupings <- bat %>% dplyr::filter(grepl("1265|1266|1267",run_chan)) %>% distinct(Cond_lab,run_chan)

batch.covs <-bat$Label[match(colnames(d_plexMS1), bat$Cond_lab)]
matrix.sc.batch <- ComBat(d_plexMS1, batch=batch.covs)
matrix.sc.batch[na_matrix] <- NA
prots_m <- melt(matrix.sc.batch)
prots_m <- prots_m %>% left_join(bat_groupings, by =c("Var2" = "Cond_lab"))
prots_m$Run <- sub('\\_.*', '', prots_m$run_chan)

prots_m <- prots_m  %>%
  dplyr::group_by(Var1,Run) %>% dplyr::mutate(norm = value-mean(value,na.rm=T)) %>% ungroup() %>%
  dplyr::group_by(Var2) %>% 
  dplyr::mutate(norm = norm-median(norm,na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::group_by(Var1) %>% dplyr::mutate(norm = norm-mean(norm,na.rm=T)) %>% ungroup()

prots_m <- prots_m %>% left_join(prot_info, by =c("Var1" = "Protein.Group"))

write.table(prots_m, "pD_batchCorLabel_withinrunNorm_05172023_CL7496.txt", sep = "\t", row.names = FALSE)

#CL7506 analysis
DA <- "/eJD1262_1267/Report.tsv"#path the report.tsv
batches_pagelab <- "batches_pagelab.tsv"

dat <- data.frame(fread(DA))
prot_info <- dat %>% dplyr::distinct(Protein.Group,.keep_all=T) %>% dplyr::select("Protein.Group","Genes")

plex_MS1 <- dat[grepl("1262|1263|1264",dat$Run),] #select 7506 cell lines
#plex_MS1 <- dat[grepl("1265|1266|1267",dat$Run),] #select 7496 cell lines
plex_MS1 <- plex_MS1[which(plex_MS1$Lib.PG.Q.Value<0.01),]  #1% protein FDR
plex_MS1 <- plex_MS1[plex_MS1$Translated.Q.Value<0.01,]
plex_MS1 <- plex_MS1 %>% dplyr::rename("Intensity" = "Ms1.Area")
plex_MS1 <- pD_channel(plex_MS1)
plex_MS1 <- pD_seqcharge(plex_MS1)
plex_MS1$precRun <- paste0(plex_MS1$seqcharge, plex_MS1$Run)
plex_MS1$chanRun <- paste0(plex_MS1$Run, "_",plex_MS1$channel_name)
plex_MS1$File.Name <- plex_MS1$chanRun

prots <- data.frame(diann_maxlfq(plex_MS1, id.header = "seqcharge", 
                                 quantity.header = "Intensity", group.header = "Protein.Group"))
prots$PG <- row.names(prots)
######### normalize within a run and batch correct label:
prots_m <- melt(prots)
prots_m$Run <- sub('\\_.*', '', prots_m$variable)
prots_m$log2_value <- log2(prots_m$value)

#normalize each sample for loading amounts, then relative within a run
prots_m <- prots_m  %>%
  dplyr::group_by(PG,Run) %>% dplyr::mutate(norm = log2_value-mean(log2_value,na.rm=T)) %>% ungroup() %>%
  dplyr::group_by(variable) %>% 
  dplyr::mutate(norm = norm-median(norm,na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::group_by(PG) %>% dplyr::mutate(norm = norm-mean(norm,na.rm=T)) %>% ungroup()

prots_m <- prots_m %>% left_join(prot_info, by =c("PG"="Protein.Group"))
bat <- read.delim(batches_pagelab)
bat$run_chan <- paste0(bat$Run, "_", bat$Label) 
bat <- bat %>% dplyr::select(-c("Run"))
prots_m <- prots_m %>% inner_join(bat, by=c("variable" = "run_chan"))
prots_m$Cond_lab <- paste0(prots_m$Condition,"_",prots_m$Label)
d_plexMS1 <- dcast(prots_m, PG~Cond_lab, value.var = "norm")
row.names(d_plexMS1) <- d_plexMS1$PG
d_plexMS1 <- d_plexMS1[,-c(1)]

#### batch correction:
na_matrix <- is.na(d_plexMS1)
d_plexMS1 <- data.frame(hknn(as.matrix(d_plexMS1),k=3))
bat$Cond_lab <- paste0(bat$Condition, "_",meta$Label)
bat_groupings <- bat %>% dplyr::filter(grepl("1262|1263|1264",run_chan)) %>% distinct(Cond_lab,run_chan)

batch.covs <-bat$Label[match(colnames(d_plexMS1), bat$Cond_lab)]
matrix.sc.batch <- ComBat(d_plexMS1, batch=batch.covs)
matrix.sc.batch[na_matrix] <- NA
prots_m <- melt(matrix.sc.batch)
prots_m <- prots_m %>% left_join(bat_groupings, by =c("Var2" = "Cond_lab"))
prots_m$Run <- sub('\\_.*', '', prots_m$run_chan)

prots_m <- prots_m  %>%
  dplyr::group_by(Var1,Run) %>% dplyr::mutate(norm = value-mean(value,na.rm=T)) %>% 
  ungroup() %>%
  dplyr::group_by(Var2) %>% 
  dplyr::mutate(norm = norm-median(norm,na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::group_by(Var1) %>% dplyr::mutate(norm = norm-mean(norm,na.rm=T)) %>% ungroup()

prots_m <- prots_m %>% left_join(prot_info, by =c("Var1" = "Protein.Group"))

write.table(prots_m, "pD_batchCorLabel_withinrunNorm_05172023_CL7506.txt", sep = "\t", row.names = FALSE)


#############
#############
#############
#############

CL7506 <- read.delim("pD_batchCorLabel_withinrunNorm_05172023_CL7506.txt")
DX3XY_7506 <- CL7506[(CL7506$Genes=="DDX3X;DDX3Y"),]
DX3X_Y_7506 <- CL7506[(CL7506$Genes=="DDX3X"|CL7506$Genes=="DDX3Y"),]
CL7506$CL <- "CL7506"

CL7496 <- read.delim("pD_batchCorLabel_withinrunNorm_05172023_CL7496.txt")
DX3XY_7496 <- CL7496[(CL7496$Genes=="DDX3X;DDX3Y"),]
DX3X_Y_7496 <- CL7496[(CL7496$Genes=="DDX3X"|CL7496$Genes=="DDX3Y"),]
CL7496$CL <- "CL7496"

both_CL <- rbind(CL7496,CL7506)
########
both_CL <- both_CL %>% dplyr::mutate("Cond" = ifelse(grepl("Control", Var2), "Control",
                                                     ifelse(grepl("DDX3X_kd",Var2),"DDX3X KD", "DDX3Y KD")))
head(both_CL)
med_pD_control_DDX3X <- both_CL[grepl("Control",both_CL$Var2)&both_CL$Genes=="DDX3X",]
med_pD_control_DDX3X <- median(med_pD_control_DDX3X$norm,na.rm=T)
med_pD_DDX3X <- both_CL[both_CL$Genes=="DDX3X",]
med_pD_DDX3X$norm  <- med_pD_DDX3X$norm-med_pD_control_DDX3X

med_pD_control_DDX3Y <- both_CL[grepl("Control",both_CL$Var2)&both_CL$Genes=="DDX3Y",]
med_pD_control_DDX3Y <- median(med_pD_control_DDX3Y$norm,na.rm=T)
med_pD_DDX3Y <- both_CL[both_CL$Genes=="DDX3Y",]
med_pD_DDX3Y$norm <- med_pD_DDX3Y$norm-med_pD_control_DDX3Y

med_pD_control_DDX3Y_X <- both_CL[grepl("Control",both_CL$Var2)&both_CL$Genes=="DDX3X;DDX3Y",]
med_pD_control_DDX3Y_X <- median(med_pD_control_DDX3Y_X$norm,na.rm=T)
med_pD_DDX3Y_X <- both_CL[both_CL$Genes=="DDX3X;DDX3Y",]
med_pD_DDX3Y_X$norm <- med_pD_DDX3Y_X$norm-med_pD_control_DDX3Y_X

both <- rbind(med_pD_DDX3X,med_pD_DDX3Y,med_pD_DDX3Y_X)

ggplot(both, aes(x=Cond, y=norm, color=CL)) + geom_boxplot() + 
  geom_point(size=1,position=position_dodge(width=0.75)) +
  facet_grid(~Genes) + labs(x="Condition", y="Relative protein abundance, Log2", color="Cell line") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("DDX3X_DDX3Y_05162023.png",width=6,height=5)

###### DA prots

#ncol(both_CL)
both_CL_fin <- both_CL[!is.na(both_CL$norm),]
both_CL_fin <- both_CL_fin[both_CL_fin$norm!=0,] #to remove proteins which had no other sample to perform relative quant to (these would have quant ==0)
prots_keep <- both_CL_fin %>% add_count(Cond, Var1) %>% filter(n>1) %>% distinct(Cond, Var1,.keep_all=T)#atleast 2 observations of protein per condition
prots_keep <- prots_keep %>% add_count(Var1) %>% filter(nn==3) %>% distinct(Var1)#make sure the 3 conditions are represented
both_CL_fin <- both_CL[both_CL$Var1%in%prots_keep$Var1,]
unique_prots <- as.vector(unique(prots_keep$Var1))
#differential abundance testing with ANOVA across conditions
prots_keep$pval <- 1000
ncol(prots_keep)

for(i in 1:length(unique_prots)){
  temp_interest <- as.character(unique_prots[i])  #for  a given protein
  prot_int <- both_CL_fin[both_CL_fin$Var1==temp_interest,]
  aav_test <- aov(norm~Cond, data = prot_int)
  sum_aav1 <- summary(aav_test)
  sum_aav1 <- data.frame(sum_aav1[[1]])
  pval_temp <- paste0(sum_aav1[1,5]) 
  prots_keep[i,2] <- as.numeric(pval_temp)
}

prots_keep$pval_adj <-  p.adjust(prots_keep$pval, method = 'BH') # qvalue
prots_sig <- prots_keep[which(prots_keep$pval_adj<0.05),]

prot_info <- both_CL_fin %>% dplyr::select("Var1","Genes") %>% distinct(Var1,.keep_all=T)
both_CL_fin2 <- both_CL_fin %>% left_join(prots_keep, by =c("Var1"="Var1"))
both_CL_fin2 <- both_CL_fin2 %>% dplyr::group_by(Cond, Var1) %>%
  dplyr::mutate("log2_mean_value" = mean(norm, na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::rename("log2_value"="norm")
#dplyr::distinct(Cond, Var1,.keep_all=T)

both_CL_fin3 <- both_CL_fin2 %>% dplyr::rename("Uniprot_PG"="Var1") %>%
  dplyr::select("Uniprot_PG", "Genes", "Cond","CL","log2_value","log2_mean_value","pval_adj")

write.table(both_CL_fin3, "ANOVA_CL7506_CL7496_combined_protein_Abudances_05192023.txt", sep = "\t", row.names = FALSE)

############# Control vs DDX3X KD
###### DA prots
both_CL_fin <- both_CL[!is.na(both_CL$norm),]
both_CL_fin <- both_CL_fin[grepl("Control|DDX3X KD",both_CL_fin$Cond),]
both_CL_fin <- both_CL_fin[both_CL_fin$norm!=0,]
prots_keep <- both_CL_fin %>% add_count(Cond, Var1) %>% filter(n>1) %>% distinct(Cond, Var1,.keep_all=T)#atleast 2 per condition
prots_keep <- prots_keep %>% add_count(Var1) %>% filter(nn==2) %>% distinct(Var1)
both_CL_fin <- both_CL_fin[both_CL_fin$Var1%in%prots_keep$Var1,]
both_CL_fin <- both_CL_fin %>% dplyr::group_by(Var1) %>% dplyr::mutate("norm" = norm-mean(norm,na.rm=T)) %>% ungroup()
#ggplot(both_CL_fin, aes(x=norm,fill=Cond)) + geom_density(alpha=0.5)
unique_prots <- as.vector(unique(prots_keep$Var1))
#differential abundance testing with ANOVA across conditions
prots_keep$pval <- 1000
ncol(prots_keep)

for(i in 1:length(unique_prots)){
  temp_interest <- as.character(unique_prots[i])  #for  a given protein
  prot_int <- both_CL_fin[both_CL_fin$Var1==temp_interest,]
  aav_test <- aov(norm~Cond, data = prot_int)
  sum_aav1 <- summary(aav_test)
  sum_aav1 <- data.frame(sum_aav1[[1]])
  pval_temp <- paste0(sum_aav1[1,5]) 
  prots_keep[i,2] <- as.numeric(pval_temp)
}

prots_keep$pval_adj <-  p.adjust(prots_keep$pval, method = 'BH') # qvalue
prots_sig <- prots_keep[which(prots_keep$pval_adj<0.05),]

prot_info <- both_CL_fin %>% dplyr::select("Var1","Genes") %>% distinct(Var1,.keep_all=T)
both_CL_fin2 <- both_CL_fin %>% left_join(prots_keep, by =c("Var1"="Var1"))
both_CL_fin2 <- both_CL_fin2 %>% dplyr::group_by(Cond, Var1) %>%
  dplyr::mutate("log2_mean_value" = mean(norm, na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::rename("log2_value"="norm")
  #dplyr::distinct(Cond, Var1,.keep_all=T)

both_CL_fin3 <- both_CL_fin2 %>% dplyr::rename("Uniprot_PG"="Var1") %>%
  dplyr::select("Uniprot_PG", "Genes", "Cond","CL","log2_value","log2_mean_value","pval_adj")

write.table(both_CL_fin3, "Control_vs_DDX3X_KD__CL7506_CL7496_combined_protein_Abudances_05192023.txt", sep = "\t", row.names = FALSE)


############# Control vs DDX3Y KD
###### DA prots
ncol(both_CL)
both_CL_fin <- both_CL[!is.na(both_CL$norm),]
both_CL_fin <- both_CL_fin[grepl("Control|DDX3Y KD",both_CL_fin$Cond),]
both_CL_fin <- both_CL_fin[both_CL_fin$norm!=0,]
prots_keep <- both_CL_fin %>% add_count(Cond, Var1) %>% filter(n>1) %>% distinct(Cond, Var1,.keep_all=T)#atleast 2 per condition
prots_keep <- prots_keep %>% add_count(Var1) %>% filter(nn==2) %>% distinct(Var1)
both_CL_fin <- both_CL_fin[both_CL_fin$Var1%in%prots_keep$Var1,]
both_CL_fin <- both_CL_fin %>% dplyr::group_by(Var1) %>% 
  dplyr::mutate("norm" = norm-mean(norm,na.rm=T)) %>% ungroup()
#ggplot(both_CL_fin, aes(x=norm,fill=Cond)) + geom_density(alpha=0.5)
unique_prots <- as.vector(unique(prots_keep$Var1))
#differential abundance testing with ANOVA across conditions
prots_keep$pval <- 1000
ncol(prots_keep)

for(i in 1:length(unique_prots)){
  temp_interest <- as.character(unique_prots[i])  #for  a given protein
  prot_int <- both_CL_fin[both_CL_fin$Var1==temp_interest,]
  aav_test <- aov(norm~Cond, data = prot_int)
  sum_aav1 <- summary(aav_test)
  sum_aav1 <- data.frame(sum_aav1[[1]])
  pval_temp <- paste0(sum_aav1[1,5]) 
  temp_interest <- temp_interest
  prots_keep[i,2] <- as.numeric(pval_temp)
}

prots_keep$pval_adj <-  p.adjust(prots_keep$pval, method = 'BH') # qvalue
prots_sig <- prots_keep[which(prots_keep$pval_adj<0.05),]

prot_info <- both_CL_fin %>% dplyr::select("Var1","Genes") %>% distinct(Var1,.keep_all=T)
both_CL_fin2 <- both_CL_fin %>% left_join(prots_keep, by =c("Var1"="Var1"))
both_CL_fin2 <- both_CL_fin2 %>% dplyr::group_by(Cond, Var1) %>%
  dplyr::mutate("log2_mean_value" = mean(norm, na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::rename("log2_value"="norm")
#dplyr::distinct(Cond, Var1,.keep_all=T)

both_CL_fin3 <- both_CL_fin2 %>% dplyr::rename("Uniprot_PG"="Var1") %>%
  dplyr::select("Uniprot_PG", "Genes", "Cond","CL","log2_value","log2_mean_value","pval_adj")

write.table(both_CL_fin3, "Control_vs_DDX3Y_KD__CL7506_CL7496_combined_protein_Abudances_05192023.txt", sep = "\t", row.names = FALSE)


############# Control vs DDX3X KD
###### DA prots
ncol(both_CL)
both_CL_fin <- both_CL[!is.na(both_CL$norm),]
both_CL_fin <- both_CL_fin[grepl("DDX3X KD|DDX3Y KD",both_CL_fin$Cond),]
both_CL_fin <- both_CL_fin[both_CL_fin$norm!=0,]
prots_keep <- both_CL_fin %>% add_count(Cond, Var1) %>% filter(n>1) %>% distinct(Cond, Var1,.keep_all=T)#atleast 2 per condition
prots_keep <- prots_keep %>% add_count(Var1) %>% filter(nn==2) %>% distinct(Var1)
both_CL_fin <- both_CL_fin[both_CL_fin$Var1%in%prots_keep$Var1,]
both_CL_fin <- both_CL_fin %>% dplyr::group_by(Var1) %>% 
  dplyr::mutate("norm" = norm-mean(norm,na.rm=T)) %>% ungroup()
#ggplot(both_CL_fin, aes(x=norm,fill=Cond)) + geom_density(alpha=0.5)
unique_prots <- as.vector(unique(prots_keep$Var1))
#differential abundance testing with ANOVA across conditions
prots_keep$pval <- 1000
ncol(prots_keep)

for(i in 1:length(unique_prots)){
  temp_interest <- as.character(unique_prots[i])  #for  a given protein
  prot_int <- both_CL_fin[both_CL_fin$Var1==temp_interest,]
  aav_test <- aov(norm~Cond, data = prot_int)
  sum_aav1 <- summary(aav_test)
  sum_aav1 <- data.frame(sum_aav1[[1]])
  pval_temp <- paste0(sum_aav1[1,5]) 
  temp_interest <- temp_interest
  prots_keep[i,2] <- as.numeric(pval_temp)
}

prots_keep$pval_adj <-  p.adjust(prots_keep$pval, method = 'BH') # qvalue
prots_sig <- prots_keep[which(prots_keep$pval_adj<0.05),]

prot_info <- both_CL_fin %>% dplyr::select("Var1","Genes") %>% distinct(Var1,.keep_all=T)
both_CL_fin2 <- both_CL_fin %>% left_join(prots_keep, by =c("Var1"="Var1"))
both_CL_fin2 <- both_CL_fin2 %>% dplyr::group_by(Cond, Var1) %>%
  dplyr::mutate("log2_mean_value" = mean(norm, na.rm=T)) %>% dplyr::ungroup() %>%
  dplyr::rename("log2_value"="norm")
#dplyr::distinct(Cond, Var1,.keep_all=T)

both_CL_fin3 <- both_CL_fin2 %>% dplyr::rename("Uniprot_PG"="Var1") %>%
  dplyr::select("Uniprot_PG", "Genes", "Cond","CL","log2_value","log2_mean_value","pval_adj")

write.table(both_CL_fin3, "DDX3X_KD_vs_DDX3Y_KD__CL7506_CL7496_combined_protein_Abudances_05192023.txt", sep = "\t", row.names = FALSE)

