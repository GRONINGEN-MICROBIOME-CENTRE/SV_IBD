setwd("/groups/umcg-tifn/tmp01/users/umcg-dwang/project/IBD_SV")

info    <- read.table("01.cleanData/Info/Informative_species_information.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)

load("01.cleanData/All/all_phen.RData")
load("01.cleanData/All/all_fmetab.RData")
load("01.cleanData/All/all_abun.RData")
load("01.cleanData/All/distMat/all_msv_dist_std.RData")

load("01.cleanData/LLD/lld_phen.RData")
load("01.cleanData/LLD/lld_fmetab.RData")
load("01.cleanData/LLD/lld_abun.RData")
load("01.cleanData/LLD/distMat/lld_msv_dist_std.RData")

load("01.cleanData/IBD/ibd_phen.RData")
load("01.cleanData/IBD/ibd_fmetab.RData")
load("01.cleanData/IBD/ibd_abun.RData")
load("01.cleanData/IBD/distMat/ibd_msv_dist_std.RData")

if (!dir.exists("06.Species_level_assoc")) {dir.create("06.Species_level_assoc")}
if (!dir.exists("06.Species_level_assoc/RData")) {dir.create("06.Species_level_assoc/RData")}


## Prepare covariate table
# all
all_abun_clr<-abundances(x=as.data.frame(na.omit(all_abun)), transform="clr") %>% as.data.frame
all_abun_clr <- all_abun_clr[match(rownames(all_abun), rownames(all_abun_clr)),]
rownames(all_abun_clr) <- rownames(all_abun)
all_covar<-cbind(all_phen,all_abun_clr)

# lld
lld_abun_clr<-abundances(x=as.data.frame(na.omit(lld_abun)), transform="clr") %>% as.data.frame
lld_abun_clr <- lld_abun_clr[match(rownames(lld_abun), rownames(lld_abun_clr)),]
rownames(lld_abun_clr) <- rownames(lld_abun)
lld_covar<-cbind(lld_phen,lld_abun_clr)

# ibd
ibd_abun_clr<-abundances(x=as.data.frame(na.omit(ibd_abun)), transform="clr") %>% as.data.frame
ibd_abun_clr <- ibd_abun_clr[match(rownames(ibd_abun), rownames(ibd_abun_clr)),]
rownames(ibd_abun_clr) <- rownames(ibd_abun)
ibd_covar<-cbind(ibd_phen,ibd_abun_clr)

covar1 <- c("host_Sex", "host_Age", "host_BMI", "Read_count")
covar2 <- c("host_Sex", "host_Age", "host_BMI", "Read_count","IBD")

##all
all_adonis_res <- my_adonis_terms_adjAbun(all_msv_dist_std[c(1:2)], all_fmetab[,c(1:3)], all_covar, covar2, info)
save(all_adonis_res, file = "06.Species_level_assoc/RData/all_adonis_res.RData")

lld_adonis_res <- my_adonis_terms_adjAbun(lld_msv_dist_std[c(1:2)], lld_fmetab[,c(1:3)], lld_covar, covar1, info)
save(lld_adonis_res, file = "06.Species_level_assoc/RData/lld_adonis_res.RData")

ibd_adonis_res <- my_adonis_terms_adjAbun(ibd_msv_dist_std[c(1:2)], ibd_fmetab[,c(1:3)], ibd_covar, covar1, info)
save(ibd_adonis_res, file = "06.Species_level_assoc/RData/ibd_adonis_res.RData")


## meta-analysis
#load("06.Species_level_assoc/RData/all_adonis_res.RData")
#load("06.Species_level_assoc/RData/lld_adonis_res.RData")
#load("06.Species_level_assoc/RData/ibd_adonis_res.RData")

all_adonis_res.table <- all_adonis_res$table
lld_adonis_res.table <- lld_adonis_res$table
ibd_adonis_res.table <- ibd_adonis_res$table

cbind_adonis_res.table<-cbind(all_adonis_res.table,lld_adonis_res.table,ibd_adonis_res.table)[,-c(8,9,15,16)]
colnames(cbind_adonis_res.table)<-c("Species","BA",
                                    paste("All.",colnames(cbind_adonis_res.table)[3:7],sep = ""),
                                    paste("LLD.",colnames(cbind_adonis_res.table)[3:7],sep = ""),
                                    paste("IBD",colnames(cbind_adonis_res.table)[3:7],sep = ""))
adonis_res <- my_batch_meta_p(cbind_adonis_res.table, c("LLD", "300OB"), c(10,15), c(9,14))
save(adonis_res, file = "06.Species_level_assoc/RData/adonis_res.RData")