## libraty packages
library(tidyverse)
library(ggpubr)
library(GGally)
library(ggmosaic)
library(cowplot)
library(ggsci)
library(ggthemes)
library(gplots)
library(VennDiagram)
library(circlize)
library(viridis)
library(wesanderson)
library(ggalluvial)
library(factoextra)

library(reshape2)
library(Hmisc)
library(vegan)
library(microbiome)
library(mediation)
library(karyoploteR)

mycolor2_blue_red = c("#28B6F6FF","#EE5250FF")
mycolor2_green_blue  <- c("#2EC4B6","#235789")
mycolor3 <- c("#2EC4B6","#235789", "grey50") # c("#4f94cd","#ff4040", "grey50")

## Convert SV name
changeSVname<-function(SVrawid){
  testname     <- SVrawid
  species_name <- as.character(taxonomy$organism[match(str_replace_all(testname, "\\..*",""), taxonomy$X)])
  region       <- str_replace_all(testname, ".*\\:","") 
  region_list  <- str_split(region,";") 
  
  region_suf   <- NULL
  i<-1
  for (i in c(1:length(region_list))){
    if(length(region_list[[i]]) <= 2){
      region_suf[i] <- paste(":",region[i],sep = "")
    }else{
      region_suf[i] <- paste(":",region_list[[i]][1]," and ",length(region_list[[i]])-1," segments", sep = "")
    }
    i <- i+1
  }
  paste(species_name,region_suf,sep = "")
}



## Calculate SV size
calcSVSize <- function(SVrawid){
  testname <- SVrawid
  region   <- str_replace_all(testname, ".*\\:","") 
  region_list <- str_split(region,";") 
  
  sv_size  <- NULL
  for (i in c(1:length(region_list))) {
    frag_df    <- str_split(region_list[[i]], "_") %>% unlist %>% as.numeric  %>% matrix(byrow = T, ncol = 2) %>% as.data.frame 
    sv_size[i] <- sum(as.numeric(frag_df$V2)-as.numeric(frag_df$V1))
  }
  sv_size
}



## pie chart
my_pie<-function(in1,item,mycol = c("#F24D4D","#4D7EAE")){
  require(tibble)
  require(showtext)
  require(Cairo)
  require(ggsci)
  require(tibble)
  require(scales)
  require(ggrepel)
  require(forcats)
  require(scatterpie)
  
  showtext_auto()
  
  in1.tibble            <- as_tibble(subset(in1, in1$items%in%item))
  in1.tibble$categories <- fct_reorder(in1.tibble$categories, in1.tibble$value)
  in1.tibble            <- in1.tibble[order(in1.tibble$value, decreasing = TRUE), ]
  piepercent            <- round(100*in1.tibble$value/sum(in1.tibble$value), 2)
  my_labels             <- tibble(x.breaks = seq(1.3, 1.3, length.out = length(piepercent)),
                                  y.breaks = cumsum(in1.tibble$value) - in1.tibble$value/2,
                                  labels = paste(in1.tibble$categories, "\n",in1.tibble$value,", ",piepercent, "%", sep = ""),
                                  categories = in1.tibble$categories)
  
  pdfName    <- paste(item, ".pie.pdf", sep = "")
  ggplot(in1.tibble, aes(x = 1, y = value, fill = categories)) +
    ggtitle(paste(item)) +
    geom_bar(stat="identity", color='white') + 
    coord_polar(theta='y') + 
    theme(legend.position = "None",
          axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    #    scale_fill_brewer(palette = "Set3", direction = -1)+
    scale_fill_manual(values=mycol) + # set fill color manually
    geom_text_repel(data = my_labels, 
                    aes(x = x.breaks, y = y.breaks, label = labels),
                    size = 7,        # label text size
                    show.legend = FALSE,
                    inherit.aes = FALSE,
                    arrow = arrow(length=unit(0.01, "npc")),
                    force = 1,
                    segment.color = 'grey50'
    )
}



## Calculate SE
se <- function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}


## Summary statistics
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N            = length2(xx[[col]], na.rm=na.rm),
                     nonZero_N    = sum(xx[[col]]!=0,na.rm = na.rm),
                     nonZero_rate = sum(xx[[col]]!=0,na.rm = na.rm)/length2(xx[[col]], na.rm=na.rm),
                     mean         = mean   (xx[[col]], na.rm=na.rm),
                     sd           = sd     (xx[[col]], na.rm=na.rm),
                     mindata      = min(xx[[col]], na.rm=na.rm),
                     maxdata      = max(xx[[col]], na.rm=na.rm),
                     quant1_  = quantile(xx[[col]], na.rm=na.rm)[2] ,
                     quant2_  = quantile(xx[[col]], na.rm=na.rm)[4] ,
                     quantmid = quantile(xx[[col]], na.rm=na.rm)[3] ,
                     mea      = mean   (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  
  ciMult   <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  datac$xd <- datac$mea-datac$se
  datac$xu<- datac$mea+datac$se
  return(datac)
}

shared_sv_dis<-function(inDf){
  #inDf<-all_sv
  inList <- split(as.matrix(inDf), seq(nrow(inDf)))
  
  shared_n_func<-function(inVec,Vec_i){
    #inVec<-inList[[3]]
    Vec_i<-Vec_i
    insvdf<-data.frame(inVec,Vec_i) %>% na.omit
    sv_dis<- vegdist(t(insvdf), method = "canberra")
    #length(na.omit(Vec_i+inVec))
    return(sv_dis)
  }
  
  marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))
  for (i in 1:nrow(inDf)) { #nrow(inDf)
    #i<-2
    Vec_i<-inList[[i]]
    shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)
    marker_n_mat[i,]<-shared_n_i
  }
  rownames(marker_n_mat)<-rownames(inDf)
  colnames(marker_n_mat)<-rownames(inDf)
  
  return(marker_n_mat)
}


## Remove samples with NA in distance matrix
dist_rmna<-function(inDist){
  while(sort(colSums(is.na(inDist)),decreasing = T)[1] > 0){
    rmid<-names(sort(colSums(is.na(inDist)),decreasing = T)[1])
    inDist<-inDist[-match(rmid, rownames(inDist)),-match(rmid, colnames(inDist))]
  }
  return(inDist)
}


get_PCs<-function(inDist, eig.cutoff = 0.6, pc.cutoff = 10){
  #inDist <- all_msv_dist_std[[1]]
  #eig.cutoff <- 0.6
  #eig.cutoff<-as.numeric(eig.cutoff)
  
  pcoa_res <- cmdscale(inDist, k=30, eig = T)
  
  total_eig <- 0
  i <- 0
  
  while (total_eig<eig.cutoff & i< pc.cutoff) {
    i<-i+1
    total_eig <- total_eig + pcoa_res$eig[i]/100
  }
  
  if(i<=1){
    pcoa<-data.frame(X1 = pcoa_res$points[, 1])
  }else{
    pcoa <- data.frame(pcoa_res$points)[,c(1:i)]
  }
  
  return(list(PCoA= pcoa,PC_num = i,Total_eig = total_eig))
  
}



## Batch adonis
my_adonis_terms<-function(inDist, inDf, covDf=NULL, covar=NULL,n_parallel = 4,n_permutations = 999){
  require(vegan)
  ## test
  #inDist <- distList[[1]]
  #inDf   <- lld_ba
  #covDf  <- ibd_covar
  #covar  <- inCovar
  ## test
  
  inDf   <- as.data.frame(inDf)
  rownames(inDf) <- rownames(covDf)
  covDf<- covDf[,covar]

  inDist_rmna <- dist_rmna(inDist)
  covDf_rmna   <- na.omit(covDf)
  
  my_adonis<-function(inMat, inVec, covdf = NULL){
    ## test
    #inMat<-inDist_rmna
    #inVec<-inDf[,976]
    #covdf<-covDf_rmna
    ## test
    inMat_covdf_inter <- intersect(rownames(inMat),rownames(covdf))
    inMat_covdf_inVec_inter <- intersect(rownames(inDf)[!is.na(inVec)], inMat_covdf_inter)
    
    inMat_rmna <- inMat[match(inMat_covdf_inVec_inter,rownames(inMat)),
                        match(inMat_covdf_inVec_inter,colnames(inMat))]
    inVec_rmna <- inVec[match(inMat_covdf_inVec_inter, rownames(inDf))]
    covdf_rmna <- covdf[match(inMat_covdf_inVec_inter, rownames(covdf)),]
    
    in_cov <- data.frame(covdf_rmna,inVec_rmna)
    
    sample_size<-nrow(in_cov)
    
    if(sample_size>10){
      adonis_res <- adonis2(as.dist(inMat_rmna)~.,data = in_cov,parallel = n_parallel,permutations = n_permutations)
      return(list(adonis_res$aov.tab[ncol(in_cov),5],adonis_res$aov.tab[ncol(in_cov),6], sample_size))
    }else{
      return(list(NA,NA, sample_size))
    }
  }
  
  adonis_res_batch <- apply(inDf, 2, my_adonis,inMat = inDist_rmna, covdf = covDf_rmna)
  adonis_res_batch.unlist <- matrix(unlist(adonis_res_batch), ncol = 3, byrow = T)
  colnames(adonis_res_batch.unlist)<-c("R2", "P", "N")
  
  return(adonis_res_batch.unlist)
}


my_adonis_terms_adjAbun<-function(distList, inDf, covDf, covar,info,n_parallel = 4,n_permutations = 999){
  #distList<-ibd_msv_dist_std[c(20)]
  #inDf<-ibd_fmetab
  #covDf<-ibd_covar
  #covar<-covar1
  #info<-info

  res_table<-NULL
  for (i in 1:length(distList)) { # 
    #i<-1
    cat(paste(i,"\n")) 
    
    n<-names(distList)[i] %>% 
      str_replace_all("msv_","") %>% 
      match(., info$organism)
    
    inCovar<-c(covar, info$organism[n])
    
    
    adonis_res<-my_adonis_terms(distList[[i]], inDf, covDf = covDf, covar = inCovar,n_parallel = n_parallel,n_permutations = n_permutations)
    
    adonis_res_df<-data.frame(Species = rep(info$Short_name[n], dim(adonis_res)[1]),
                              BA = colnames(inDf),
                              as.data.frame(adonis_res))
    
    res_table<-rbind(res_table, adonis_res_df)
  }
  
  res_table$Species<-as.character(res_table$Species)
  res_table$BA<-as.character(res_table$BA)
  
  res_table$fdr<-p.adjust(res_table$P,method = 'fdr')
  res_table$bonferroni<-p.adjust(res_table$P,method = 'bonferroni')
  
  # R2
  r2_mat<-matrix(data  = res_table$R2,
                 nrow  = length(unique(res_table$Species)),
                 ncol  = length(unique(res_table$BA)),
                 byrow = T)
  
  rownames(r2_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(r2_mat)<-colnames(inDf)
  
  # P
  p_mat<-matrix(data  = res_table$P,
                nrow  = length(unique(res_table$Species)),
                ncol  = length(unique(res_table$BA)),
                byrow = T)
  
  rownames(p_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(p_mat)<-colnames(inDf)
  
  # fdr
  fdr_mat<-matrix(data  = res_table$fdr,
                  nrow  = length(unique(res_table$Species)),
                  ncol  = length(unique(res_table$BA)),
                  byrow = T)
  
  rownames(fdr_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(fdr_mat)<-colnames(inDf)
  
  return(list(table = res_table,
              r2 = r2_mat,
              P = p_mat,
              FDR = fdr_mat))
}

my_meta_p <- function(inVec, study_name, n_col, p_col) {
  require(metap)
  #require(metafor)
  
  #inVec<-cbind_adonis_res.table[20,]
  #study_name<-c("LLD", "300OB")
  #n_col<-c(10,15)
  #_col<-c(9,14)
  
  study_n <- inVec[n_col] %>% as.numeric
  study_p <- inVec[p_col] %>% as.numeric
  
  inDf <- data.frame(study_name, study_n, study_p)
  
  Wi=sqrt(study_n)
  #Convert p-values to Z-scores
  Zi=qnorm(1-(study_p/2))
  # Meta-zscore  
  Z=(sum(Zi*Wi)/sqrt(sum(Wi^2)))
  # Convert Z-score to p-value
  MetaP=2*pnorm(-abs(Z))
  
  #Cochran Q-test (test heterogeneity of your meta-analysis)
  
  WeightedZ= sum(sqrt(study_n)*Zi)
  totalSample=sum(study_n)
  
  #Calculate expected Z
  expected_Z= sqrt(study_n)*WeightedZ/totalSample
  het_Sum= sum((Zi-expected_Z)*(Zi-expected_Z))
  
  #Get p-value of heterogeneity test!
  my_pvalue_co=pchisq(het_Sum, lower.tail = F, df=length(study_p)-1)
  
  return(
    list(
      Meta.p = MetaP,
      Meta.hetero = het_Sum,
      Meta.hetero.p = my_pvalue_co
    )
  )
}
## Batch meta-analysis
my_batch_meta_p <- function(inDf,study_name,n_col,p_col,row_var_col = 1,col_var_col = 2) {
  #inDf<-cbind_adonis_res.table[c(1:10),]
  #study_name<-c("LLD", "300OB")
  #n_col<-c(10,15)
  #p_col<-c(9,14)
  #row_var_col<-1
  #col_var_col<-2
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_p,
    study_name = study_name,
    n_col = n_col,
    p_col = p_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 3, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.p',  "Meta.hetero", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  rowName <- unique(inDf[, row_var_col])
  colName <- unique(inDf[, col_var_col])
  
  N_row <- length(rowName)
  N_col <- length(colName)
  
  
  batch_res.p <- matrix(batch_res.unlist[, 1], ncol = N_col, byrow = F)
  colnames(batch_res.p) <- colName
  rownames(batch_res.p) <- rowName
  batch_res.p[is.na(batch_res.p)] <- 0
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  # fdr matrix
  batch_res.fdr <- matrix(
    data  = batch_res_edge$Meta.fdr.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  colnames(batch_res.fdr) <- colName
  rownames(batch_res.fdr) <- rowName
  batch_res.fdr[is.na(batch_res.fdr)] <- 1
  
  # bonferroni matrix
  batch_res.bon <- matrix(
    data  = batch_res_edge$Meta.bonferroni.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  
  colnames(batch_res.bon) <- colName
  rownames(batch_res.bon) <- rowName
  batch_res.bon[is.na(batch_res.bon)] <- 1
  
  return(
    list(
      table      = batch_res_edge,
      p          = batch_res.p,
      fdr        = batch_res.fdr,
      bonferroni = batch_res.bon
    )
  )
}


## The Normal Quantile Transformation
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}


## linear model
lm_btw_mats<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #y_mat<- y_mat_i #lld_intri[,c(1:3)]
  #x_mat<- x_mat #lld_vsv[,c(1:3)]
  #cov_mat<- cov_mat# lld_covar
  #covar<- covar_i #covar
  ## test block
  
  my_lm<-function(y,x){
    # y<-y_mat[,1]
    #  x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

lm_btw_mats_v2<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #  y_mat<- all_abun_clr[,info$organism] #lld_intri[,c(1:3)]
  # x_mat<- as.data.frame(all_phen[,c("IBD")]) #lld_vsv[,c(1:3)]
  # cov_mat<- all_covar# lld_covar
  #  covar<- covar1 #covar
  ## test block
  
  my_lm<-function(y,x){
    #y<-y_mat[,1]
    #x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    try(anv_res <- anova(lm(Y~.,data = lm_input)), silent = T)
    r_sq<-anv_res$`Sum Sq`/sum(anv_res$`Sum Sq`)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    try(R_sq    <- r_sq[match(indv, rownames(anv_res))])
    try(anv_P   <- anv_res[match(indv, rownames(anv_res)), 5])
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    R_sq = R_sq,
                    anv_P = anv_P,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 12, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"),
                       anv_fdr.p = p.adjust(y_x.unlist[,5], method = "fdr"),
                       anv_bonferroni.p = p.adjust(y_x.unlist[,5], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","R_sq","anv_p",
                        "N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate",
                        "fdr.p","bonferroni.p", "anv_fdr.p", "anv_bonferroni.p")
  
  return(y_x_edge)
}





## Meta-analysis
my_meta_lm <- function(inVec, study_name, beta_col, se_col) {
  require(meta)
  require(metafor)
  
  ## test start
  #inVec<-vsv_crass_lm_res.edge[100,]
  #study_name<-c("LLD", "300OB", "IBD")
  #beta_col<-c(3,10,17)
  #se_col<-c(4,11,18)
  
  ## test end
  
  study_beta <- inVec[beta_col] %>% as.numeric
  study_se <- inVec[se_col] %>% as.numeric
  
  #study_beta<-c(0.7, 0.65)
  #study_se<-c(0.07, 0.08)
  #stydy_n<-c(1000, 300)
  
  
  inDf <- data.frame(study_name, study_beta, study_se)
  
  m.hksj <- metagen(
    study_beta,
    study_se,
    data = inDf,
    studlab = study_name
  )
  
  m.hksj.res <- summary(m.hksj)
  
  return(
    list(
      Meta.beta = m.hksj.res$random$TE,
      Meta.se = m.hksj.res$random$seTE,
      Meta.p = m.hksj.res$random$p,
      Meta.I2 = m.hksj.res$I2$TE,
      Meta.hetero.p = m.hksj$pval.Q
    )
  )
}



## Batch meta-analysis
my_batch_meta_lm <- function(inDf,study_name,beta_col,se_col) {
  ## test start
  #inDf<-vsv_crass_lm_res.edge[c(1:10),]
  #study_name<-c("LLD", "300OB","IBD")
  #beta_col<-c(3,10,17)
  #se_col<-c(4,11,18)
  ## test end
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_lm,
    study_name = study_name,
    beta_col = beta_col,
    se_col = se_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 5, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.beta', 'Meta.SE', "Meta.p", "Meta.I2", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  return(batch_res_edge)
}


lm_btw_mats_adjAbun<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col){
  #y_mat<-ob_contin[,c(1:10)]
  #x_mat<-ob_vsv[,c(1:100)]
  #cov_mat<-ob_basic
  #covar<-covar
  #abun<-ob_abun_clr
  #info<-info
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lm_btw_mats_adjAbun2<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col){
# y_mat<-all_vsv.sig
# x_mat<-all_exp.sig
# cov_mat<-all_covar
# covar<-covar
# abun<-all_abun_clr
# info<-info
# abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
#    i<-7
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]
    
    if(!is.null(dim(y_mat_i)[2])){
      if(dim(y_mat_i)[2]>0){
        if(info[i,abun_col]%in%colnames(abun)){
          covar_i<-c(covar,info[i,abun_col])
          y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
          y_x_i$Adjust_Abundance<-"Yes"
          y_x.edge<-rbind(y_x.edge,y_x_i)
        }else{
          covar_i<-covar
          y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
          y_x_i$Adjust_Abundance<-"No"
          y_x.edge<-rbind(y_x.edge,y_x_i)
        }
      }
    }
    
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}


lm_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat, info, abun_col){
  #  y_mat<-lld_ba[,c(1:10)]
  #  x_mat<-vsgv_lld[,c(1:100)]
  #  cov_mat<-lld_basic
  #  covar<-covar
  #  abun<-lld_abun_clr
  #  pc_mat<-lld_msv_pc_cum0.6
  #  info<-info
  
  cov_mat<-cbind(cov_mat, abun,pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #    i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))] 
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}



## logistic model
lr_btw_mats<-function(y_mat,x_mat,cov_mat,covar, direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  
  require(reshape2)
  require(glmnet)
  require(R.utils)
  
  ## test block
  y_mat <- y_mat
  #x_mat <- x_mat_i
  #cov_mat <- cov_mat
  #covar   <- covar #covar
  #direction<-c(1,1)
  ## test block
  
  my_lm<-function(y,x){
    #y<-y_mat[,2]
    #x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    y_0_N<-NA
    y_1_N<-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    y_0_N<-sum(lm_input[,"Y"]==0)
    y_1_N<-sum(lm_input[,"Y"]==1)
    
    lm_input_tmp <- apply(lm_input[,c(2:ncol(lm_input))], 2, qtrans) %>% as.data.frame
    lm_input<-data.frame(Y = lm_input[,1],lm_input_tmp)
    
    try(lm_res <- summary(glm(Y~., data = lm_input, family = 'binomial')))
    
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate,
                    y_0_N = y_0_N,
                    y_1_N = y_1_N)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 12, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(as.data.frame(y_mat)), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2:1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","y_0_N","y_1_N","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}


lr_btw_mats_adjAbun<-function(y_mat, x_mat, cov_mat, covar, abun,info, abun_col){
#  y_mat<-all_phen[,c("IBD")]
#  x_mat<-all_vsv
#  cov_mat<-all_covar
#  covar<-covar
#  abun<-all_abun_clr
#  info<-info
#  abun_col = 9
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
 #   i<-1
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}


lr_btw_mats_adjAbun2<-function(y_mat, x_mat, cov_mat, covar, abun,info, abun_col){
  #y_mat<-lld_intri[,c(1:10)]
  #x_mat<-lld_vsv[,c(1:100)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun
  #info<-info
  #direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]
    
    if(!is.null(dim(y_mat_i)[2])){
      if(dim(y_mat_i)[2]>0){
        if(info[i,abun_col]%in%colnames(abun)){
          covar_i<-c(covar,info[i,abun_col])
          y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
          y_x.edge<-rbind(y_x.edge,y_x_i)
        }else{
          covar_i<-covar
          y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
          y_x.edge<-rbind(y_x.edge,y_x_i)
        }
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}





## Clustering analysis
my_cluster<-function(inDist, ps.cutoff=0.8,my_seed = 1){
  require(NbClust)
  require(fpc)
  require(tsne)
  require(ggsci)
  
  # test
  #inDist<-all_msv_dist_std[[5]]
  # test
  
  clu_n<-prediction.strength(as.dist(inDist), Gmin=2, Gmax=10, M=50,
                             clustermethod=claraCBI ,usepam=T,diss = T,
                             classification="centroid",cutoff=ps.cutoff,
                             distances=T,count=F)
  clu<-claraCBI(as.dist(inDist), clu_n$optimalk, usepam=T,diss=T)
  clu_df<-as.data.frame(clu$partition)
  
  rownames(clu_df)<-rownames(inDist)
  
  pcoa_res<-cmdscale(inDist, k=5, eig = T)
  pcoa <- data.frame(pcoa_res$points)
  
  set.seed(my_seed)
  tsne_res<-Rtsne::Rtsne(inDist, is_distance = TRUE,
                         perplexity = as.integer((nrow(inDist)-1)/3),
                         theta = 0, pca = T,
                         eta=as.integer(nrow(inDist)/12))
  
  tsne_df <- tsne_res$Y %>%
    data.frame() %>%
    setNames(c("X", "Y"))
  tsne_df <- data.frame(Cluster = as.factor(clu$partition),tsne_df)
  
  return(list(pcoa=pcoa, pcoa_res = pcoa_res, clu_n=clu_n,tsne_df=tsne_df))
}



## kruskal-wallis test
kw_btw_mats<-function(mat0,mat1,direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  require(reshape2)
  require(rstatix)
  
  #mat0<-all_ba[,c(1:10)]
  #mat1<-all_msv_cluster_sub[,c(1:10)]
  #mat2<-covar
  #direction<-c(1,1)
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  col_var<-mat0
  row_var<-mat1
  
  my_kw<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,1]
    
    effSize  <- NA
    p.value  <- NA
    N        <- NA
    N_minor  <- NA
    uniq_N   <- NA
    
    kw_input <-data.frame(phen = a,
                          mbio = b)
    kw_input<-na.omit(kw_input)
    
    res.kw <- kruskal_test(phen ~ mbio, data = kw_input)
    res.kw.eff<-kruskal_effsize(phen ~ mbio, data = kw_input)
    
    effSize   <- res.kw.eff$effsize
    p.value   <- res.kw$p
    N         <- nrow(kw_input)
    N_minor   <- table(kw_input$mbio)[(length(table(kw_input$mbio))+1-rank(table(kw_input$mbio)))==2]
    uniq_N    <- length(unique(kw_input$phen))
    
    try(return(list(effSize = effSize, p.value = p.value,N=N,N_minor=N_minor,uniq_N=uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_kw(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 5, byrow = T)
  
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,col_var_row_var.unlist)[,-c(3)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "effectSize", "p","N","N_minor","uniq_N")
  
  
  # add p adjust
  col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
                                   fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
                                   bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$fdr.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  # bonferroni matrix
  col_var_row_var.bon<-matrix(data  = col_var_row_var_edge$bonferroni.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  
  colnames(col_var_row_var.bon)<-colnames(col_var)
  rownames(col_var_row_var.bon)<-colnames(row_var)
  col_var_row_var.bon[is.na(col_var_row_var.bon)]<-1
  
  return(list(table      = col_var_row_var_edge,
              effectSize = col_var_row_var.beta,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr,
              bonferroni = col_var_row_var.bon))
  
}

## Permutation kruskal-wallis test
permKW_btw_mats<-function(mat0,mat1,direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  require(reshape2)
  require(rstatix)
  require(coin)
  
  #mat0<-all_ba[,c(1:10)]
  #mat1<-all_msv_cluster_sub[,c(1:10)]
  #direction<-c(1,1)
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  col_var<-mat0
  row_var<-mat1
  
  my_kw<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,3]
    
    effSize  <- NA
    p.value  <- NA
    perm.fdr <- NA
    N        <- NA
    N_minor  <- NA
    uniq_N   <- NA
    
    kw_input <-data.frame(phen = a,
                          mbio = b)
    kw_input<-na.omit(kw_input)
    
    if(length(unique(kw_input$mbio))>=2){
      res.kw <- kruskal_test(phen ~ as.factor(mbio), data = kw_input)
      res.kw.eff<-kruskal_effsize(phen ~ as.factor(mbio), data = kw_input)
      res.kw.perm <- oneway_test(phen ~ as.factor(mbio), data=kw_input)
      
      effSize   <- res.kw.eff$effsize
      p.value   <- pvalue(res.kw)
      perm.fdr  <- pvalue(res.kw.perm)
      N         <- nrow(kw_input)
      N_minor   <- table(kw_input$mbio)[(length(table(kw_input$mbio))+1-rank(table(kw_input$mbio),ties.method = "last"))==2]
      uniq_N    <- length(unique(kw_input$phen))
      
    }
    
    try(return(list(effSize = effSize,
                    p.value = p.value,
                    perm.fdr=perm.fdr , 
                    N=N,
                    N_minor=N_minor,
                    uniq_N=uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_kw(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 6, byrow = T)
  
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,col_var_row_var.unlist)[,-c(3)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "effectSize", "p","perm.fdr","N","N_minor","uniq_N")
  
  
  # add p adjust
  #col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
  #                                 fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
  #                                 bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$perm.fdr,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  
  return(list(table      = col_var_row_var_edge,
              effectSize = col_var_row_var.beta,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr))
  
}



lm_btw_mats_interaction<-function(y_mat,x_mat,cov_mat,covar,int,transf){
  require(reshape2)
  require(R.utils)
  require(broom)
  
  ## test block
  #  y_mat<- all_fmetab[,c(1:10)] #lld_intri[,c(1:3)]
  #  x_mat<- all_vsv[,c(1:10)] #lld_vsv[,c(1:3)]
  #  cov_mat<- all_covar # lld_covar
  #  covar<- covar #covar
  #  int<-"IBD"
  #  transf<-c("Y", "X1", "host_Age", "host_BMI", "Read_count")
  ## test block
  int_vector<-cov_mat[,int]
  
  my_lm_inter<-function(y,x1,x2){
    #y<-y_mat[,1]
    #x1<-x_mat[,1]
    #x2<-int_vector
    
    add.X1.beta <- NA
    add.X1.se <- NA
    add.X1.p<- NA
    add.X2.beta <- NA
    add.X2.se <- NA
    add.X2.p <- NA
    add.R_sq <- NA
    add.anv_P <- NA
    
    int.X1.beta <- NA
    int.X1.se <- NA
    int.X1.p <- NA
    int.X2.beta <- NA
    int.X2.se <- NA
    int.X2.p <- NA
    int.X1xX2.beta <- NA
    int.X1xX2.se <- NA
    int.X1xX2.p <- NA
    int.R_sq <- NA
    int.anv_P <- NA
    
    N       <- NA
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X1 = x1, X2 = x2, cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X1"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X1"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    for (i in 1:length(transf)) {
      lm_input[,match(transf[i], colnames(lm_input))]<-qtrans(lm_input[,match(transf[i], colnames(lm_input))])
    }
    lm_input<-as.data.frame(lm_input)
    #transf_index<-match(transf, colnames(lm_input))
    #lm_input[,transf_index] <- apply(lm_input[,transf_index], 2, qtrans) %>% as.data.frame
    
    # Additive model
    try(lm_add_res <- summary(lm(Y~.,data = lm_input)),silent = T)
    #    anv_add_res <- anova(lm(Y~.,data = lm_input))
    # Interaction model
    try(lm_int_res <- summary(lm(Y~.+X1*X2,data = lm_input)),silent = T)
    #    anv_int_res <- anova(lm(Y~.+X1*X2,data = lm_input))
    
    indv1<-'X1'
    indv2<-'X2'
    int_term<-'X1:X2'
    
    ##
    try(add.X1.beta    <- lm_add_res$coefficients[match(indv1,rownames(lm_add_res$coefficients)),1],silent = T)
    try(add.X1.se      <- lm_add_res$coefficients[match(indv1,rownames(lm_add_res$coefficients)),2], silent = T)
    try(add.X1.p       <- lm_add_res$coefficients[match(indv1,rownames(lm_add_res$coefficients)),4],silent = T)
    try(add.X2.beta    <- lm_add_res$coefficients[match(indv2,rownames(lm_add_res$coefficients)),1],silent = T)
    try(add.X2.se      <- lm_add_res$coefficients[match(indv2,rownames(lm_add_res$coefficients)),2], silent = T)
    try(add.X2.p       <- lm_add_res$coefficients[match(indv2,rownames(lm_add_res$coefficients)),4],silent = T)
    try(add.R_sq       <- lm_add_res$adj.r.squared, silent = T)
    try(add.anv_P      <- glance(lm_add_res)$p.value,silent = T)
    
    try(int.X1.beta    <- lm_int_res$coefficients[match(indv1,rownames(lm_int_res$coefficients)),1],silent = T)
    try(int.X1.se      <- lm_int_res$coefficients[match(indv1,rownames(lm_int_res$coefficients)),2], silent = T)
    try(int.X1.p       <- lm_int_res$coefficients[match(indv1,rownames(lm_int_res$coefficients)),4],silent = T)
    try(int.X2.beta    <- lm_int_res$coefficients[match(indv2,rownames(lm_int_res$coefficients)),1],silent = T)
    try(int.X2.se      <- lm_int_res$coefficients[match(indv2,rownames(lm_int_res$coefficients)),2], silent = T)
    try(int.X2.p       <- lm_int_res$coefficients[match(indv2,rownames(lm_int_res$coefficients)),4],silent = T)
    try(int.X1xX2.beta    <- lm_int_res$coefficients[match(int_term,rownames(lm_int_res$coefficients)),1],silent = T)
    try(int.X1xX2.se      <- lm_int_res$coefficients[match(int_term,rownames(lm_int_res$coefficients)),2], silent = T)
    try(int.X1xX2.p       <- lm_int_res$coefficients[match(int_term,rownames(lm_int_res$coefficients)),4],silent = T)
    try(int.R_sq       <- lm_int_res$adj.r.squared, silent = T)
    try(int.anv_P      <- glance(lm_int_res)$p.value,silent = T)
    
    try(return(list(add.X1.beta = add.X1.beta,
                    add.X1.se   = add.X1.se,
                    add.X1.p    = add.X1.p,
                    add.X2.beta = add.X2.beta,
                    add.X2.se   = add.X2.se,
                    add.X2.p    = add.X2.p,
                    add.R_sq    = add.R_sq,
                    add.anv_P   = add.anv_P,
                    
                    int.X1.beta = int.X1.beta,
                    int.X1.se   = int.X1.se,
                    int.X1.p    = int.X1.p,
                    int.X2.beta = int.X2.beta,
                    int.X2.se   = int.X2.se,
                    int.X2.p    = int.X2.p,
                    int.X1xX2.beta = int.X1xX2.beta,
                    int.X1xX2.se = int.X1xX2.se,
                    int.X1xX2.p = int.X1xX2.p,
                    int.R_sq = int.R_sq,
                    int.anv_P = int.anv_P,
                    
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm_inter(a,b,int_vector),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 26, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       rep(int, nrow(y_x_edge_beta)),
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("Y", "X1","X2",
                        "add.X1.beta","add.X1.se", "add.X1.p","add.X2.beta","add.X2.se", "add.X2.p",
                        "add.R_sq","add.anv_P",
                        "int.X1.beta","int.X1.se", "int.X1.p","int.X2.beta","int.X2.se", "int.X2.p","int.X1xX2.beta","int.X1xX2.se", "int.X1xX2.p",
                        "int.R_sq","int.anv_P",
                        "N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate")
  
  y_x_edge$add.anv_fdr.p<-p.adjust(y_x_edge$add.anv_P, method = "fdr")
  y_x_edge$add.anv_bonferroni.p<-p.adjust(y_x_edge$add.anv_P, method = "bonferroni")
  y_x_edge$int.anv_fdr.p<-p.adjust(y_x_edge$int.anv_P, method = "fdr")
  y_x_edge$int.anv_bonferroni.p<-p.adjust(y_x_edge$int.anv_P, method = "bonferroni")
  y_x_edge$int.X1xX2.fdr.p<-p.adjust(y_x_edge$int.X1xX2.p, method = "fdr")
  y_x_edge$int.X1xX2.bonferroni.p<-p.adjust(y_x_edge$int.X1xX2.p, method = "bonferroni")
  y_x_edge$int_add_R_sq_diff<-y_x_edge$int.R_sq-y_x_edge$add.R_sq
  
  return(y_x_edge)
}



lm_btw_mats_interaction_adjAbun<-function(y_mat, x_mat, cov_mat, covar,int,transf, abun, info, abun_col){
  #  y_mat<- all_fmetab[,c(1:10)] #lld_intri[,c(1:3)]
  #  x_mat<- all_vsv #lld_vsv[,c(1:3)]
  #  cov_mat<- all_covar # lld_covar
  #  covar<- covar #covar
  # int<-"IBD"
  # transf<-c("Y", "X1", "host_Age", "host_BMI", "Read_count")
  #  abun<-all_abun_clr
  #  info<-info
  #  abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #    i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats_interaction(y_mat, x_mat_i, cov_mat, covar_i,int,transf)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i,int,transf)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  y_x.edge$add.anv_fdr.p<-p.adjust(y_x.edge$add.anv_P,method = 'fdr')
  y_x.edge$add.anv_bonferroni.p<-p.adjust(y_x.edge$add.anv_P,method = 'bonferroni')
  y_x.edge$int.anv_fdr.p<-p.adjust(y_x.edge$int.anv_P,method = 'fdr')
  y_x.edge$int.anv_bonferroni.p<-p.adjust(y_x.edge$int.anv_P,method = 'bonferroni')
  y_x.edge$int.X1xX2.fdr.p<-p.adjust(y_x.edge$int.X1xX2.p,method = 'fdr')
  y_x.edge$int.X1xX2.bonferroni.p<-p.adjust(y_x.edge$int.X1xX2.p,method = 'bonferroni')
  
  return(y_x.edge)
}

lm_btw_mats_interaction_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar,int,transf, abun,pc_mat, info, abun_col){
  #  y_mat<- all_fmetab[,c(1:10)] #lld_intri[,c(1:3)]
  #  x_mat<- all_vsv #lld_vsv[,c(1:3)]
  #  cov_mat<- all_covar # lld_covar
  #  covar<- covar #covar
  # int<-"IBD"
  # transf<-c("Y", "X1", "host_Age", "host_BMI", "Read_count")
  #  abun<-all_abun_clr
  #  pc_mat<-lld_msv_pc_cum0.6
  #  info<-info
  #  abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun,pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #    i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats_interaction(y_mat, x_mat_i, cov_mat, covar_i,int,transf)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i,int,transf)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  y_x.edge$add.anv_fdr.p<-p.adjust(y_x.edge$add.anv_P,method = 'fdr')
  y_x.edge$add.anv_bonferroni.p<-p.adjust(y_x.edge$add.anv_P,method = 'bonferroni')
  y_x.edge$int.anv_fdr.p<-p.adjust(y_x.edge$int.anv_P,method = 'fdr')
  y_x.edge$int.anv_bonferroni.p<-p.adjust(y_x.edge$int.anv_P,method = 'bonferroni')
  y_x.edge$int.X1xX2.fdr.p<-p.adjust(y_x.edge$int.X1xX2.p,method = 'fdr')
  y_x.edge$int.X1xX2.bonferroni.p<-p.adjust(y_x.edge$int.X1xX2.p,method = 'bonferroni')
  
  return(y_x.edge)
}

my_uniq<-function(x){
  length(unique(x))
}

my_mediation<-function(iv, m1, m2, dv, covDf){
  ## test
  #  iv<-all_phen$host_SmokeCurrentSmoker
  #  m1<-all_vsv$`Blautia wexlerae DSM 19850:2081_2082`
  # m2<-all_fmetab$cholate
  #  dv<-all_phen$IBD
  #  covDf<-input_cov
  
  require(lavaan)
  
  ## get input data frame
  input.df<-data.frame(iv,m1,m2,dv,covDf) %>% na.omit
  for (i in 1:ncol(input.df)) {
    if(length(unique(input.df[,i]))>2){
      input.df[,i]<-qtrans(input.df[,i])
    }
  }
  input.df<-as.data.frame(input.df)
  
  cov_str<-colnames(covDf)%>% paste(collapse = '+')
  
  ## Serial linear model
  model=paste(
    "
  #Regressions
  m1 ~ a*iv
  m2 ~ b*m1 + iv
  dv ~ c*m2 + m1 + d*iv + ",cov_str,
    "
  #Defined Parameters:
  ie := a*b*c
  de := d
  ", sep = ""
  )
  
  fit=sem(model,input.df)
  capture.output(fit.summ<-summary(fit),file='NUL')
  
  res_list <- list(indirect_effect.est = fit.summ$PE$est[grep("ie", fit.summ$PE$lhs)],
                   indirect_effect.se = fit.summ$PE$se[grep("ie", fit.summ$PE$lhs)],
                   indirect_effect.p = fit.summ$PE$pvalue[grep("ie", fit.summ$PE$lhs)],
                   direct_effect.est = fit.summ$PE$est[grep("de", fit.summ$PE$lhs)],
                   direct_effect.se = fit.summ$PE$se[grep("de", fit.summ$PE$lhs)],
                   direct_effect.p = fit.summ$PE$pvalue[grep("de", fit.summ$PE$lhs)])
  
  return(res_list)
}


serial_mediation<-function(inVec, iv_Df, m1_Df, m2_Df, dv_Df, covDf, covar ){
  ## test block
  #  inVec<-exp_sv_fmetab_df[i,]
  # iv_Df<-all_exp.sig
  # m1_Df<-all_sv.sig
  #  m2_Df<-all_fmetab.sig
  # dv_Df<-all_phen
  # covDf<-all_covar
  #  covar<-covar
  # test block
  
  if(is.na(inVec[5])){
    covar<-covar
  }else{
    covar<-c(covar, colnames(covDf)[grep(inVec[5],colnames(covDf))])
  }
  input_cov<-covDf[,covar]
  
  colnames(input_cov)[ncol(input_cov)]<-"abun"
  
  iv<-iv_Df[,match(inVec[1], colnames(iv_Df))]
  m1<-m1_Df[,match(inVec[2], colnames(m1_Df))]
  m2<-m2_Df[,match(inVec[3], colnames(m2_Df))]
  dv<-dv_Df[,match(inVec[4], colnames(dv_Df))]
  
  dir_res <- my_mediation(iv, m1, m2, dv, input_cov)
  
  variables<-list(iv=inVec[1], m1=inVec[2], m2=inVec[3],dv=inVec[4] )
  
  final_res<-c(variables,dir_res)
  
  return(final_res)
}



