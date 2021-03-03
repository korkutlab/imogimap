#' Finds combinatorial association of immunotherapy co-targets with tumor intrinsic features as listed in TCGA_immune_features_list.
#'
#' @import dplyr
#' @import curatedTCGAData
#' @param onco_gene A character vector of gene Hugo symbols.
#' @param icp_gene An optional character vector of immune checkpoint gene Hugo symbols.
#' @param cohort A list of TCGA diseases
#' @param sample_list An optional character vector of TCGA samples barcodes indicating a subset of samples within a cohort.
#' @param method A character string indicating which synergy score to be used. one of "max" or "independence". Default is "max".
#' @param data_feature  An optional numeric matrix or data frame containing normalized immune features.
#' @param pvalue An optional logical indicating if a random bootstrapping value should be calculated. Default is FALSE.
#' @param N_iteration Number of iterations for random bootstrapping for pvalue calculation and sensitivity analysis. Default is 1000.
#' @param sensitivity An optional logical indicating if a sensitivity analysis should be done. Default is FALSE.
#' @keywords Synergy scoring, immune feature, immune checkpoint, bootstrapping, TCGA
#' @return A dataframe of synergy scores and bootstrapping pvalues.
#' @details
#'
#' im_syng_tcga uses RNASeq2GeneNorm data from curatedTCGAData to find combinatorial association of immunotherapy co-targets and immune checkpoints with immuno-oncology features as listed in TCGA_immune_features_list
#'
#' For synergy score calculations all features are normalized to be on [0,1] range. For details of synergy score calculations see get_syng_score function.
#'
#' By default (if no icp_gene is specified), icp_gene_list will be used.
#'
#' All barcodes in sample_list must be 15 character long and belong to the same cohort. When sample_list is provided, cohort should be the disease cohort that they belong to, otherwise only the first element of the cohort list will be used.
#'
#' The optional data_feature must be normalized to have a range between [0,1].
#'
#' A p.value is computed using random bootstrapping with replacement from the distribution of synergy scores for each immune checkpoint and immune feature pair.
#'
#' Sensitivity analysis reports a sensitivity index which is a variance in synergy score using 90% of samples.
#' @examples im_syng_tcga(onco_gene=c("TP53","TGFB1"), cohort=c("acc","gbm"),add_pvalue=TRUE, N_iteration=1000)
#' @export

im_syng_tcga <- function(onco_gene, icp_gene, cohort, sample_list, method, feature, add_pvalue, N_iteration, sensitivity){

  df_syng <- data.frame(Disease=character(),
                        agent1=character(),
                        agent2=character(),
                        Immune_feature=character(),
                        Synergy_score=numeric(),
                        Pvalue=numeric(),
                        variance=numeric())

  #Check input parameters------------------------
  #----------------------------------------------

  if( missing( method ) ){
    method <- "max"
  }else{
    method <- tolower(method)
    if(!(method=="max" || method=="independence")){
      stop("ERROR: Method is not found. Please choose a method from: max or independence.")
    }
  }
  if(!missing(feature)){
    df_min <- min(feature,na.rm=T)
    df_max <- max(feature,na.rm=T)
    if(df_min < 0 | df_max > 1 ){
      stop("ERROR: feature is out of range. Normalize data_feature to [0,1].")
    }
    feature <- as.data.frame(feature)
    feature$Tumor_Sample_ID <- rownames(feauture)
  }
  if(!missing(sample_list)){
    cohort <- cohort[1]
  }
  cohort <- tolower(cohort)


  #Loop over diseases-----------------------------
  #-----------------------------------------------

  for(cohortID in 1:length(cohort)){

    #Get expression data----------------------------
    #-----------------------------------------------
    disease <- cohort[cohortID]
    message(disease,"\n")
    df <-curatedTCGAData::curatedTCGAData( diseaseCode = disease,
                                           assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]

    data_expression <- df@assays$data@listData[[1]]
    colnames(data_expression)<-  substr(colnames(data_expression), 1, 15)
    if(!missing(sample_list)){
      data_expression<-data_expression[,sample_list ]
      if(ncol(data_expression==0)){
        stop("ERROR: barcodes not found.")
      }
    }

    message("Ranking Gene expressions...")

    #Check for co-target expressions----------------
    #-----------------------------------------------
    if(length(onco_gene)==1){
      df_selected <- as.data.frame((data_expression[rownames(data_expression ) %in% onco_gene,]))
      colnames(df_selected) <- onco_gene
    }else{
      df_selected <- t(data_expression[rownames(data_expression ) %in% onco_gene,])
    }
    if(nrow(df_selected)==0){
      stop("ERROR: No Hugo symbols found for onco-genes.")
    }
    df_selected <- as.data.frame(df_selected[, colSums(df_selected != 0) > 0,drop=FALSE])
    if(ncol(df_selected)==0){
      warning("All onco_gene's have zero expression in ", disease )
      next
    }
    onco_gene_sub <- colnames(df_selected)

    #Check for immune checkpoint expressions--------
    #-----------------------------------------------
    if(missing(icp_gene)){
      icp_gene <- icp_gene_list
    }
    if(length(icp_gene)==1){
      df_icp <- as.data.frame(data_expression[rownames(data_expression) %in% icp_gene,])
      colnames(df_icp) <- icp_gene
    }else{
      df_icp <- t(data_expression[rownames(data_expression) %in% icp_gene,])
    }
    if(nrow(df_icp)==0){
      stop("ERROR: No Hugo symbols found for icp_genes.")
    }
    df_icp <- as.data.frame(df_icp[, colSums(df_icp != 0) > 0,drop=FALSE])
    if(ncol(df_icp)==0){
      warning("All icp_gene's have zero expression in ", disease )
      next
    }
    icp_gene_sub <- colnames(df_icp)

    #Construct quantile ranking matrices for each sample
    #-----------------------------------------------
    df_selected <- scale(log2(df_selected+1),center = T,scale = T)
    df_icp <- scale(log2(df_icp+1),center = T,scale = T)
    df_select_qr <- get_quantile_rank(df_selected)
    df_icp_qr <- get_quantile_rank(df_icp)
    df_all <- merge(df_select_qr,df_icp_qr)

    #Construct quantile ranking matrices for each patient.
    #For immune cell count features as calculated by CIBERSORT,
    #PATIENT_BARCODE is used instead of Tumor_Sample_ID
    #-----------------------------------------------
    df_select_qr2 <- df_select_qr
    df_select_qr2$PATIENT_BARCODE <- substr(df_select_qr2$Tumor_Sample_ID, 1, 12)
    df_select_qr2$Tumor_Sample_ID <- NULL
    df_select_qr2 <- df_select_qr2 %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
    df_select_qr2 <- df_select_qr2[,c("PATIENT_BARCODE",
                                      c(setdiff(colnames(df_select_qr2), "PATIENT_BARCODE")))]

    df_icp_qr2 <- df_icp_qr
    df_icp_qr2$PATIENT_BARCODE <- substr(df_icp_qr2$Tumor_Sample_ID, 1, 12)
    df_icp_qr2$Tumor_Sample_ID <- NULL
    df_icp_qr2 <- df_icp_qr2 %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
    df_icp_qr2 <- df_icp_qr2[,c("PATIENT_BARCODE",
                                c(setdiff(colnames(df_icp_qr2), "PATIENT_BARCODE")))]
    df_all2 <- merge(df_select_qr2,df_icp_qr2)

    message("Calculating features...")

    #Get features-----------------------------------
    #-----------------------------------------------

    cohort_EMT <- get_emt_score(data_expression)
    cohort_IFNG <- get_ifng_score(data_expression)
    cohort_AG <- get_angio_score(data_expression)

    #Normalize features-----------------------------
    #-----------------------------------------------
    sd_EMT <- sd(cohort_EMT$EMTscore,na.rm = T)
    sd_IFNG <- sd(cohort_IFNG$IFNGscore,na.rm = T)
    sd_AG <- sd(cohort_AG$AGscore,na.rm = T)
    cohort_EMT$EMTscore <- ( tanh( cohort_EMT$EMTscore/sd_EMT ) + 1 ) / 2
    cohort_IFNG$IFNGscore <- ( tanh( cohort_IFNG$IFNGscore/sd_IFNG ) + 1 ) / 2
    cohort_AG$AGscore <- ( tanh( cohort_AG$AGscore/sd_AG ) + 1 ) / 2
    df_lf <- TCGA_TMB
    df_lf$TMB_Non.silent_per_Mb <- tanh(  df_lf$TMB_Non.silent_per_Mb/10 )
    df_lf$TMB_Silent_per_Mb <- tanh(  df_lf$TMB_Silent_per_Mb/10 )

    #Merge all features-----------------------------
    #-----------------------------------------------
    dft <- as.data.frame(merge(cohort_EMT , cohort_AG, by="Tumor_Sample_ID"))
    dft <- as.data.frame(merge(dft , cohort_IFNG, by="Tumor_Sample_ID"))
    dft<- as.data.frame(merge(dft , TCGA_Leukocyte_fraction, by="Tumor_Sample_ID"))
    dft <- as.data.frame(merge(dft , df_lf, by="Tumor_Sample_ID"))

    #Add optional user provided features------------
    #-----------------------------------------------
    if(!missing(feature)){
      dft <- as.data.frame(merge(dft , feature, by="Tumor_Sample_ID"))
    }
    data_feature <- dft
    rm(dft)

    #Build unique permutations----------------------
    #-----------------------------------------------
    all_genes <- unique(c(onco_gene_sub,icp_gene_sub))
    all_perms <- t(combn(all_genes,m = 2))
    N_perm <- as.numeric(nrow(all_perms))
    N_genes <- as.numeric(length(all_genes))
    N_feature <- as.numeric(ncol(data_feature))-1
    N_imcell_feature <- as.numeric(ncol(TCGA_IMCell_fraction))-1


    #Calculate synergy scores for features ---------
    #-----------------------------------------------
    message("Calculating synerrgy scores for immune features ...")
    for(pair_ID in 1:N_perm){
      gene_ID1 <- which(colnames(df_all)==all_perms[pair_ID,1])
      gene_ID2 <- which(colnames(df_all)==all_perms[pair_ID,2])
      dft <- df_all[ , c(1 , gene_ID1,gene_ID2)]
      dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
      dft <- dft[dft[ , 3] %in% c(1 , 4) , ]

      for(im_ID in 2:ncol(data_feature)){

        dft2 <- merge(data_feature[ , c(1 , im_ID)] , dft , by="Tumor_Sample_ID")
        dft2 <- dft2[complete.cases(dft2),]
        dft2 <- dft2[,c(2,3,4)]
        dfts <- get_syng_score(dft2,method)
        dfts$Disease <- disease
        dfts <- dfts[ , c("Disease" , c(setdiff(colnames(dfts) , "Disease")))]
        dfts$Pvalue <- NA
        dfts$variance <- NA
        df_syng <- rbind(df_syng , dfts)

      }
    }
    message("Calculating synerrgy scores for immune cell fractions ...")
    #Calculate synergy scores for immune cell features
    #-----------------------------------------------
    for(pair_ID in 1:N_perm){
      gene_ID1 <- which(colnames(df_all2)==all_perms[pair_ID,1])
      gene_ID2 <- which(colnames(df_all2)==all_perms[pair_ID,2])
      dft <- df_all2[ , c(1 , gene_ID1,gene_ID2)]
      dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
      dft <- dft[dft[ , 3] %in% c(1 , 4) , ]

      for(if_ID in 2:ncol(TCGA_IMCell_fraction)){
        dft2 <- merge(TCGA_IMCell_fraction[,c(1,if_ID)],dft,by="PATIENT_BARCODE")
        dft2 <- dft2[complete.cases(dft2),]
        dft2<- dft2[,c(2,3,4)]
        dfts <- get_syng_score(dft2,method)
        dfts$Disease <- disease
        dfts <- dfts[ , c("Disease" , c(setdiff(colnames(dfts) , "Disease")))]
        dfts$Pvalue <- NA
        dfts$variance <- NA
        df_syng <- rbind(df_syng,dfts)
      }
    }
    message("Synergy calculation completed! ")

    #Check if pvalue should be calculated
    #-----------------------------------------------
    if(missing(add_pvalue)){
      add_pvalue <- FALSE
    }
    if(add_pvalue){

      if(missing(N_iteration)){
        N_iteration<-1000
      }

      message("Building bootstrapping distribution from ", N_iteration," genes")

      #Build unique permutations----------------------
      #-----------------------------------------------
      df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
      df_syng_complete <- df_syng_complete[ which(df_syng_complete$Disease==disease), ]
      all_genes <- unique( c( df_syng_complete$agent1,df_syng_complete$agent2))
      N_genes <- as.numeric( length( all_genes))
      sub_feature_list <- unique(df_syng_complete[ df_syng_complete$Immune_feature %in% colnames( data_feature),]$Immune_feature)
      N_feature <- as.numeric( length( sub_feature_list))
      sub_imcell_feature_list <- unique(df_syng_complete[ df_syng_complete$Immune_feature %in% colnames( TCGA_IMCell_fraction ),]$Immune_feature)
      N_imcell_feature <-as.numeric( length( sub_imcell_feature_list))

      #Build random bank of expression values---------
      #-----------------------------------------------
      df_bank <- data.frame()
      while(ncol(df_bank)==0){
        df_bank <- as.data.frame(t(data_expression[
          c(sample(nrow(data_expression),N_iteration,replace = T),
            which(rownames(data_expression) %in% all_genes)),]))
        df_bank <-  df_bank[, colSums(abs(df_bank)!= 0 ) > 0]
        df_bank <- df_bank[!duplicated(as.list(df_bank))]
      }
      df_bank <- df_bank[,c(all_genes, setdiff(colnames(df_bank),all_genes))]
      #Construct quantile ranking matrices for bank---
      #-----------------------------------------------
      df_bank <- scale(log2(df_bank+1),center = T,scale = T)
      df_bank_qr <- get_quantile_rank(df_bank)
      df_bank_qr2 <- df_bank_qr
      df_bank_qr2$PATIENT_BARCODE <- substr(df_bank_qr2$Tumor_Sample_ID, 1, 12)
      df_bank_qr2$Tumor_Sample_ID <- NULL
      df_bank_qr2 <- df_bank_qr2 %>% group_by(PATIENT_BARCODE) %>%
        mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
      df_bank_qr2 <- df_bank_qr2[,c("PATIENT_BARCODE",
                                    c(setdiff(colnames(df_bank_qr2), "PATIENT_BARCODE")))]


      #Calculate p.values for each feature
      #-----------------------------------------------
      message("Calculating pvalues...")

      for( i in 1:N_feature){

        #Select feature
        select_feature <- sub_feature_list[i]
        mark_feature <- as.integer( which( colnames( data_feature) == select_feature))
        df_feature <- data_feature[ , c( 1, mark_feature)]
        message( "          ...",select_feature)

        #Make a list of all genes with non-NA synergy score value
        df_syng_t <- df_syng_complete[ which( df_syng_complete$Immune_feature == select_feature), ]
        sub_pairs <- unique( df_syng_t[ c( "agent1" ,"agent2")])
        sub_genes <- unique( c( sub_pairs$agent1, sub_pairs$agent2))

        syng_dist <- vector( "list", length( sub_genes))
        names(syng_dist) <- sub_genes

        df_bank_sub1 <- merge( df_feature, df_bank_qr, by="Tumor_Sample_ID")
        bank_size <- as.numeric(ncol(df_bank_sub1))
        bank_mark <- N_genes+3

        #Build bootstrapping distribution --------------
        #-----------------------------------------------
        for( j in 1 :length( sub_genes)){

          gene_mark <- as.numeric( which( colnames( df_bank_sub1) == sub_genes[ j]))
          df_bank_sub2 <- df_bank_sub1[ , c(1, gene_mark, bank_mark : bank_size)]
          df_bank_sub2 <- df_bank_sub2[ df_bank_sub2[ , ncol( df_bank_sub2) ] %in% c( 1, 4), ]
          bank_size2 <- as.numeric( ncol( df_bank_sub2))

          syng_dist_t <- vector()
          for(k in 3:(bank_size2-1)){
            dft <- df_bank_sub2[,c(1,2,k,bank_size2)]
            dft <- dft[dft[,3] %in% c(1,4),]
            dft <- dft[complete.cases(dft),]
            dft<- dft[,c(2,3,4)]
            syng_dist_t[k-2] <- get_syng_score(dft,method)$Synergy_score
          }
          syng_dist_t <- syng_dist_t[complete.cases(syng_dist_t)]
          syng_dist[[j]] <- syng_dist_t
        }

        # Loop over all pairs---------------------------
        #-----------------------------------------------
        for(pair_ID in 1:nrow(df_syng_t)){
          gene1 <- df_syng_t$agent1[pair_ID]
          gene2 <- df_syng_t$agent2[pair_ID]
          gene_ID1 <- which(colnames(df_all)==gene1)
          gene_ID2 <- which(colnames(df_all)==gene2)
          dft <- df_all[ , c(1 , gene_ID1,gene_ID2)]
          dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
          dft <- dft[dft[ , 3] %in% c(1 , 4) , ]
          dft <- merge(df_feature , dft , by="Tumor_Sample_ID")
          dft <- dft[complete.cases(dft),]
          dft <- dft[,c(2,3,4)]

          myscore <- df_syng_t$Synergy_score[pair_ID]
          mysign <- sign(myscore)
          if( mysign==0) mysign <- 1

          #Get pvalue-------------------------------------
          #-----------------------------------------------
          j1 <- which( names( syng_dist) == gene1)
          tmp_dist <- syng_dist[[j1]]
          Pvalue1 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
          j2 <- which(names(syng_dist)==gene2)
          tmp_dist <- syng_dist[[j2]]
          Pvalue2 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
          cc <- which( df_syng$Immune_feature == select_feature &
                         df_syng$agent1==gene1 &
                         df_syng$agent2 == gene2 &
                         df_syng$Disease==disease)
          df_syng$Pvalue[cc] <- max(Pvalue1, Pvalue2)

        }
      }


      #Calculate p.value for each immune cell counts
      #-----------------------------------------------
      for( i in 1:N_imcell_feature){

        #Select feature
        select_feature <- sub_imcell_feature_list[i]
        mark_feature <- as.integer( which( colnames( TCGA_IMCell_fraction) == select_feature))
        df_feature <- TCGA_IMCell_fraction[ , c( 1, mark_feature)]
        message( "          ...",select_feature)

        #Make a list of all genes with non-NA synergy score value
        df_syng_t <- df_syng_complete[ which( df_syng_complete$Immune_feature == select_feature), ]
        sub_pairs <- unique( df_syng_t[ c( "agent1" ,"agent2")])
        sub_genes <- unique( c( sub_pairs$agent1, sub_pairs$agent2))

        syng_dist <- vector( "list", length( sub_genes))
        names(syng_dist) <- sub_genes

        df_bank_sub1 <- merge( df_feature, df_bank_qr2, by="PATIENT_BARCODE")
        bank_size <- as.numeric(ncol(df_bank_sub1))
        bank_mark <- N_genes+3

        #Build bootstrapping distribution --------------
        #-----------------------------------------------
        for( j in 1 :length( sub_genes)){

          gene_mark <- as.numeric( which( colnames( df_bank_sub1) == sub_genes[ j]))
          df_bank_sub2 <- df_bank_sub1[ , c(1, gene_mark, bank_mark : bank_size)]
          df_bank_sub2 <- df_bank_sub2[ df_bank_sub2[ , ncol( df_bank_sub2) ] %in% c( 1, 4), ]
          bank_size2 <- as.numeric( ncol( df_bank_sub2 ))

          syng_dist_t <- vector()
          for(k in 3:(bank_size2-1)){
            dft <- df_bank_sub2[,c(1,2,k,bank_size2)]
            dft <- dft[dft[,3] %in% c(1,4),]
            dft <- dft[complete.cases(dft),]
            dft<- dft[,c(2,3,4)]
            syng_dist_t[k-2] <- get_syng_score(dft,method)$Synergy_score
          }
          syng_dist_t <- syng_dist_t[complete.cases(syng_dist_t)]
          syng_dist[[j]] <- syng_dist_t
        }

        # Loop over all pairs---------------------------
        #-----------------------------------------------
        for(pair_ID in 1:nrow(df_syng_t)){
          gene1 <- df_syng_t$agent1[pair_ID]
          gene2 <- df_syng_t$agent2[pair_ID]
          gene_ID1 <- which(colnames(df_all)==gene1)
          gene_ID2 <- which(colnames(df_all)==gene2)
          dft <- df_all2[ , c(1 , gene_ID1,gene_ID2)]
          dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
          dft <- dft[dft[ , 3] %in% c(1 , 4) , ]
          dft <- merge(df_feature , dft , by="PATIENT_BARCODE")
          dft <- dft[complete.cases(dft),]
          dft <- dft[,c(2,3,4)]

          myscore <- df_syng_t$Synergy_score[pair_ID]
          mysign <- sign(myscore)
          if( mysign==0) mysign <- 1

          #Get pvalue-------------------------------------
          #-----------------------------------------------
          j1 <- which( names( syng_dist) == gene1)
          tmp_dist <- syng_dist[[j1]]
          Pvalue1 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
          j2 <- which(names(syng_dist)==gene2)
          tmp_dist <- syng_dist[[j2]]
          Pvalue2 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
          cc <- which( df_syng$Immune_feature == select_feature &
                         df_syng$agent1==gene1 &
                         df_syng$agent2 == gene2 &
                         df_syng$Disease==disease)
          df_syng$Pvalue[cc] <- max(Pvalue1, Pvalue2)

        }
      }
    }

    #Check if sensitivity analysis should be done---
    #-----------------------------------------------
    if(missing(sensitivity)){
      sensitivity <- FALSE
    }
    if(sensitivity){
      if(missing(N_iteration)){
        N_iteration <- 1000
      }
      #Define parameters
      df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
      df_syng_complete$sum <- 0
      df_syng_complete$sum2 <- 0
      df_syng_complete$N <- 0


      message("Calculating sensitivity variance. Iteratiing ", N_iteration, " times.")


      dft1 <- as.data.frame(df_selected)
      dft2 <- as.data.frame(df_icp)
      dft1$ID <- rownames(dft1)
      dft2$ID <- rownames(dft2)
      df_comb <- merge(dft1,dft2)
      rownames(df_comb) <- df_comb$ID
      df_comb$ID <- NULL
      df_comb <- as.matrix(df_comb)
      N_sample <- as.numeric(nrow(df_comb))

      N_sub <- floor(N_sample*0.7)

      df_syng_complete1 <- df_syng_complete[ which(df_syng_complete$Disease==disease &
                                                     df_syng_complete$Immune_feature %in%
                                                     colnames(data_feature)),]
      df_syng_complete2 <- df_syng_complete[ which(df_syng_complete$Disease==disease &
                                                     df_syng_complete$Immune_feature %in%
                                                     colnames(TCGA_IMCell_fraction)),]
      N_syng_complete1 <- nrow(df_syng_complete1)
      N_syng_complete2 <- nrow(df_syng_complete2)


      for(n_sns in 1:N_iteration){

        message("Iteration ", n_sns)

        df_sub <- df_comb[sample(N_sample,N_sub,replace = F),]
        df_sub_qr <- get_quantile_rank(df_sub)

        df_sub_qr2 <-  df_sub_qr
        df_sub_qr2$PATIENT_BARCODE <- substr(df_sub_qr2$Tumor_Sample_ID, 1, 12)
        df_sub_qr2$Tumor_Sample_ID <- NULL
        df_sub_qr2 <- df_sub_qr2 %>% group_by(PATIENT_BARCODE) %>%
          mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
        df_sub_qr2 <- df_sub_qr2[,c("PATIENT_BARCODE",
                                      c(setdiff(colnames(df_sub_qr2), "PATIENT_BARCODE")))]
        df_sub_qr2 <- as.data.frame(df_sub_qr2)

        for(pair_ID in 1:N_syng_complete1 ){
          gene_ID1 <- df_syng_complete1$agent1[pair_ID]
          gene_ID2 <- df_syng_complete1$agent2[pair_ID]
          feature_ID <- df_syng_complete1$Immune_feature[pair_ID]
          mark_feature <- as.integer( which( colnames( data_feature) == feature_ID))
          df_feature <- data_feature[ , c( 1, mark_feature)]

          dft <- df_sub_qr[ , c(1 ,
                                which(colnames(df_sub_qr)==gene_ID1),
                                which(colnames(df_sub_qr)==gene_ID2))]
          dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
          dft <- dft[dft[ , 3] %in% c(1 , 4) , ]
          dft <- merge(df_feature , dft , by="Tumor_Sample_ID")
          dft <- dft[complete.cases(dft),]
          dft <- dft[,c(2,3,4)]
          dfts <- get_syng_score(dft,method)$Synergy_score

          df_syng_complete1$sum[pair_ID] <- sum(df_syng_complete1$sum[pair_ID], dfts, na.rm = T)
          df_syng_complete1$sum2[pair_ID] <-sum(df_syng_complete1$sum2[pair_ID], dfts*dfts, na.rm = T)
          df_syng_complete1$N[pair_ID] <- sum(df_syng_complete1$N[pair_ID], abs(sign(dfts)), na.rm = T)
        }
        for(pair_ID in 1:N_syng_complete2 ){
          gene_ID1 <- df_syng_complete2$agent1[pair_ID]
          gene_ID2 <- df_syng_complete2$agent2[pair_ID]
          feature_ID <- df_syng_complete2$Immune_feature[pair_ID]
          mark_feature <- as.integer( which( colnames(TCGA_IMCell_fraction) == feature_ID))
          df_feature <- TCGA_IMCell_fraction[ , c( 1, mark_feature)]

          dft <- df_sub_qr2[ , c(1 ,
                                which(colnames(df_sub_qr2)==gene_ID1),
                                which(colnames(df_sub_qr2)==gene_ID2))]
          dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
          dft <- dft[dft[ , 3] %in% c(1 , 4) , ]
          dft <- merge(df_feature , dft , by="PATIENT_BARCODE")
          dft <- dft[complete.cases(dft),]
          dft <- dft[,c(2,3,4)]
          dfts <- get_syng_score(dft,method)$Synergy_score
          df_syng_complete2$sum[pair_ID] <- sum(df_syng_complete2$sum[pair_ID], dfts, na.rm = T)
          df_syng_complete2$sum2[pair_ID] <-sum(df_syng_complete2$sum2[pair_ID], dfts*dfts, na.rm = T)
          df_syng_complete2$N[pair_ID] <- sum(df_syng_complete2$N[pair_ID], abs(sign(dfts)), na.rm = T)

        }
      }
      df_syng_complete <- rbind(df_syng_complete1,df_syng_complete2)
      df_syng_complete$variance <-
      df_syng_complete$sum2 / df_syng_complete$N - (df_syng_complete$sum/df_syng_complete$N)^2
      df_syng_complete$sum <- NULL
      df_syng_complete$sum2 <- NULL
      df_syng_complete$N <- NULL
      rownames(df_syng_complete) <- NULL
      df_syng$variance <- NULL
      df_syng <- merge(df_syng,df_syng_complete,all=T)
    }

    message(disease," completed\n")
  }
  message("All analysis done!")
  return(df_syng)
}
