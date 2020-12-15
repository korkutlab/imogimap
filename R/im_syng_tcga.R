#' Find combinatorial association of immunotherapy co-targets with all tumor intrinsic features.
#'
#' @import dplyr
#' @import cBioPortalData
#' @param cotarget A character vector of gene Hugo symbols.
#' @param checkpoint An optional character vector of immune checkpoint gene Hugo symbols.
#' @param cohort a list of TCGA diseases
#' @param method a charachter string indicating which synergy score to be used. one of "HSA" or "Bliss". Default is HSA.
#' @param data_feature  An optional numeric matrix or data frame containing normalized immune features.
#' @param add_pvalue An optional logical indicating if a random bootstrapping value should be calculated. Default is FALSE.
#' @param N_iteration Number of iterrations for random bootstrapping
#' @keywords Synergy scoring, immune feature, immune checkpoint, bootstrapping, TCGA
#' @return a dataframe of synergy scores and bootstrapping pvalues
#' @details
#'
#' im_syng_tcga uses gene expressions from cbioportal data, 2018  tcga pancan atlas to find combinatorial association of immunotherapy co-targets and immune checkpoints with immuno-oncology features as listed in TCGA_immune_features_list
#'
#' For details of synergy score calculations see get_syng_score function.
#'
#' By default (if no checkpoint is specified), icp_gene_list will be used.
#'
#' The optional data_feature must be normalized to have a range between [0,1].
#'
#' A p.value is computed using random bootstraping with replacement from the distibution of synergy scores for each immune checkpointand-immune feature pair. The default values of N_iteration is 1000.
#'
#' @examples im_syng_tcga(cotarget=c("TP53","TGFB1"),checkpoint=c(),cohort=c("acc","gbm"),add_pvalue=TRUE, N_iteration=1000)
#' @export

im_syng_tcga<-function(cotarget, checkpoint, cohort, method, feature, add_pvalue, N_iteration){

  #Check inputs
  if(missing(method)){
    method <- "HSA"
  }else{
    method <- toupper(method)
    if(!(method=="HSA" || method=="BLISS")){
      stop("ERROR: Method is not HSA or BLISS.")
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

  cohort <- tolower(cohort)

  for(cohortID in 1:length(cohort)){

    #Read data------------------------------------
    disease <- cohort[cohortID]
    cohort_study <- paste0(disease,"_tcga_pan_can_atlas_2018")
    df <- cBioPortalData::cBioDataPack(cohort_study,ask = F)@ExperimentList@listData
    df <- df$RNA_Seq_v2_expression_median
    df2 <- df@elementMetadata@listData
    data_expression <- df@assays@data@listData[[1]]
    rownames(data_expression)<- plyr::mapvalues(rownames(data_expression),df2$Entrez_Gene_Id,df2$Hugo_Symbol,
      warn_missing = F)
    rm(df,df2)

    #Check for co-target expressions---------------
    if(length(cotarget)==1){
      df_selected <- as.data.frame((data_expression[rownames(data_expression ) %in% cotarget,]))
      colnames(df_selected) <- cotarget
    }else{
      df_selected <- t(data_expression[rownames(data_expression ) %in% cotarget,])
    }
    if(nrow(df_selected)==0){
      stop("ERROR: No valid Hugo symbol found")
    }

    #Check for immune checkpoint expressions--------
    if(missing(checkpoint)){
      checkpoint <- icp_gene_list
    }
    if(length(checkpoint)==1){
      df_icp <- as.data.frame(data_expression[rownames(data_expression) %in% checkpoint,])
      colnames(df_icp) <- checkpoint
    }else{
      df_icp <- t(data_expression[rownames(data_expression) %in% checkpoint,])
    }
    #Construct quantile ranking matrices for each sample -----------------------
    df_selected <- scale(log2(df_selected+1),center = T,scale = T)
    df_icp <- scale(log2(df_icp+1),center = T,scale = T)
    df_select_qr <- get_qunatile_rank(df_selected)
    df_icp_qr <- get_qunatile_rank(df_icp)

    #Construct quantile ranking matrices for each patient -------------------------------------
    #PATIENT_BARCODE is used instead of Tumor_Sample_ID for immune cell count features as calculated by CIBERSORT
    df_select_qr2 <- df_select_qr
    df_select_qr2$PATIENT_BARCODE <- substr(df_select_qr2$Tumor_Sample_ID, 1, 12)
    df_select_qr2$Tumor_Sample_ID <- NULL
    df_select_qr2 <- df_select_qr2 %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
    df_select_qr2 <- df_select_qr2[,c("PATIENT_BARCODE",
      c(setdiff(colnames(df_select_qr2), "PATIENT_BARCODE")))]

    df_icp_qr2 <- df_icp_qr
    df_icp_qr2$PATIENT_BARCODE <- substr(df_icp_qr2$Tumor_Sample_ID, 1, 12)
    df_icp_qr2$Tumor_Sample_ID <- NULL
    df_icp_qr2 <- df_icp_qr2 %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
    df_icp_qr2 <- df_icp_qr2[,c("PATIENT_BARCODE",
      c(setdiff(colnames(df_icp_qr2), "PATIENT_BARCODE")))]


    #Get features--------------------------------------
    cohort_EMT <- get_emt_score(data_expression)
    cohort_IFNG <- get_ifng_score(data_expression)
    cohort_AG <- get_angio_score(data_expression)

    #Normalize features--------------------------------------
    cohort_EMT$EMTscore <- ( tanh( cohort_EMT$EMTscore ) + 1 ) / 2
    cohort_IFNG$IFNGscore <- ( tanh( cohort_IFNG$IFNGscore ) + 1 ) / 2
    cohort_AG$AGscore <- ( tanh( cohort_AG$AGscore ) + 1 ) / 2
    df_lf <- TCGA_TMB
    df_lf$TMB_Non.silent_per_Mb <- tanh(  df_lf$TMB_Non.silent_per_Mb )
    df_lf$TMB_Silent_per_Mb <- tanh(  df_lf$TMB_Silent_per_Mb )

    #Merge all features-------------------------------
    dft <- as.data.frame(merge(cohort_EMT , cohort_AG, by="Tumor_Sample_ID"))
    dft <- as.data.frame(merge(dft , cohort_IFNG, by="Tumor_Sample_ID"))
    dft<- as.data.frame(merge(dft , TCGA_Leukocyte_fraction, by="Tumor_Sample_ID"))
    dft <- as.data.frame(merge(dft , df_lf, by="Tumor_Sample_ID"))

    #Add optional user provided features---------------
    if(!missing(feature)){
      dft <- as.data.frame(merge(dft , feature, by="Tumor_Sample_ID"))
    }
    data_feature <- dft
    rm(dft)

    #Check if bootstrapping pvalue should be calculated----------------------
    if(missing(add_pvalue)){
      add_pvalue <- FALSE
    }
    if(add_pvalue){
      if(missing(N_iteration)){
        N_iteration<-1000
      }
      df_syng <- data.frame(Disease=character(), Co_target=character(),
        Immune_checkpoint=character(),
        Immune_feature=character(),
        Synergy_score=numeric(),
        pvalue=numeric())

      #Build random bank of expression values-----------------------------------
      df_bank <- as.data.frame(t(data_expression[
        c(sample(nrow(data_expression),N_iteration,replace = T),
          which(rownames(data_expression) %in% checkpoint)),]))
      df_bank <- df_bank[, -which(colSums(abs(df_bank)) == 0 )]

      #Construct quantile ranking matrices for bank---------------------
      df_bank <- scale(log2(df_bank+1),center = T,scale = T)
      df_bank_qr <- get_qunatile_rank(df_bank)
      df_bank_qr2 <- df_bank_qr
      df_bank_qr2$PATIENT_BARCODE <- substr(df_bank_qr2$Tumor_Sample_ID, 1, 12)
      df_bank_qr2$Tumor_Sample_ID <- NULL
      df_bank_qr2 <- df_bank_qr2 %>% group_by(PATIENT_BARCODE) %>%
        mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
      df_bank_qr2 <- df_bank_qr2[,c("PATIENT_BARCODE",
        c(setdiff(colnames(df_bank_qr2), "PATIENT_BARCODE")))]

      #Merge bank with feature data---------------------------------
      df_bank_merged <- merge(data_feature,df_bank_qr,by="Tumor_Sample_ID")
      df_bank_merged <- as.data.frame(df_bank_merged)
      df_bank_merged2 <- merge(TCGA_IMCell_fraction,df_bank_qr2,by="PATIENT_BARCODE")
      df_bank_merged2 <- as.data.frame(df_bank_merged2)


      #Loop over checkpoints-------------------------------------------------
      N_checkpoint <- as.numeric(length(checkpoint))
      N_feature <- as.numeric(ncol(data_feature))-1
      N_imcell_feature <- as.numeric(ncol(TCGA_IMCell_fraction))-1

      for( i in 1:N_checkpoint){
        checkpoint_mark <- as.numeric(ncol(df_bank_merged))-N_checkpoint+i
        df_bank_sub1 <- df_bank_merged[,c(1:(checkpoint_mark -i),checkpoint_mark )]
        df_bank_sub1 <- df_bank_sub1[df_bank_sub1[,ncol(df_bank_sub1) ] %in% c(1,4),]

        #Loop over features--------------------------------
        for( j in 1:N_feature){
          df_bank_sub2 <- df_bank_sub1[,c(1,(j+1),(N_feature+2):ncol(df_bank_sub1))]
          bank_size <- as.numeric(ncol(df_bank_sub2))

          #Build bootstrapping distribution-----------
          syng_dist <- vector()
          for(k in 3:(bank_size-1)){
            dft <- df_bank_sub2[,c(1,2,k,bank_size)]
            dft <- dft[dft[,3] %in% c(1,4),]
            dft <- dft[complete.cases(dft),]
            dft<- dft[,c(2,3,4)]
            syng_dist[k-2] <- get_syng_score(dft,method)$Synergy_score
          }
          syng_dist <- syng_dist[complete.cases(syng_dist)]

          # Loop over co-targets---------
          for(cotarget_ID in 2:ncol(df_select_qr)){
            dft <- merge(df_select_qr[,c(1,cotarget_ID)],df_icp_qr[,c(1,i+1)],by="Tumor_Sample_ID")
            dft <- merge(data_feature[,c(1,j+1)],dft,by="Tumor_Sample_ID")
            dft <- dft[dft[,3] %in% c(1,4),]
            dft <- dft[dft[,4] %in% c(1,4),]
            dft <- dft[complete.cases(dft),]
            dft<- dft[,c(2,3,4)]
            cc <- get_syng_score(dft,method)
            #Get synergy score--------
            myscore <- cc$Synergy_score
            if(is.na(myscore)){
              cc$Pvalue <- NA
            }else{
              mysign <- sign(myscore)
              if(mysign==0) mysign<-1
              #Get pvalue---------
              cc$Pvalue <- sum(mysign*syng_dist > mysign*myscore, syng_dist==myscore )/length(syng_dist)
            }
            cc$Disease <- disease
            cc <- cc[ , c("Disease" , c(setdiff(colnames(cc) , "Disease")))]
            df_syng <- rbind(df_syng,cc)
          }
        }

        checkpoint_mark <- as.numeric(ncol(df_bank_merged2))-N_checkpoint+i
        df_bank_sub1 <- df_bank_merged2[,c(1:(checkpoint_mark -i),checkpoint_mark )]
        df_bank_sub1 <- df_bank_sub1[df_bank_sub1[,ncol(df_bank_sub1) ] %in% c(1,4),]

        #Loop over immune cell counts------------
        for( j in 1:N_imcell_feature){
          df_bank_sub2 <- df_bank_sub1[,c(1,(j+1),(N_imcell_feature+2):ncol(df_bank_sub1))]
          bank_size <- as.numeric(ncol(df_bank_sub2))

          #Build bootstrapping distribution-----------
          syng_dist <- vector()
          for(k in 3:(bank_size-1)){
            dft <- df_bank_sub2[,c(1,2,k,bank_size)]
            dft <- dft[dft[,3] %in% c(1,4),]
            dft <- dft[complete.cases(dft),]
            dft<- dft[,c(2,3,4)]
            syng_dist[k-2] <- get_syng_score(dft,method)$Synergy_score
          }
          syng_dist <- syng_dist[complete.cases(syng_dist)]

          #Loop over co-targets--------
          for(cotarget_ID in 2:ncol(df_select_qr2)){
            dft <- merge(df_select_qr2[,c(1,cotarget_ID)],df_icp_qr2[,c(1,i+1)],by="PATIENT_BARCODE")
            dft <- merge(TCGA_IMCell_fraction[,c(1,j+1)],dft,by="PATIENT_BARCODE")
            dft <- dft[dft[,3] %in% c(1,4),]
            dft <- dft[dft[,4] %in% c(1,4),]
            dft <- dft[complete.cases(dft),]
            dft<- dft[,c(2,3,4)]
            #Get synergy score--------
            cc <- get_syng_score(dft,method)
            myscore <- cc$Synergy_score
            if(is.na(myscore)){
              cc$Pvalue <- NA
            }else{
              mysign <- sign(myscore)
              if(mysign==0) mysign<-1
              #Get pvalue----------
              cc$Pvalue <- sum(mysign*syng_dist > mysign*myscore, syng_dist==myscore )/length(syng_dist)
            }
            cc$Disease <- disease
            cc <- cc[ , c("Disease" , c(setdiff(colnames(cc) , "Disease")))]
            df_syng <- rbind(df_syng,cc)
          }
        }
      }
    }else{

      df_syng <- data.frame(Disease=character(),
        Co_target=character(),
        Immune_checkpoint=character(),
        Immune_feature=character(),
        Synergy_score=numeric())

      #Calculate synergy scores for features -----------------------
      for(gene_ID in 2:ncol(df_select_qr)){
        for(icp_ID in 2:ncol(df_icp_qr)){

          dft <- merge(df_select_qr[ , c(1 , gene_ID)] , df_icp_qr[ , c(1,icp_ID)] ,
            by="Tumor_Sample_ID")

          for(im_ID in 2:ncol(data_feature)){

            dft2 <- merge(data_feature[ , c(1 , im_ID)] , dft , by="Tumor_Sample_ID")
            dft2 <- dft2[dft2[ , 3] %in% c(1 , 4) , ]
            dft2 <- dft2[dft2[ , 4] %in% c(1 , 4) , ]
            dft2 <- dft2[complete.cases(dft2),]

            dft2<- dft2[,c(2,3,4)]
            dfts <- get_syng_score(dft2,method)
            dfts$Disease <- disease
            dfts <- dfts[ , c("Disease" , c(setdiff(colnames(dfts) , "Disease")))]
            df_syng <- rbind(df_syng , dfts)

          }
        }
      }

      #Calculate synergy scores for immune cell features -----------------------
      for(gene_ID in 2:ncol(df_select_qr2)){
        for(icp_ID in 2:ncol(df_icp_qr2)){
          dft <- merge(df_select_qr2[,c(1,gene_ID)],df_icp_qr2[,c(1,icp_ID)],by="PATIENT_BARCODE")
          for(if_ID in 2:ncol(TCGA_IMCell_fraction)){
            dft2 <- merge(TCGA_IMCell_fraction[,c(1,if_ID)],dft,by="PATIENT_BARCODE")
            dft2 <- dft2[dft2[,3] %in% c(1,4),]
            dft2 <- dft2[dft2[,4] %in% c(1,4),]
            dft2 <- dft2[complete.cases(dft2),]
            dft2<- dft2[,c(2,3,4)]
            dfts <- get_syng_score(dft2,method)
            dfts$Disease <- disease
            dfts <- dfts[ , c("Disease" , c(setdiff(colnames(dfts) , "Disease")))]
            df_syng <- rbind(df_syng,dfts)
          }
        }
      }
    }
  }
  return(df_syng)
}
