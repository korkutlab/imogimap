#' Find combinatorial association of immunotherapy co-targets with tumor intrinsic features.
#'
#' @param cotarget A charachter vector of gene IDs.
#' @param checkpoint An optional charachter vector of immune checkpoint gene/protein IDs.
#' @param data_expression A non-negative numeric matrix or data frame containing gene/protein expressions in linear scale.
#' @param data_feature An optional numeric matrix or data frame containing immune features.
#' @param add_features An optional logical indicating if EMT score, angiogenesis score and IFNG expression should be added to immune features. Default is TRUE.
#' @param add_pvalue An optional logical indicating if a prandom bootstrapping value should be calculated. Default is FALSE.
#' @param N_iteration Number of iterrations for random bootstrapping
#' @keywords Synergy, immune feature, immune checkpoint, bootstrapping
#' @return a dataframes of synergy scores and bootstrapping pvalues
#' @details
#'
#' For details of synergy score calculations see get_syng_score function.
#'
#' By default (if no checkpoint is specified), icp_gene_list will be used.
#'
#' data_expression is formatted with genes/proteins as rows and samples/patients as columns.
#' For data_expression sample formats see sample_mRNA_data.
#'
#' data_feature is formated with samples/patients as rows and immune features as columns.
#' For data_feature sample format see sample_immune_cell_fraction_data.
#'
#' A p.value is computed using random bootstraping with replacement from the distibution of synergy scores for each immune checkpointand-immune feature pair. The default values of N_iteration is 1000.
#'
#' @examples im_syng(cotarget =  c("BRAF"),
#'                   checkpoint= c("CD274","CTLA4"),
#'                   data_expression =  sample_mRNA_data,
#'                   data_feature = sample_Leukocyte_fraction_data,
#'                   add_features=T,add_pvalue=T,
#'                   N_iteration=1000)
#' @export


im_syng<-function(cotarget,checkpoint,data_expression,data_feature, add_features,add_pvalue,N_iteration){


  #Check for co-target expressions-----------------------
  data_expression <- as.data.frame(data_expression)
  if(length(cotarget)==1){
    df_selected <- as.data.frame(t(data_expression[rownames(data_expression) %in% cotarget,]))
    colnames(df_selected) <- cotarget
  }else{
    df_selected <- t(data_expression[rownames(data_expression) %in% cotarget,])
  }
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a co-target gene from expression data")
  }

  #Check for immune checkpoint expressions----------------
  if(missing(checkpoint)){
    checkpoint <- icp_gene_list
  }
  missing_icp <- checkpoint[-which(checkpoint %in% rownames(data_expression))]
  if(length(missing_icp)>0){
    warning("Missing immune checkpoints:   ",lapply(missing_icp, function(x)paste0(x,"  ")))
  }
  if(length(checkpoint)==1){
    df_icp <- as.data.frame(data_expression[rownames(data_expression) %in% checkpoint,])
    colnames(df_icp) <- checkpoint
  }else{
    df_icp <- t(data_expression[rownames(data_expression) %in% checkpoint,])
  }
  if(nrow(df_icp)==0){
    stop("No immune checkpoint found in expression data. expression data")
  }

  #Check for immune features------------------------------
  if(missing(data_feature)){
    data_feature <- data.frame(row.names = colnames(data_expression))
    if(missing(add_features)){
      add_features <- TRUE
    }else{
      if(add_features==FALSE){
        warning("No feature data is specified. Using optional features. \n")
        add_features <- TRUE
      }
    }
  }else{
    data_feature <- as.data.frame(data_feature)
    if(missing(add_features)){
      add_features <- TRUE
    }
  }
  data_feature$Tumor_Sample_ID <- rownames(data_feature)

  #Check for additional immune features---------------------
  if(add_features==T){
    df_EMT <-  get_emt_score(data_expression)
    if(nrow(df_EMT)==0){
      warning("No EMT signature marker found.\n")
    }else{
      data_feature <- merge(df_EMT,data_feature,by="Tumor_Sample_ID")
    }

    df_ang <-  get_angio_score(data_expression)
    if(nrow(df_ang)==0){
      warning("No angiogenesis signature marker found.\n")
    }else{
      data_feature <- merge(df_ang,data_feature,by="Tumor_Sample_ID")
    }

    df_ifng <-  as.data.frame(t(data_expression[rownames(data_expression)=="IFNG",]))
    if(nrow(df_ang)==0){
      warning("No IFNG expression found.\n")
    }else{
      df_ifng$Tumor_Sample_ID <- rownames(df_ifng)
      data_feature <- merge(df_ifng,data_feature,by="Tumor_Sample_ID")
    }
  }


  #Construct quantile ranking matrices--------------
  df_selected <- scale(log2(df_selected+1),center = T,scale = T)
  df_icp <- scale(log2(df_icp+1),center = T,scale = T)
  df_select_qr <- get_quantile_rank(df_selected)
  df_icp_qr <- get_quantile_rank(df_icp)



  #Build random bootstrapping distribution----------
  if(missing(add_pvalue)){
    add_pvalue <- FALSE
  }
  if(add_pvalue){
    if(missing(N_iteration)){
      N_iteration<-1000
    }
    df_syng <- data.frame(Co_target=character(),
      Immune_checkpoint=character(),
      Immune_feature=character(),
      Synergy_score=numeric(),
      pvalue=numeric())

    #Randomly Select N_iteration rows from expression data
    df_bank <- as.data.frame(t(data_expression[
      c(sample(nrow(data_expression),N_iteration,replace = T),
        which(rownames(data_expression) %in% checkpoint)),]))
    df_bank <- df_bank[, -which(colSums(abs(df_bank)) == 0 )]
    #Construct quantile ranking matrices
    df_bank <- scale(log2(df_bank+1),center = T,scale = T)
    df_bank_qr <- get_quantile_rank(df_bank)

    #Merge with feature data
    df_bank_merged <- merge(data_feature,df_bank_qr,by="Tumor_Sample_ID")
    df_bank_merged <- as.data.frame(df_bank_merged)

    N_checkpoint <- as.numeric(length(checkpoint))
    N_feature <- as.numeric(ncol(data_feature))-1

    #Loop over checkpoints
    for( i in 1:N_checkpoint){
      checkpoint_mark <- as.numeric(ncol(df_bank_merged))-N_checkpoint+i
      df_bank_sub1 <- df_bank_merged[,c(1:(checkpoint_mark -i),checkpoint_mark )]
      df_bank_sub1 <- df_bank_sub1[df_bank_sub1[,ncol(df_bank_sub1) ] %in% c(1,4),]

      #Loop over features
      for( j in 1:N_feature){
        df_bank_sub2 <- df_bank_sub1[,c(1,(j+1),(N_feature+2):ncol(df_bank_sub1))]
        bank_size <- as.numeric(ncol(df_bank_sub2))

        #Build distribution
        syng_dist <- vector()
        for(k in 3:(bank_size-1)){
          dft <- df_bank_sub2[,c(1,2,k,bank_size)]
          dft <- dft[dft[,3] %in% c(1,4),]
          dft <- dft[complete.cases(dft),]
          dft<- dft[,c(2,3,4)]
          syng_dist[k-2] <- get_syng_score(dft)$Synergy_score
        }
        syng_dist <- syng_dist[complete.cases(syng_dist)]

        #Calculate synergy scores and probaility values for co-targets------------------
        for(cotarget_ID in 2:ncol(df_select_qr)){
          dft <- merge(df_select_qr[,c(1,cotarget_ID)],df_icp_qr[,c(1,i+1)],by="Tumor_Sample_ID")
          dft <- merge(data_feature[,c(1,j+1)],dft,by="Tumor_Sample_ID")
          dft <- dft[dft[,3] %in% c(1,4),]
          dft <- dft[dft[,4] %in% c(1,4),]
          dft<- dft[,c(2,3,4)]
          dft <- dft[complete.cases(dft),]
          cc <- get_syng_score(dft)
          myscore <- cc$Synergy_score
          if(is.na(myscore)){
            cc$Pvalue <- NA
          }else{
            mysign <- sign(myscore)
            if(mysign==0) mysign<-1
            cc$Pvalue <- sum(mysign*syng_dist > mysign*myscore)/length(syng_dist)
          }
          df_syng <- rbind(df_syng,cc)
        }
      }
    }
  }else{

    df_syng <- data.frame(Co_target=character(),
      Immune_checkpoint=character(),
      Immune_feature=character(),
      Synergy_score=numeric())
    #Calculate synergy scores -------------------------------------
    for(cotarget_ID in 2:ncol(df_select_qr)){
      for(checkpoint_ID in 2:ncol(df_icp_qr)){
        dft <- merge(df_select_qr[,c(1,cotarget_ID)],df_icp_qr[,c(1,checkpoint_ID)],by="Tumor_Sample_ID")
        for(feature_ID in 2:ncol(data_feature)){
          dft2 <- merge(data_feature[,c(1,feature_ID)],dft,by="Tumor_Sample_ID")
          dft2 <- dft2[dft2[,3] %in% c(1,4),]
          dft2 <- dft2[dft2[,4] %in% c(1,4),]
          dft2<- dft2[,c(2,3,4)]
          dft2 <- dft2[complete.cases(dft2),]
          df_syng <- rbind(df_syng,get_syng_score(dft2))
        }
      }
    }

  }

  return(df_syng)
}


