#' Find combinatorial association of immunotherapy co-targets with tumor intrinsic features.
#'
#' @param onco_gene A charachter vector of gene IDs.
#' @param icp_gene An optional charachter vector of immune checkpoint gene/protein IDs.
#' @param data_expression A non-negative numeric matrix or data frame containing gene/protein expressions in linear scale.
#' @param data_feature An optional numeric matrix or data frame containing immune features.
#' @param add_features An optional logical indicating if EMT score, angiogenesis score and IFNG expression should be added to immune features. Default is TRUE.
#' @param method a charachter string indicating which synergy score to be used. one of "max" or "independence". Default is "max".
#' @param add_pvalue An optional logical indicating if a prandom bootstrapping value should be calculated. Default is FALSE.
#' @param N_iteration Number of iterrations for random bootstrapping
#' @keywords Synergy, immune feature, immune checkpoint, bootstrapping
#' @return a dataframes of synergy scores and bootstrapping pvalues
#' @details
#'
#' For details of synergy score calculations see get_syng_score function.
#'
#' By default (if no icp_gene is specified), icp_gene_list will be used.
#'
#' data_expression is formatted with genes/proteins as rows and samples/patients as columns.
#' For data_expression sample formats see sample_mRNA_data.
#'
#' data_feature is formated with samples/patients as rows and immune features as columns.
#' For data_feature sample format see sample_immune_cell_fraction_data.
#'
#' A p.value is computed using random bootstraping with replacement from the distibution of synergy scores for each immune checkpoint and-immune feature pair. The default values of N_iteration is 1000.
#'
#' @examples im_syng(onco_gene =  c("BRAF"),
#'                   icp_gene= c("CD274","CTLA4"),
#'                   data_expression =  sample_mRNA_data,
#'                   data_feature = sample_Leukocyte_fraction_data,
#'                   add_features=T,add_pvalue=T,
#'                   N_iteration=1000)
#' @export


im_syng<-function(onco_gene,icp_gene,data_expression,data_feature, add_features,method,add_pvalue,N_iteration){

  #Check inputs
  if(missing(method)){
    method <- "max"
  }else{
    method <- tolower(method)
    if(!(method=="max" || method=="independence")){
      stop("ERROR: Method is not found. Please choose a method from: max or independence.")
    }
  }
  #Check for co-target expressions-----------------------
  data_expression <- as.data.frame(data_expression)
  if(length(onco_gene)==1){
    df_selected <- as.data.frame(t(data_expression[rownames(data_expression) %in% onco_gene,]))
    colnames(df_selected) <- onco_gene
  }else{
    df_selected <- t(data_expression[rownames(data_expression) %in% onco_gene,])
  }
  onco_gene_sub <- colnames(df_selected)
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a co-target gene from expression data")
  }

  #Check for immune checkpoint expressions----------------
  if(missing(icp_gene)){
    icp_gene <- icp_gene_list
  }
  missing_icp <- icp_gene[-which(icp_gene %in% rownames(data_expression))]
  if(length(missing_icp)>0){
    warning("Missing immune checkpoints:   ",lapply(missing_icp, function(x)paste0(x,"  ")))
  }
  if(length(icp_gene)==1){
    df_icp <- as.data.frame(data_expression[rownames(data_expression) %in% icp_gene,])
    colnames(df_icp) <- icp_gene
  }else{
    df_icp <- t(data_expression[rownames(data_expression) %in% icp_gene,])
  }
  icp_gene_sub <- colnames(df_icp)
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
      sd_EMT <- sd(df_EMT$EMTscore,na.rm = T)
      df_EMT$EMTscore <- ( tanh( df_EMT$EMTscore/sd_EMT ) + 1 ) / 2
      data_feature <- merge(df_EMT,data_feature,by="Tumor_Sample_ID")
    }

    df_ang <-  get_angio_score(data_expression)
    if(nrow(df_ang)==0){
      warning("No angiogenesis signature marker found.\n")
    }else{
      sd_AG <- sd(df_ang$AGscore,na.rm = T)
      df_ang$AGscore <- ( tanh( df_ang$AGscore/sd_AG ) + 1 ) / 2
      data_feature <- merge(df_ang,data_feature,by="Tumor_Sample_ID")
    }

    df_ifng <-get_ifng_score(data_expression)
    if(nrow(df_ang)==0){
      warning("No IFNG expression found.\n")
    }else{
      sd_IFNG <- sd(df_ifng$IFNGscore,na.rm = T)
      df_ifng$IFNGscore <- ( tanh( df_ifng$IFNGscore/sd_IFNG ) + 1 ) / 2
      data_feature <- merge(df_ifng,data_feature,by="Tumor_Sample_ID")
    }
  }


  #Construct quantile ranking matrices--------------
  df_selected <- scale(log2(df_selected+1),center = T,scale = T)
  df_icp <- scale(log2(df_icp+1),center = T,scale = T)
  df_select_qr <- get_quantile_rank(df_selected)
  df_icp_qr <- get_quantile_rank(df_icp)
  df_all <- merge(df_select_qr,df_icp_qr)

  #Build unique permutations
  all_genes <- unique(c(onco_gene_sub,icp_gene_sub))
  all_perms <- t(combn(all_genes,m = 2))
  N_perm <- as.numeric(nrow(all_perms))

  #Build random bootstrapping distribution----------
  if(missing(add_pvalue)){
    add_pvalue <- FALSE
  }
  if(add_pvalue){
    if(missing(N_iteration)){
      N_iteration<-1000
    }

    message("Building bootstrapping distribution ...")

    #Randomly Select N_iteration rows from expression data
    df_bank <- as.data.frame(t(data_expression[
      c(sample(nrow(data_expression),N_iteration,replace = T),
        which(rownames(data_expression) %in% all_genes)),]))
    df_bank <- df_bank[, colSums(abs(df_bank)) != 0 ]
    #Construct quantile ranking matrices
    df_bank <- scale(log2(df_bank+1),center = T,scale = T)
    df_bank_qr <- get_quantile_rank(df_bank)

    #Merge with feature data
    df_bank_merged <- merge(data_feature,df_bank_qr,by="Tumor_Sample_ID")
    df_bank_merged <- as.data.frame(df_bank_merged)

    #Define parameters
    df_syng <- data.frame(agent1=character(),
      agent2=character(),
      Immune_feature=character(),
      Synergy_score=numeric(),
      pvalue=numeric())

    N_genes <- as.numeric(length(all_genes))
    N_feature <- as.numeric(ncol(data_feature))-1

    message("Calculating scynergy scores for feature:")

    #Get synergy score and p.value for each feature -----------
    for( i in 1:N_feature){
      message(colnames(df_bank_merged)[i+1]," ...")
      syng_dist <- vector("list",N_genes)
      names(syng_dist) <- all_genes
      df_bank_sub1 <- df_bank_merged[,c(1,(i+1),(N_feature+2):ncol(df_bank_merged))]

      #Build bootstrapping distribution ----------
      for( j in 1:N_genes){
        gene_mark <- as.numeric(ncol(df_bank_sub1))-N_genes+j
        df_bank_sub2 <- df_bank_sub1[,c(1:(gene_mark -j),gene_mark )]
        df_bank_sub2 <- df_bank_sub2[df_bank_sub2[,ncol(df_bank_sub2) ] %in% c(1,4),]
        bank_size <- as.numeric(ncol(df_bank_sub2))
        syng_dist_t <- vector()
        for(k in 3:(bank_size-1)){
          dft <- df_bank_sub2[,c(1,2,k,bank_size)]
          dft <- dft[dft[,3] %in% c(1,4),]
          dft <- dft[complete.cases(dft),]
          dft<- dft[,c(2,3,4)]
          syng_dist_t[k-2] <- get_syng_score(dft,method)$Synergy_score
        }
        syng_dist_t <- syng_dist_t[complete.cases(syng_dist_t)]
        syng_dist[[j]] <- syng_dist_t
      }
      # Loop over all pairs--------
      for(pair_ID in 1:N_perm){
        gene1 <- all_perms[pair_ID,1]
        gene2 <- all_perms[pair_ID,2]
        gene_ID1 <- which(colnames(df_all)==gene1)
        gene_ID2 <- which(colnames(df_all)==gene2)
        dft <- df_all[ , c(1 , gene_ID1,gene_ID2)]
        dft <- dft[dft[ , 2] %in% c(1 , 4) , ]
        dft <- dft[dft[ , 3] %in% c(1 , 4) , ]
        dft <- merge(data_feature[ , c(1 , (i+1))] , dft , by="Tumor_Sample_ID")
        dft <- dft[complete.cases(dft),]
        dft <- dft[,c(2,3,4)]
        #Get synergy score--------
        cc <- get_syng_score(dft,method)
        myscore <- cc$Synergy_score
        if(is.na(myscore)){
          cc$Pvalue1 <- NA
          cc$Pvalue2 <- NA
        }else{
          mysign <- sign(myscore)
          if(mysign==0) mysign<-1
          #Get pvalue---------
          j1 <- which(names(syng_dist)==gene1)
          tmp_dist <- syng_dist[[j1]]
          cc$Pvalue1 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
          j2 <- which(names(syng_dist)==gene2)
          tmp_dist <- syng_dist[[j2]]
          cc$Pvalue2 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
        }
        cc$Disease <- disease
        cc <- cc[ , c("Disease" , c(setdiff(colnames(cc) , "Disease")))]
        df_syng <- rbind(df_syng,cc)
      }
      message(cc$Immune_feature," ...done!\n")
    }

  }else{

    df_syng <- data.frame(agent1=character(),
      agent2=character(),
      Immune_feature=character(),
      Synergy_score=numeric())

    #Calculate synergy scores -------------------------------------

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
        df_syng <- rbind(df_syng , dfts)

      }
    }
  }

  return(df_syng)
}


