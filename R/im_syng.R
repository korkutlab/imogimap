#' Find combinatorial association of immunotherapy co-targets with tumor intrinsic features.
#' @importFrom dplyr bind_rows across mutate group_by everything
#' @importFrom magrittr %>%
#' @importFrom data.table setkey as.data.table
#' @importFrom data.table ":=" ".SD"
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom stats complete.cases median
#' @param onco_gene A character vector of gene IDs.
#' @param icp_gene An optional character vector of immune checkpoint gene/protein IDs.
#' @param data_expression A non-negative numeric matrix or data frame containing gene/protein expressions in linear scale.
#' @param data_feature An optional numeric matrix or data frame containing immune features.
#' @param add_features An optional logical indicating if EMT score, angiogenesis score and IFNG expression should be added to immune features. Default is TRUE.
#' @param add_receptor_ligand An optional logical indicating whether receptor_ligands pair should be added. Default is FALSE.
#' @param method A character string indicating which synergy score to be used. one of "max" or "independence". Default is "max".
#' @param ndatamin minimum number of samples. Synergy score calculation will be skipped for matrices with number of rows less than ndatamin
#' @param specificity An optional logical indicating if specificity analysis should be done. Default is FALSE.
#' @param N_iteration_specificity Number of iterations for random sampling for specificity p.value calculation.
#' Default is 1000.
#' @param sensitivity An optional logical indicating if a sensitivity analysis should be done. Default is FALSE.
#' @param N_iteration_sensitivity Number of iterations for random sampling for sensitivity analysis.
#' Default is 1000.
#' @keywords Synergy, immune feature, immune checkpoint, bootstrapping
#' @return a data.frames of synergy scores and bootstrapping p.values
#' @details
#'
#' For details of synergy score calculations see get_syng_score function.
#'
#' By default (if no icp_gene is specified), icp_gene_list will be used.
#'
#' data_expression is formatted with genes/proteins as rows and samples/patients as columns.
#' For data_expression sample formats see sample_mRNA_data.
#'
#' data_feature is formatted with samples/patients as rows and immune features as columns.
#' For data_feature sample format see sample_immune_cell_fraction_data.
#'
#' For synergy score calculations all features are normalized to be on \code{[0,1]} range. For details of synergy score and significance pvalue calculations see \code{find_a_synergy} function.
#'
#' A specificity p.value is computed using random sampling with replacement from two null models, generated from one of the two genes against a set of genes randomly selected from the genome. Two P-values are calculated for the synergistic interaction of the pair against the two null models. The highest of the two P-values is used to assess the specificity of the interaction against the whole genome. The number of randomly selected genes in each null model is determined by N_iteration_specificity.
#'
#' Sensitivity (Robustness) score defined as normalized root mean square deviation of scores calculated over 70% of samples, selected via random sampling. The number of sub-sample iterations is determined by N_iteration_sensitivity.
#'
#'
#' @examples im_syng(onco_gene =  c("BRAF"),
#'                   icp_gene= c("CD274","CTLA4"),
#'                   data_expression =  sample_mRNA_data,
#'                   data_feature = sample_Leukocyte_fraction_data)
#' @export


im_syng<-function(onco_gene,icp_gene,data_expression,data_feature, ndatamin=8, add_features,method,specificity, N_iteration_specificity, sensitivity, N_iteration_sensitivity,add_receptor_ligand=FALSE){

  agent1 <- agent2 <- Immune_feature <- agent1_expression <- agent2_expression <- NULL
  specificity_pvalue <- sensitivity_R <- i.sensitivity_R <- NULL
  df_syng <- data.frame(agent1=character(),
                        agent2=character(),
                        Immune_feature=character(),
                        Synergy_score=numeric(),
                        agent1_expression=character(),
                        agent2_expression=character(),
                        wilcox_pvalue=numeric(),
                        specificity_pvalue=numeric(),
                        sensitivity_R=numeric())
  #Check inputs
  if(missing(method)){
    method <- "max"
  }else{
    method <- tolower(method)
    if(!(method=="max" || method=="independence")){
      stop("ERROR: Method is not found. Please choose a method from: max or independence.")
    }
  }

  #Add interacting genes-----------
  if(add_receptor_ligand==TRUE){
    lgn1 <- lgn_receptor_ligand[lgn_receptor_ligand$Gene1 %in% onco_gene,]$Gene2
    lgn2 <- lgn_receptor_ligand[lgn_receptor_ligand$Gene2 %in% onco_gene,]$Gene1
    lgn <- unique(c(lgn1,lgn2))
    if(length(lgn)>0){
      onco_gene<- unique(c(onco_gene,lgn))
    }
  }
  if(missing(icp_gene)){
    icp_gene <- icp_gene_list
  }
  if(add_receptor_ligand==TRUE){
    lgn1 <- lgn_receptor_ligand[lgn_receptor_ligand$Gene1 %in% icp_gene,]$Gene2
    lgn2 <- lgn_receptor_ligand[lgn_receptor_ligand$Gene2 %in% icp_gene,]$Gene1
    lgn <- unique(c(lgn1,lgn2))
    if(length(lgn)>0){
      icp_gene<- unique(c(icp_gene,lgn))
    }
  }

  #Check for co-target expressions-----------------------
  #onco_gene[onco_gene == "VSIR"]<- "C10orf54"
  #onco_gene[onco_gene == "NCR3LG1"]<- "DKFZp686O24166"
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
  #icp_gene[icp_gene=="VSIR"]<- "C10orf54"
  #icp_gene[icp_gene=="NCR3LG1"]<- "DKFZp686O24166"

  missing_icp <- icp_gene[-which(icp_gene %in% rownames(data_expression))]
  if(length(missing_icp)>0){
    warning("Missing immune checkpoints:   ",lapply(missing_icp, function(x)paste0(x,"  ")))
  }
  if(length(icp_gene)==1){
    df_icp <- as.data.frame(t(data_expression[rownames(data_expression) %in% icp_gene,]))
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
        warning("No feature data is specified. Optional features are being added.. \n")
        add_features <- TRUE
      }
    }
  }else{
    df_min <- min(data_feature,na.rm=TRUE)
    df_max <- max(data_feature,na.rm=TRUE)
    if(df_min < 0 | df_max > 1 ){
      stop("ERROR: feature is out of range. Normalize data_feature to [0,1].")
    }
    data_feature <- as.data.frame(data_feature)
    if(missing(add_features)){
      add_features <- TRUE
    }
  }
  data_feature$Tumor_Sample_ID <- rownames(data_feature)
  data_feature <-data_feature[,c("Tumor_Sample_ID",setdiff(colnames(data_feature), "Tumor_Sample_ID")),drop=F]
  data_feature0 <- data_feature

  #Calculate additional immune features---------------------
  if(add_features==TRUE){
    message("Calculating additional features...")
    df_tmp <- get_features(data_expression)
    df_tmp <- df_tmp[,colSums(is.na(df_tmp))<nrow(df_tmp)]
    df_tmp$Leukocyte_fraction<-NULL
    df_tmp$TMB_Non.silent_per_Mb<-NULL
    df_tmp$TMB_Silent_per_Mb<-NULL
    data_feature <- merge(df_tmp,data_feature,by="Tumor_Sample_ID")
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
  N_genes <- as.numeric(length(all_genes))
  N_feature <- as.numeric(ncol(data_feature)-1)

  #Calculate synergy scores -------------------------------------
  if(N_feature>0){
    message("Calculating synergy scores for immune features ...")

    pb <- txtProgressBar(min = 0, max = N_perm, char="-",style = 3)

    for(pair_ID in 1:N_perm){
      gene_ID1 <- which(colnames(df_all)==all_perms[pair_ID,1])
      gene_ID2 <- which(colnames(df_all)==all_perms[pair_ID,2])
      dft <- df_all[ , c(1 , gene_ID1,gene_ID2)]
      dft <- dft[dft[ , 2] %in% c(1 , 4) , ,drop=F]
      dft <- dft[dft[ , 3] %in% c(1 , 4) , ,drop=F ]
      if(nrow(dft) > 4){
        df_helper <-  data.frame(agent1=character(),
                                 agent2=character(),
                                 Immune_feature=character(),
                                 Synergy_score=numeric(),
                                 agent1_expression=character(),
                                 agent2_expression=character(),
                                 wilcox_pvalue=numeric())

        for(im_ID in 2:ncol(data_feature)){

          dft2 <- merge(data_feature[ , c(1 , im_ID)] , dft , by="Tumor_Sample_ID")
          dft2 <- dft2[complete.cases(dft2),]
          dft2 <- dft2[,c(2,3,4)]
          dfts <- find_a_synergy(dft2,method = method,ndatamin = ndatamin)
          df_helper <- dplyr::bind_rows(df_helper , dfts)
        }
      }else{
        df_helper <-  data.frame(agent1=colnames(dft)[1],
                                 agent2=colnames(dft)[2],
                                 Immune_feature=colnames(data_feature),
                                 Synergy_score=NA,
                                 agent1_expression=NA,
                                 agent2_expression=NA,
                                 wilcox_pvalue=NA)
      }
      if(nrow(df_helper)>0){
        df_helper$specificity_pvalue <- NA
        df_helper$sensitivity_R <- NA
        df_syng <- dplyr::bind_rows(df_syng , df_helper)
      }
      setTxtProgressBar(pb, pair_ID)
    }
    df_syng <- as.data.table(df_syng)
    data.table::setkey(df_syng,agent1,agent2,Immune_feature)

    message("\nSynergy calculation completed! ")
  }else{
    message("No feature data found.\n")
  }

  #rownames(data_feature)<- data_feature$Tumor_Sample_ID
  #data_feature<- as.matrix(data_feature[,-1,drop=F])
  #data_feature <- data_feature[,colSums(is.na(data_feature))<nrow(data_feature),drop=F]


  #Check if Specificity should be calculated---------------------------

  if(missing(specificity)){
    specificity <- FALSE
  }
  if(specificity){

    message("\nStarting Specificity analysis:\n ")
    #Build unique permutations----------------------
    df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]

    if(nrow(df_syng_complete)>0){

      all_genes <- unique( c( df_syng_complete$agent1,df_syng_complete$agent2))
      N_genes <- as.numeric( length( all_genes))
      sub_feature_list <- unique(df_syng_complete[
        df_syng_complete$Immune_feature %in% colnames( data_feature),]$Immune_feature)
      N_feature <- as.numeric( length( sub_feature_list))
      if(missing(N_iteration_specificity)){
        N_iteration_specificity <- 1000
      }
      message("  Building bootstrapping distribution from ", N_iteration_specificity ," genes")
      df_syng_complete<-as.data.table(df_syng_complete)
      #Build random bank of expression values---------
      df_bank <- data.frame()
      while(ncol(df_bank)==0){
        df_bank <- as.data.frame(t(data_expression[
          c(sample(nrow(data_expression),N_iteration_specificity ,replace = T),
            which(rownames(data_expression) %in% all_genes)),]))
        df_bank <-  df_bank[, colSums(abs(df_bank)!= 0 ) > 0]
        df_bank <- df_bank[!duplicated(as.list(df_bank))]
      }
      df_bank <- df_bank[,c(all_genes, setdiff(colnames(df_bank),all_genes))]

      #Construct quantile ranking matrices for bank---
      df_bank <- scale(log2(df_bank+1),center = T,scale = T)
      df_bank_qr <- get_quantile_rank(df_bank)
      df_bank_qr <- as.matrix(df_bank_qr[,-1])

      rownames(data_feature)<-data_feature$Tumor_Sample_ID
      #Calculate p.values for each feature----
      message(" Calculating specificity pvalues...\n")
      if(N_feature > 0){
        for( i in 1:N_feature){

          #Select feature
          select_feature <- sub_feature_list[i]
          mark_feature <- as.integer( which( colnames( data_feature) == select_feature))
          df_feature <- as.matrix(data_feature[ , mark_feature, drop=F ])
          message( "          ...",select_feature)

          #Make a list of all genes with non-NA synergy score value
          df_syng_t <- df_syng_complete[ df_syng_complete$Immune_feature == select_feature,]
          sub_genes1 <- unique( df_syng_t[ ,.( agent1 , agent1_expression)])
          sub_genes2 <- unique( df_syng_t[ ,.( agent2 , agent2_expression)])
          colnames(sub_genes1)<-c("gene","effect")
          colnames(sub_genes2)<-c("gene","effect")
          sub_genes <- unique( rbind(sub_genes1,sub_genes2))

          syng_dist <- vector( "list", nrow( sub_genes ))
          names(syng_dist) <- paste0(sub_genes$gene,"_",sub_genes$effect)

          df_bank_sub1 <- cbind(df_feature,
                                df_bank_qr[match(rownames(df_feature),rownames(df_bank_qr)),])
          df_bank_sub1 <- df_bank_sub1[complete.cases(df_bank_sub1),]
          bank_size <- as.numeric(ncol(df_bank_sub1))
          bank_mark <- N_genes+2
          df_bank_sub1_cols <- colnames( df_bank_sub1)

          #Build bootstrapping distribution
          tmpn <- nrow(sub_genes)
          pb <- txtProgressBar(min = 0, max = tmpn, char="-",style = 3)

          for( j in 1 :tmpn){
            gene_mark <- as.numeric( which(  df_bank_sub1_cols  == sub_genes$gene[ j]))
            df_bank_sub2 <- df_bank_sub1[ , c(1, gene_mark, bank_mark : bank_size)]
            df_bank_sub2 <- df_bank_sub2[ df_bank_sub2[ , 2 ] %in% c( 1, 4), ]
            bank_size2 <- as.numeric( ncol( df_bank_sub2))
            gene_effect <- sub_genes$effect[j]

            syng_dist_t <- vector()
            for(k in 3:bank_size2){
              dft <- df_bank_sub2[,c(1,2,k)]
              dft <- dft[dft[,3] %in% c(1,4),,drop=F]
              if(nrow(dft)>0){
                syng_dist_t[k-2] <- find_a_synergy(fdata = dft,
                                                   method = method,ndatamin = ndatamin,
                                                   oncogene1 = gene_effect)$Synergy_score
              }
            }
            syng_dist_t <- syng_dist_t[complete.cases(syng_dist_t)]
            syng_dist[[j]] <- syng_dist_t
            setTxtProgressBar(pb, j)
          }

          #Calculate p.values for all pairs

          for(pair_ID in 1:nrow(df_syng_t)){

            tmp_row <- df_syng_t[pair_ID,]
            gene1 <- tmp_row$agent1[1]
            gene2 <- tmp_row$agent2[1]
            effect1 <- tmp_row$agent1_expression[1]
            effect2 <- tmp_row$agent2_expression[1]
            myscore <- tmp_row$Synergy_score[1]
            mysign <- sign(myscore)
            if( mysign==0) mysign <- 1

            j1 <- which( names( syng_dist) == paste0(gene1,"_",effect1))
            tmp_dist <- syng_dist[[j1]]
            Pvalue1 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
            j2 <- which(names(syng_dist) == paste0(gene2,"_",effect2))
            tmp_dist <- syng_dist[[j2]]
            Pvalue2 <- sum(mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
            df_syng <- df_syng[.(gene1,gene2,select_feature,effect1,effect2),
                               specificity_pvalue :=max(Pvalue1, Pvalue2)]
          }
        }
      }
    }
  }

  #Check if robustness analysis should be done---------------------------

  if(missing(sensitivity)){
    sensitivity <- FALSE
  }
  if(sensitivity){
    if(missing(N_iteration_sensitivity)){
      N_iteration_sensitivity <- 1000
    }
    message("\nCalculating robustness R. Iterating ", N_iteration_sensitivity, " times:")

    #Define parameters
    df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
    if( nrow(df_syng_complete) > 0 ){
      df_syng_complete$sum <- 0
      df_syng_complete$sum2 <- 0
      df_syng_complete$N <- 0

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


      df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
      df_syng_complete1 <- df_syng_complete[df_syng_complete$Immune_feature %in%
                                              colnames(data_feature),]
      N_syng_complete1 <- nrow(df_syng_complete1)


      pb <- txtProgressBar(min = 0, max = N_iteration_sensitivity, char="-",style = 3)

      for(n_sns in 1:N_iteration_sensitivity){

        df_sub <- df_comb[sample(N_sample,N_sub,replace = F),]
        df_sub_qr <- get_quantile_rank(df_sub)
        df_sub_qr <- as.matrix(df_sub_qr[,-1])

        dft <- data_expression[,colnames(data_expression) %in% rownames(df_sub)]
        tmp_feature <- data_feature0[rownames(data_feature0) %in% rownames(df_sub),]
        if(add_features==TRUE){
          df_tmp <- get_features(dft)
          df_tmp <- df_tmp[,colSums(is.na(df_tmp))<nrow(df_tmp)]
          df_tmp$Leukocyte_fraction<-NULL
          df_tmp$TMB_Non.silent_per_Mb<-NULL
          df_tmp$TMB_Silent_per_Mb<-NULL
          tmp_feature <- merge(df_tmp,tmp_feature,by="Tumor_Sample_ID")
          rownames(tmp_feature) <- tmp_feature$Tumor_Sample_ID
        }
        #rownames(tmp_feature) <- tmp_feature$Tumor_Sample_ID
        #tmp_feature <- as.matrix(tmp_feature[,-1])

        for(pair_ID in 1:N_syng_complete1 ){
          tmp_df <- df_syng_complete1[pair_ID,]
          gene_ID1 <- tmp_df$agent1
          gene_ID2 <- tmp_df$agent2
          base_score <- tmp_df$Synergy_score
          effect1 <- tmp_df$agent1_expression
          effect2 <- tmp_df$agent2_expression
          feature_ID <- tmp_df$Immune_feature
          mark_feature <- as.integer( which( colnames( tmp_feature) == feature_ID))
          df_feature <- as.matrix(tmp_feature[ , mark_feature,drop=F])

          dft <- df_sub_qr[ , c(which(colnames(df_sub_qr)==gene_ID1),
                                which(colnames(df_sub_qr)==gene_ID2))]
          dft <- dft[dft[ , 1] %in% c(1 , 4) ,,drop=F ]
          dft <- dft[dft[ , 2] %in% c(1 , 4) ,,drop=F ]

          dft <- cbind(df_feature[match(rownames(dft),rownames(df_feature)),,drop=F], dft)
          dft <- dft[complete.cases(dft),,drop=F]

          if(nrow(dft)>0){
            dfts <- find_a_synergy(fdata = dft,
                                   method = method,ndatamin = ndatamin,
                                   oncogene1 = effect1,
                                   oncogene2 = effect2)$Synergy_score
          }else{
            dfts <- NA
          }

          df_syng_complete1$sum[pair_ID] <- sum(df_syng_complete1$sum[pair_ID],
                                                dfts, na.rm = T)
          dfts <- dfts -base_score
          df_syng_complete1$sum2[pair_ID] <- sum(df_syng_complete1$sum2[pair_ID],
                                                 dfts*dfts, na.rm = T)
          df_syng_complete1$N[pair_ID] <- sum(df_syng_complete1$N[pair_ID],
                                              abs(sign(dfts)), na.rm = T)
        }
        setTxtProgressBar(pb, n_sns)
      }

        df_syng_complete <- df_syng_complete1

      df_syng_complete$sensitivity_R <- sqrt(df_syng_complete$sum2/df_syng_complete$N) / abs(df_syng_complete$Synergy_score)
      df_syng_complete$sum <- NULL
      df_syng_complete$sum2 <- NULL
      df_syng_complete$N <- NULL
      rownames(df_syng_complete) <- NULL
      df_syng<- df_syng[df_syng_complete, sensitivity_R := i.sensitivity_R]
    }
  }
  df_syng <- as.data.frame(df_syng)
  message("\n Done\n")
df_syng[df_syng=="DKFZp686O24166"]<-"NCR3LG1"
df_syng[df_syng=="C10orf54"]<-"VSIR"
colnames(df_syng)[1:6] <- c("Gene1","Gene2",
                                "IAP","Synergy_score",
                                "Gene1_expression","Gene2_expression")

df_syng <-merge(df_syng,lgn_receptor_ligand,by=c("Gene1","Gene2"),all.x=T)
df_syng <-merge(df_syng,lgn_receptor_ligand,by.x=c("Gene1","Gene2"),by.y=c("Gene2","Gene1"),all.x=T)
df_syng$ligand_receptor_interaction <- FALSE
df_syng$ligand_receptor_interaction[df_syng$ligand_receptor_interaction.x==T] <- TRUE
df_syng$ligand_receptor_interaction[df_syng$ligand_receptor_interaction.y==T] <- TRUE
df_syng$ligand_receptor_interaction.x <- NULL
df_syng$ligand_receptor_interaction.y <- NULL

return(df_syng)
}

