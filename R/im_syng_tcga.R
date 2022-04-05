#'  Calculation and statistical assessment of synergistic associations using TCGA data.
#'
#' @param onco_gene A character vector of gene Hugo symbols.
#' @param icp_gene An optional character vector of immune checkpoint gene Hugo symbols.
#' @param cohort A list of TCGA diseases
#' @param select_iap  An optional  character vector or numeric matrix or data.frame.
#' @param method A character string indicating which synergy score to be used. one of "max" or "independence".
#' Default is "max".
#'@param ndatamin minimum number of samples. Synergy score calculation will be skipped for matrices with number of rows less than ndatamin
#' @param specificity An optional logical indicating if specificity analysis should be done. Default is FALSE.
#' @param N_iteration_specificity Number of iterations for random sampling for specificity p.value calculation.
#' Default is 1000.
#' @param sensitivity An optional logical indicating if a sensitivity analysis should be done. Default is FALSE.
#' @param N_iteration_sensitivity Number of iterations for random sampling for sensitivity analysis.
#' Default is 1000.
#' @param sample_list An optional character vector of TCGA samples barcodes indicating a subset of samples within a cohort.
#'
#' @keywords Synergy scoring, immune feature, immune checkpoint, bootstrapping, TCGA
#' @return A data.frame of synergy scores and bootstrapping p.values.
#' @description
#'
#' Takes a list of Tumor-intrinsic pathway(TIP) genes and returns their combinatorial association with immune checkpoint(ICP) genes by evaluating their synergistic impact on immune-associated phenotypes(IAP) using RNASeq2GeneNorm expressions as provided by \pkg{curatedTCGAData}.
#'
#' @details
#'
#   Output is a dataframe, with each row containing two gene Hugo IDs, an IAP name, the inferred impact of each single gene on the IAP, the calculated synergy score and the statistical assessment parameters to measure robustness, significance, and specificity.
#
# IAP names are listed in TCGA_immune_features_list.
#'If no icp_gene is specified, the default icp_gene_list will be used.
#'If select_iap is a character vector, it must be any sub-list of IAP names as listed in TCGA_immune_features_list. If a numeric matrix or data.frame, each column represents a user-defined IAP and must have a range between \code{[0,1]}. If select_iap is missing all IAPs listed in TCGA_immune_features_list will be considered for analysis.
#'
#' For synergy score calculations all features are normalized to be on \code{[0,1]} range. For details of synergy score and significance pvalue calculations see \code{find_a_synergy} function.
#'
#' A specificity p.value is computed using random sampling with replacement from two null models, generated from one of the two genes against a set of genes randomly selected from the genome. Two P-values are calculated for the synergistic interaction of the pair against the two null models. The highest of the two P-values is used to assess the specificity of the interaction against the whole genome. The number of randomly selected genes in each null model is determined by N_iteration_specificity.
#'
#' Sensitivity (Robustness) score defined as normalized root mean square deviation of scores calculated over 70% of samples, selected via random sampling. The number of sub-sample iterations is determined by N_iteration_sensitivity.
#'
#'All barcodes in sample_list must be 15 character long and belong to the same cohort. When sample_list is provided, cohort should be the disease cohort that they belong to, otherwise only the first element of the cohort list will be used.
#'
#' @examples im_syng_tcga(onco_gene=c("TGFB1","SERPINB9"), cohort=c("ucec"))
#'
#' @importFrom dplyr bind_rows across mutate group_by everything distinct
#' @importFrom magrittr %>%
#' @importFrom data.table setkey as.data.table
#' @import curatedTCGAData
#' @importFrom data.table ":=" ".SD"
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom stats complete.cases median
#'
#' @export

im_syng_tcga <- function(onco_gene, icp_gene, cohort, select_iap, method, ndatamin=8,specificity, N_iteration_specificity=1000, sensitivity, N_iteration_sensitivity=1000, sample_list){
  
  PATIENT_BARCODE <- Disease <- agent1 <- agent2 <- Immune_feature <- NULL
  agent1_expression <- agent2_expression <- specificity_pvalue <- NULL
  sensitivity_R <- i.sensitivity_R <- NULL
  
  df_syng_all <- data.frame(Disease=character(),
                            agent1=character(),
                            agent2=character(),
                            Immune_feature=character(),
                            Synergy_score=numeric(),
                            agent1_expression=character(),
                            agent2_expression=character(),
                            wilcox_pvalue=numeric(),
                            specificity_pvalue=numeric(),
                            sensitivity_R=numeric())
  
  #Check input parameters------------------------
  if(missing(ndatamin)){
    ndatamin <- 8
  }
  
  if(missing(method)){
    method <- "max"
  }else{
    method <- tolower(method)
    if(!(method=="max" || method=="independence")){
      stop("ERROR: Method is not found. Please choose a method from: max or independence.")
    }
  }
  
  if(!missing(sample_list)){
    cohort <- cohort[1]
  }
  cohort <- tolower(cohort)
  
  
  #Loop over diseases-----------------------------
  for(cohortID in 1:length(cohort)){
    
    df_syng <- data.frame(Disease=character(),
                          agent1=character(),
                          agent2=character(),
                          Immune_feature=character(),
                          Synergy_score=numeric(),
                          agent1_expression=character(),
                          agent2_expression=character(),
                          wilcox_pvalue=numeric(),
                          specificity_pvalue=numeric(),
                          sensitivity_R=numeric())
    
    #Get expression data----------------------------
    
    disease <- cohort[cohortID]
    message("\nReading TCGA ", toupper(disease), " data\n")
    
    df <- curatedTCGAData::curatedTCGAData(diseaseCode = disease, version = "1.1.38",
                                           assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]
    
    data_expression <- df@assays$data@listData[[1]]
    colnames(data_expression)<-  substr(colnames(data_expression), 1, 15)
    if(!missing(sample_list)){
      data_expression<-data_expression[,sample_list ]
      if(ncol(data_expression)==0){
        stop("ERROR: barcodes not found.")
      }
    }
    
    message("Quantile ranking Gene expressions...")
    
    #Check for co-target expressions----------------
    onco_gene[onco_gene=="VSIR"]<- "C10orf54"
    onco_gene[onco_gene=="NCR3LG1"]<- "DKFZp686O24166"
    
    if(length(onco_gene)==1){
      df_selected <- as.data.frame((data_expression[rownames(data_expression ) %in% onco_gene,]))
      colnames(df_selected) <- onco_gene
    }else{
      df_selected <- t(data_expression[rownames(data_expression ) %in% onco_gene,])
      missing_genes <- onco_gene[-which(onco_gene %in% rownames(df_selected))]
      if(length(missing_genes)>0){
        warning(length(missing_genes)," Some oncogenes are missing:  \n  ",
                lapply(missing_genes, function(x)paste0(x,"  ")))
      }
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
    if(missing(icp_gene)){
      icp_gene <- icp_gene_list
    }
    icp_gene[icp_gene=="VSIR"]<- "C10orf54"
    icp_gene[icp_gene=="NCR3LG1"]<- "DKFZp686O24166"
    
    if(length(icp_gene)==1){
      df_icp <- as.data.frame(data_expression[rownames(data_expression) %in% icp_gene,])
      colnames(df_icp) <- icp_gene
    }else{
      df_icp <- t(data_expression[rownames(data_expression) %in% icp_gene,])
      missing_genes <- icp_gene[-which(icp_gene %in% rownames(df_icp))]
      if(length(missing_genes)>0){
        warning(length(missing_genes)," Some oncogenes are missing:  \n  ",
                lapply(missing_genes, function(x)paste0(x,"  ")))
      }
      
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
    
    #Construct quantile ranking matrices for each sample--------
    
    df_selected <- scale(log2(df_selected+1),center = T,scale = T)
    df_icp <- scale(log2(df_icp+1),center = T,scale = T)
    df_select_qr <- get_quantile_rank(df_selected)
    df_icp_qr <- get_quantile_rank(df_icp)
    df_all <- merge(df_select_qr,df_icp_qr)
    rownames(df_all)<- df_all$Tumor_Sample_ID
    df_all <- as.matrix(df_all[,-1])
    
    #Construct quantile ranking matrices for each patient---------
    #For immune cell count features as calculated by CIBERSORT,
    #PATIENT_BARCODE is used instead of Tumor_Sample_ID
    
    tmp <- as.data.frame(df_selected)
    tmp$PATIENT_BARCODE <- substr(rownames(tmp), 1, 12)
    tmp <- tmp %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
    df_selected2 <- as.matrix(tmp[,-which(colnames(tmp)=="PATIENT_BARCODE")])
    rownames(df_selected2)<- tmp$PATIENT_BARCODE
    df_select_qr2 <- get_quantile_rank(df_selected2)
    
    tmp <- as.data.frame(df_icp)
    tmp$PATIENT_BARCODE <- substr(rownames(tmp), 1, 12)
    tmp <- tmp %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
    df_icp2 <- as.matrix(tmp[,-which(colnames(tmp)=="PATIENT_BARCODE")])
    rownames(df_icp2)<- tmp$PATIENT_BARCODE
    df_icp_qr2 <- get_quantile_rank(df_icp2)
    df_all2 <- merge(df_select_qr2,df_icp_qr2)
    rownames(df_all2)<- df_all2$Tumor_Sample_ID
    df_all2 <- as.matrix(df_all2[,-1])
    
    #Construct IAPs-----------------------------------
    
    message("Constructing IAPs...")
    
    if(!missing(select_iap)){
      if(is.data.frame(select_iap) || is.matrix(select_iap)){
        df_min <- min(select_iap,na.rm=TRUE)
        df_max <- max(select_iap,na.rm=TRUE)
        if(df_min < 0 | df_max > 1 ){
          stop("ERROR: feature is out of range. Normalize data_feature to [0,1].")
        }
        select_iap <- as.data.frame(select_iap)
        data_feature <- select_iap
      }else {
        data_feature <- get_features(data_expression)
        data_feature <- data_feature[,which(colnames(data_feature) %in% c("Tumor_Sample_ID",select_iap)),drop=F]
        rownames(data_feature)<- data_feature$Tumor_Sample_ID
        if(ncol(data_feature)>1){
          data_feature<- as.matrix(data_feature[,-1,drop=F])
          data_feature <- data_feature[,colSums(is.na(data_feature))<nrow(data_feature),drop=F]
        }else{
          data_feature<-data.frame()
        }
        data_cell <- TCGA_IMCell_fraction
        data_cell <- data_cell[,which(colnames(data_cell) %in% c("PATIENT_BARCODE",select_iap)),drop=F]
        rownames(data_cell)<- data_cell$PATIENT_BARCODE
        if(ncol(data_cell)>1){
          data_cell <- as.matrix(data_cell[,-1,drop=F])
          data_cell <- data_cell[,colSums(is.na(data_cell))<nrow(data_cell),drop=F]
        }else{
          data_cell <-data.frame()
        }
      }
    } else {
      select_iap <- TCGA_immune_features_list
      data_feature <- get_features(data_expression)
      rownames(data_feature)<- data_feature$Tumor_Sample_ID
      data_feature<- as.matrix(data_feature[,-1,drop=F])
      data_feature <- data_feature[,colSums(is.na(data_feature))<nrow(data_feature),drop=F]
      
      data_cell <- TCGA_IMCell_fraction
      rownames(data_cell)<- data_cell$PATIENT_BARCODE
      data_cell <- as.matrix(data_cell[,-1])
      data_cell <- data_cell[,colSums(is.na(data_cell))<nrow(data_cell),drop=F]
    }
    
    #Build unique permutations
    
    all_genes <- unique(c(onco_gene_sub, icp_gene_sub))
    all_perms <- t(combn(all_genes,m = 2))
    N_perm <- as.numeric(nrow(all_perms))
    N_genes <- as.numeric(length(all_genes))
    N_feature <- as.numeric(ncol(data_feature))
    N_imcell_feature <- as.numeric(ncol(data_cell))
    
    
    #Calculate synergy scores for features ---------
    
    message("Calculating synergy scores ...\n")
    if(N_feature>0){
      pb <- txtProgressBar(min = 0, max = N_perm, char="-",style = 3)
      
      for(pair_ID in 1:N_perm){
        gene_ID1 <- which(colnames(df_all)==all_perms[pair_ID,1])
        gene_ID2 <- which(colnames(df_all)==all_perms[pair_ID,2])
        dft <- df_all[ , c(gene_ID1,gene_ID2)]
        dft <- dft[dft[ , 1] %in% c(1 , 4) ,,drop=F ]
        dft <- dft[dft[ , 2] %in% c(1 , 4) ,,drop=F ]
        if(nrow(dft) > ndatamin){
          df_helper <-  data.frame(agent1=character(),
                                   agent2=character(),
                                   Immune_feature=character(),
                                   Synergy_score=numeric(),
                                   agent1_expression=character(),
                                   agent2_expression=character(),
                                   wilcox_pvalue=numeric())
          
          for(im_ID in 1:ncol(data_feature)){
            dft2 <- cbind(data_feature[ , im_ID][match(rownames(dft),rownames(data_feature))], dft)
            colnames(dft2)[1]<- colnames(data_feature)[im_ID]
            dft2 <- dft2[complete.cases(dft2),,drop=F]
            
            dfts <- find_a_synergy(dft2,method = method,ndatamin = ndatamin)
            df_helper <- bind_rows(df_helper , dfts)
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
          df_helper$Disease <- disease
          df_helper <- df_helper[ , c("Disease" , c(setdiff(colnames(df_helper) , "Disease")))]
          df_helper$specificity_pvalue <- NA
          df_helper$sensitivity_R <- NA
          df_syng <- bind_rows(df_syng , df_helper)
        }
        setTxtProgressBar(pb, pair_ID)
      }
    }
    #Calculate synergy scores for immune cell features--------
    if(N_imcell_feature>0){
      pb <- txtProgressBar(min = 0, max = N_perm, char="-",style = 3)
      
      for(pair_ID in 1:N_perm){
        gene_ID1 <- which(colnames(df_all2)==all_perms[pair_ID,1])
        gene_ID2 <- which(colnames(df_all2)==all_perms[pair_ID,2])
        dft <- df_all2[ , c(gene_ID1,gene_ID2)]
        dft <- dft[dft[ , 1] %in% c(1 , 4) , ,drop=F]
        dft <- dft[dft[ , 2] %in% c(1 , 4) , ,drop=F]
        if(nrow(dft) > ndatamin){
          df_helper <-  data.frame(agent1=character(),
                                   agent2=character(),
                                   Immune_feature=character(),
                                   Synergy_score=numeric(),
                                   agent1_expression=character(),
                                   agent2_expression=character(),
                                   wilcox_pvalue=numeric())
          
          for(if_ID in 1:ncol(data_cell)){
            dft2 <- cbind(data_cell[ , if_ID][match(rownames(dft),rownames(data_cell))], dft)
            colnames(dft2)[1]<- colnames(data_cell)[if_ID]
            dft2 <- dft2[complete.cases(dft2),,drop=F]
            if(nrow(dft2)>0){
              dfts <- find_a_synergy(dft2,method = method,ndatamin = ndatamin)
              df_helper <- bind_rows(df_helper , dfts)
            }
          }
        }else{
          df_helper <-  data.frame(agent1=colnames(dft)[1],
                                   agent2=colnames(dft)[2],
                                   Immune_feature=colnames(data_cell),
                                   Synergy_score=NA,
                                   agent1_expression=NA,
                                   agent2_expression=NA,
                                   wilcox_pvalue=NA)
        }
        if(nrow(df_helper)>0){
          df_helper$Disease <- disease
          df_helper <- df_helper[ , c("Disease" , c(setdiff(colnames(df_helper) , "Disease")))]
          df_helper$specificity_pvalue <- NA
          df_helper$sensitivity_R <- NA
          df_syng <- bind_rows(df_syng , df_helper)
        }
        setTxtProgressBar(pb, pair_ID)
      }
    }
    df_syng <- data.table::as.data.table(df_syng)
    setkey(df_syng, Disease, agent1,agent2,Immune_feature)
    
    message("\nSynergy calculation completed!\n ")
    
    
    #Check if p.value should be calculated-----------------------
    if(missing(specificity)){
      specificity <- FALSE
    }
    if(specificity){
      message("\nStarting Specificity analysis:\n ")
      
      #Build unique permutations----------------------
      df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
      df_syng_complete <- df_syng_complete[ which(df_syng_complete$Disease==disease), ]
      
      if(nrow(df_syng_complete)>0){
        
        all_genes <- unique( c(df_syng_complete$agent1, df_syng_complete$agent2))
        N_genes <- as.numeric( length(all_genes))
        sub_feature_list <- unique(df_syng_complete[
          df_syng_complete$Immune_feature %in% colnames( data_feature),]$Immune_feature)
        N_feature <- as.numeric( length( sub_feature_list))
        sub_imcell_feature_list <- unique(df_syng_complete[
          df_syng_complete$Immune_feature %in% colnames( TCGA_IMCell_fraction ),]$Immune_feature)
        N_imcell_feature <-as.numeric( length( sub_imcell_feature_list))
        
        
        message("  Building bootstrapping distribution from ", N_iteration_specificity ," genes")
        
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
        
        tmp <- as.data.frame(df_bank)
        tmp$PATIENT_BARCODE <- substr(rownames(tmp), 1, 12)
        tmp <- data.table::as.data.table(tmp)
        tmp <- tmp[,lapply(.SD,median),by=PATIENT_BARCODE]
        tmp<- as.data.frame(tmp)
        df_bank2 <- as.matrix(tmp[,-which(colnames(tmp)=="PATIENT_BARCODE")])
        rownames(df_bank2)<- tmp$PATIENT_BARCODE
        df_bank_qr2 <- get_quantile_rank(df_bank2)
        
        df_bank_qr <- as.matrix(df_bank_qr[,-1])
        df_bank_qr2 <- as.matrix(df_bank_qr2[,-1])
        
        #Calculate p.values for each feature----
        
        message("Calculating pvalues...")
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
            df_bank_sub1 <- df_bank_sub1[complete.cases(df_bank_sub1),,drop=F]
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
              
              syng_dist_t <- rep(NA_real_, bank_size2-2)
              for(k in 3:bank_size2){
                dft <- df_bank_sub2[,c(1,2,k)]
                dft <- dft[dft[,3] %in% c(1,4),,drop=F]
                syng_dist_t[k-2] <- find_a_synergy(fdata = dft,
                                                   method = method,
                                                   ndatamin = ndatamin,
                                                   oncogene1 = gene_effect)$Synergy_score
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
              df_syng <- df_syng[.(disease,gene1,gene2,select_feature,effect1,effect2),
                                 specificity_pvalue :=max(Pvalue1, Pvalue2)]
            }
          }
        }
        #Calculate p.value for each immune cell count------------------------------------
        if( N_imcell_feature > 0){
          for( i in 1:N_imcell_feature){
            
            #Select feature
            select_feature <- sub_imcell_feature_list[i]
            mark_feature <- as.integer( which( colnames(data_cell) == select_feature))
            df_feature <- as.matrix(data_cell[ , mark_feature,drop=F])
            message( "          ...",select_feature)
            
            #Make a list of all genes with non-NA synergy score value
            df_syng_t <- df_syng_complete[df_syng_complete$Immune_feature==select_feature,]
            sub_genes1 <- unique( df_syng_t[ ,.( agent1 , agent1_expression)])
            sub_genes2 <- unique( df_syng_t[ ,.( agent2 , agent2_expression)])
            colnames(sub_genes1)<-c("gene","effect")
            colnames(sub_genes2)<-c("gene","effect")
            sub_genes <- unique( rbind(sub_genes1,sub_genes2))
            
            syng_dist <- vector( "list", nrow( sub_genes ))
            names(syng_dist) <- paste0(sub_genes$gene,"_",sub_genes$effect)
            
            df_bank_sub1 <- cbind(df_feature,
                                  df_bank_qr2[match(rownames(df_feature),rownames(df_bank_qr2)),])
            df_bank_sub1 <- df_bank_sub1[complete.cases(df_bank_sub1),,drop=F]
            bank_size <- as.numeric(ncol(df_bank_sub1))
            bank_mark <- N_genes+2
            df_bank_sub1_cols <- colnames( df_bank_sub1)
            
            #Build bootstrapping distribution
            tmpn <- nrow(sub_genes)
            pb <- txtProgressBar(min = 0, max = tmpn, char="-",style = 3)
            
            for( j in 1 :tmpn){
              
              gene_mark <- as.numeric( which( df_bank_sub1_cols == sub_genes$gene[ j]))
              df_bank_sub2 <- df_bank_sub1[ , c(1, gene_mark, bank_mark : bank_size)]
              df_bank_sub2 <- df_bank_sub2[ df_bank_sub2[ , 2 ] %in% c( 1, 4),,drop=F ]
              bank_size2 <- as.numeric( ncol( df_bank_sub2 ))
              gene_effect <- sub_genes$effect[j]
              
              syng_dist_t <- vector()
              for(k in 3:(bank_size2)){
                dft <- df_bank_sub2[,c(1,2,k)]
                dft <- dft[dft[,3] %in% c(1,4),,drop=F]
                syng_dist_t[k-2] <- find_a_synergy(fdata = dft,
                                                   method = method,
                                                   ndatamin = ndatamin,
                                                   oncogene1 = gene_effect)$Synergy_score
              }
              syng_dist_t <- syng_dist_t[complete.cases(syng_dist_t)]
              syng_dist[[j]] <- syng_dist_t
              setTxtProgressBar(pb, j)
            }
            
            #Calculate p.value for all pairs
            
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
              Pvalue1 <- sum( mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
              j2 <- which(names(syng_dist)==paste0(gene2,"_",effect2))
              tmp_dist <- syng_dist[[j2]]
              Pvalue2 <- sum( mysign*tmp_dist > mysign*myscore, tmp_dist==myscore )/length(tmp_dist)
              df_syng <- df_syng[.(disease,gene1,gene2,select_feature,effect1,effect2),
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
        df_syng_complete <- df_syng_complete[ which(df_syng_complete$Disease==disease), ]
        df_syng_complete1 <- df_syng_complete[df_syng_complete$Immune_feature %in%
                                                colnames(data_feature),]
        df_syng_complete2 <- df_syng_complete[df_syng_complete$Immune_feature %in%
                                                colnames(data_cell),]
        
        N_syng_complete1 <- nrow(df_syng_complete1)
        N_syng_complete2 <- nrow(df_syng_complete2)
        
        my_features1 <- unique(df_syng_complete1$Immune_feature)
        my_features_var1 <- my_features1[my_features1 %like% "score"]
        my_features_const1 <- my_features1[!(my_features1 %like% "score")]
        my_features2 <- unique(df_syng_complete2$Immune_feature)
        
        pb <- txtProgressBar(min = 0, max = N_iteration_sensitivity, char="-",style = 3)
        
        
        for(n_sns in 1:N_iteration_sensitivity){
          
          df_sub <- df_comb[sample(N_sample,N_sub,replace = F),]
          df_sub_qr <- get_quantile_rank(df_sub)
          dft <- data_expression[,colnames(data_expression) %in% rownames(df_sub)]
          
          if(N_syng_complete1>0){
            
            data_feature1 <- get_selected_features(dft,my_features_var1)
            data_feature_const1 <-data_feature[rownames(df_sub),my_features_const1,drop=F]
            if(length(data_feature_const1)>0){
              data_feature1 <- merge(data_feature1,data_feature_const1,by="")
            }
            for(pair_ID in 1:N_syng_complete1 ){
              tmp_df <- df_syng_complete1[pair_ID,]
              gene_ID1 <- tmp_df$agent1
              gene_ID2 <- tmp_df$agent2
              base_score <- tmp_df$Synergy_score
              effect1 <- tmp_df$agent1_expression
              effect2 <- tmp_df$agent2_expression
              feature_ID <- tmp_df$Immune_feature
              mark_feature <- as.integer( which( colnames( data_feature1) == feature_ID))
              df_feature <- as.matrix(data_feature1[ , mark_feature,drop=F])
              
              dft <- df_sub_qr[ , c(which(colnames(df_sub_qr)==gene_ID1),
                                    which(colnames(df_sub_qr)==gene_ID2))]
              dft <- dft[dft[ , 1] %in% c(1 , 4) ,,drop=F ]
              dft <- dft[dft[ , 2] %in% c(1 , 4) ,,drop=F ]
              
              dft <- cbind(df_feature[match(rownames(dft),rownames(df_feature)),,drop=F], dft)
              dft <- dft[complete.cases(dft),,drop=F]
              
              
              dfts <- find_a_synergy(fdata = dft,
                                     method = method,
                                     ndatamin = ndatamin,
                                     oncogene1 = effect1,
                                     oncogene2 = effect2)$Synergy_score
              
              
              df_syng_complete1$sum[pair_ID] <- sum(df_syng_complete1$sum[pair_ID],
                                                    dfts, na.rm = T)
              dfts <- dfts - base_score
              df_syng_complete1$sum2[pair_ID] <- sum(df_syng_complete1$sum2[pair_ID],
                                                     dfts*dfts, na.rm = T)
              df_syng_complete1$N[pair_ID] <- sum(df_syng_complete1$N[pair_ID],
                                                  abs(sign(dfts)), na.rm = T)
            }
          }
          if(N_syng_complete2>0){
            tmp <- as.data.frame(df_sub)
            tmp$PATIENT_BARCODE <- substr(rownames(tmp), 1, 12)
            tmp <- data.table::as.data.table(tmp)
            tmp <- tmp[,lapply(.SD,median),by=PATIENT_BARCODE]
            tmp<- as.data.frame(tmp)
            df_sub2 <- as.matrix(tmp[,-which(colnames(tmp)=="PATIENT_BARCODE")])
            rownames(df_sub2)<- tmp$PATIENT_BARCODE
            df_sub_qr2 <- get_quantile_rank(df_sub2)
            df_sub_qr <- as.matrix(df_sub_qr[,-1])
            df_sub_qr2 <- as.matrix(df_sub_qr2[,-1])
            
            
            for(pair_ID in 1:N_syng_complete2 ){
              
              tmp_df <- df_syng_complete2[pair_ID]
              gene_ID1 <- tmp_df$agent1
              gene_ID2 <- tmp_df$agent2
              base_score <- tmp_df$Synergy_score
              effect1 <- tmp_df$agent1_expression
              effect2 <- tmp_df$agent2_expression
              feature_ID <- tmp_df$Immune_feature
              mark_feature <- as.integer( which( colnames(data_cell) == feature_ID))
              df_feature <- as.matrix(data_cell[ , mark_feature,drop=F])
              
              dft <- df_sub_qr2[ , c(which(colnames(df_sub_qr2)==gene_ID1),
                                     which(colnames(df_sub_qr2)==gene_ID2))]
              dft <- dft[dft[ , 1] %in% c(1 , 4) ,,drop=F ]
              dft <- dft[dft[ , 2] %in% c(1 , 4) ,,drop=F ]
              dft <- cbind(df_feature[match(rownames(dft),rownames(df_feature)),,drop=F], dft)
              dft <- dft[complete.cases(dft),,drop=F]
              dfts <- find_a_synergy(fdata = dft,
                                     method = method,
                                     ndatamin = ndatamin,
                                     oncogene1 = effect1,
                                     oncogene2 = effect2)$Synergy_score
              
              df_syng_complete2$sum[pair_ID] <- sum(df_syng_complete2$sum[pair_ID],
                                                    dfts, na.rm = T)
              dfts <- dfts -base_score
              df_syng_complete2$sum2[pair_ID] <-sum(df_syng_complete2$sum2[pair_ID],
                                                    dfts*dfts, na.rm = T)
              df_syng_complete2$N[pair_ID] <- sum(df_syng_complete2$N[pair_ID],
                                                  abs(sign(dfts)), na.rm = T)
              
            }
          }
          setTxtProgressBar(pb, n_sns)
        }
        if(nrow(df_syng_complete1)>0 ){
          if(nrow(df_syng_complete2)>0 ){
            df_syng_complete <- rbind(df_syng_complete1,df_syng_complete2)
          }else{
            df_syng_complete <- df_syng_complete1
          }
        }else{
          df_syng_complete <- df_syng_complete2
        }
        df_syng_complete$sensitivity_R <- sqrt(df_syng_complete$sum2/df_syng_complete$N) / abs(df_syng_complete$Synergy_score)
        df_syng_complete$sum <- NULL
        df_syng_complete$sum2 <- NULL
        df_syng_complete$N <- NULL
        rownames(df_syng_complete) <- NULL
        df_syng<- df_syng[df_syng_complete, sensitivity_R := i.sensitivity_R]
      }
    }
    df_syng <- as.data.frame(df_syng)
    df_syng_all <- bind_rows(df_syng_all,df_syng )
    message("\n",disease," completed\n")
  }
  df_syng_all[df_syng_all=="DKFZp686O24166"]<-"NCR3LG1"
  df_syng_all[df_syng_all=="C10orf54"]<-"VSIR"
  colnames(df_syng_all)[1:7] <- c("Disease",
                                  "Gene1","Gene2",
                                  "IAP","Synergy_score",
                                  "Gene1_expression","Gene2_expression")
  return(df_syng_all)
}
