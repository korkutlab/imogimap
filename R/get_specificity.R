#'  Calculation and statistical assessment of synergistic associations using TCGA data.
#'
#' @param df_syng a dataframe with 10 columns as outputed by im_syng_tcga. 
#' @param method A character string indicating which synergy score to be used. one of "max" or "independence".
#' Default is "max".
#'@param ndatamin minimum number of samples. Synergy score calculation will be skipped for matrices with number of rows less than ndatamin
#' @param N_iteration_specificity Number of iterations for random sampling for specificity p.value calculation.
#' Default is 1000.
#' @keywords specificity, pvalue ,bootstrapping
#' @return Specificty pvalues for each row of dataframe
#' @description
#' A helper function to allow user calculate specifity pvalue for a subset of its dataset. To calculate specificity pvalue for all data points set specificity = TRUE in im_syng_TCGA.
#' @details
#'A specificity p.value is computed using random sampling with replacement from two null models, generated from one of the two genes against a set of genes randomly selected from the genome. Two P-values are calculated for the synergistic interaction of the pair against the two null models. The highest of the two P-values is used to assess the specificity of the interaction against the whole genome. The number of randomly selected genes in each null model is determined by N_iteration_specificity.
#'
#'
#' @examples 
#' df <- im_syng_tcga(onco_gene=c("TGFB1","SERPINB9"), cohort=c("ucec"),specificity = F)
#' df <- get_specificity(df)
#' @importFrom dplyr bind_rows across mutate group_by everything distinct
#' @importFrom magrittr %>%
#' @importFrom data.table setkey as.data.table
#' @import curatedTCGAData
#' @importFrom data.table ":=" ".SD"
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom stats complete.cases median
#'
#' @export

get_specificity <- function(df_syng,method='max',ndatamin=8,N_iteration_specificity=10000){
  
  if(missing(method)){
    method <- "max"
  }else{
    method <- tolower(method)
    if(!(method=="max" || method=="independence")){
      stop("ERROR: Method is not found. Please choose a method from: max or independence.")
    }
  }
  
  if(nrow(df_syng)>0){
    
    colnames(df_syng)[1:7] <- c("Disease","agent1","agent2",
                                "Immune_feature","Synergy_score",
                                "agent1_expression","agent2_expression")
    
    df_syng[df_syng$agent1=="VSIR",]$agent1 <- "C10orf54"
    df_syng[df_syng$agent1=="NCR3LG1",]$agent1 <- "DKFZp686O24166"
    df_syng[df_syng$agent2=="VSIR",]$agent2 <- "C10orf54"
    df_syng[df_syng$agent2=="NCR3LG1",]$agent2 <- "DKFZp686O24166"
    
    df_syng <- data.table::as.data.table(df_syng)
    setkey(df_syng, Disease, agent1,agent2,Immune_feature)
    
    cohort <- unique(df_syng$Disease)
    
    for(cohortID in 1:length(cohort)){
      
      disease <- cohort[cohortID]
      df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
      df_syng_complete <- df_syng_complete[ which(df_syng_complete$Disease==disease), ]
      
      #Get expression data----------------------------
      df <- curatedTCGAData::curatedTCGAData(diseaseCode = disease, version = "1.1.38",
                                             assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]
      
      data_expression <- df@assays$data@listData[[1]]
      colnames(data_expression)<-  substr(colnames(data_expression), 1, 15)
      
      #Construct IAPs-----------------------------------
      select_iap <- unique(df_syng_complete$Immune_feature)
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
      
      #Build unique permutations----------------------
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
          df_bank_sub1 <- df_bank_sub1[complete.cases(df_bank_sub1),]
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
  df_syng <- as.data.frame(df_syng)
  df_syng[df_syng=="DKFZp686O24166"]<-"NCR3LG1"
  df_syng[df_syng=="C10orf54"]<-"VSIR"
  colnames(df_syng)[1:7] <- c("Disease",
                                  "Gene1","Gene2",
                                  "IAP","Synergy_score",
                                  "Gene1_expression","Gene2_expression")
  return(df_syng)
  
}