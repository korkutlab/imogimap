#'  Calculation and statistical assessment of synergistic associations using TCGA data.
#'
#' @param df_syng a dataframe with 10 columns as outputed by im_syng_tcga. 
#' @param method A character string indicating which synergy score to be used. one of "max" or "independence".
#' Default is "max".
#'@param ndatamin minimum number of samples. Synergy score calculation will be skipped for matrices with number of rows less than ndatamin
#' @param N_iteration_sensitivity Number of iterations for random sampling for sensitivity analysis.
#' Default is 1000.
#' @param sample_list An optional character vector of TCGA samples barcodes indicating a subset of samples within a cohort.
#' @keywords specificity, pvalue ,bootstrapping
#' @return Specificty pvalues for each row of dataframe
#' @description
#' A helper function to allow user calculate specifity pvalue for a subset of its dataset. To calculate specifity pvalue for all data points set specificity = TRUE in im_syng_TCGA.
#' @details
#'A specificity p.value is computed using random sampling with replacement from two null models, generated from one of the two genes against a set of genes randomly selected from the genome. Two P-values are calculated for the synergistic interaction of the pair against the two null models. The highest of the two P-values is used to assess the specificity of the interaction against the whole genome. The number of randomly selected genes in each null model is determined by N_iteration_specificity.
#'
#'All barcodes in sample_list must be 15 character long and belong to the same cohort. When sample_list is provided, cohort should be the disease cohort that they belong to, otherwise only the first element of the cohort list will be used.
#'
#' @examples 
#' df <- im_syng_tcga(onco_gene=c("TGFB1","SERPINB9"), cohort=c("ucec"),specificity = F)
#' df <- get_sensitivity(df)
#' @importFrom dplyr bind_rows across mutate group_by everything distinct
#' @importFrom magrittr %>%
#' @importFrom data.table setkey as.data.table
#' @import curatedTCGAData
#' @importFrom data.table ":=" ".SD" "%like%"
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom stats complete.cases median
#'
#' @export

get_sensitivity <- function(df_syng,method='max',ndatamin=8,N_iteration_sensitivity=1000,sample_list){
  
  PATIENT_BARCODE <- Disease <- agent1 <- agent2 <- Immune_feature <- NULL
  agent1_expression <- agent2_expression <- NULL
  sensitivity_R <- i.sensitivity_R <- NULL
  
  
  if(missing(method)){
    method <- "max"
  }else{
    method <- tolower(method)
    if(!(method=="max" || method=="independence")){
      stop("ERROR: Method is not found. Please choose a method from: max or independence.")
    }
  }
  
  if(nrow(df_syng)>0){
    
    colnames(df_syng)[1:7] <- c("Disease", "agent1","agent2",
                                "Immune_feature","Synergy_score",
                                "agent1_expression","agent2_expression")
    cohort <- unique(df_syng$Disease)
    df_syng[df_syng=="NCR3LG1"]<-"DKFZp686O24166"
    df_syng[df_syng=="VSIR"]<-"C10orf54"
    
    df_syng <- data.table::as.data.table(df_syng)
    setkey(df_syng, Disease, agent1,agent2,Immune_feature)
    df_syng$sensitivity_R<- as.numeric(df_syng$sensitivity_R)
    
    
    for(cohortID in 1:length(cohort)){
      
      disease <- cohort[cohortID]
      df_syng_complete <- df_syng[ !is.na( df_syng$Synergy_score),]
      df_syng_complete <- df_syng_complete[ which(df_syng_complete$Disease==disease), ]
      df_syng_complete$sum <- 0
      df_syng_complete$sum2 <- 0
      df_syng_complete$N <- 0
      
      #Get expression data----------------------------
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
      onco_gene <- unique(c(df_syng_complete$agent1,df_syng_complete$agent2))
      onco_gene[onco_gene=="VSIR"]<- "C10orf54"
      onco_gene[onco_gene=="NCR3LG1"]<- "DKFZp686O24166"
      
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
    
      #Construct quantile ranking matrices for each sample--------
  
      df_selected <- scale(log2(df_selected+1),center = T,scale = T)
  
      df_comb <- df_selected
      N_sample <- as.numeric(nrow(df_comb))
      N_sub <- floor(N_sample*0.7)
      
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
          data_feature_const1 <- data_feature[rownames(df_sub),my_features_const1,drop=F]
          if(length(data_feature_const1) > 0 ){
            data_feature1 <- merge(data_feature1,data_feature_const1,by=0)
            rownames(data_feature1)<-data_feature1$Row.names
            data_feature1$Row.names<-NULL
            data_feature1<-as.matrix(data_feature1)
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
                                  which(colnames(df_sub_qr)==gene_ID2)),drop=F]
            dft <- dft[dft[ , 1] %in% c(1 , 4) ,,drop=F ]
            dft <- dft[dft[ , 2] %in% c(1 , 4) ,,drop=F ]
            
            dft <- cbind(df_feature[match(rownames(dft),rownames(df_feature)),,drop=F], dft)
            dft <- dft[complete.cases(dft),,drop=F]
            
            dfts <- find_a_synergy(fdata = dft,
                                   method = method,
                                   ndatamin = ndatamin,
                                   oncogene1 = effect1,
                                   oncogene2 = effect2)$Synergy_score
            
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

      df_syng_complete$sum2 <- NULL
      df_syng_complete$N <- NULL
      rownames(df_syng_complete) <- NULL
      df_syng<- df_syng[df_syng_complete, sensitivity_R := i.sensitivity_R]
    }
  }
  df_syng[df_syng=="DKFZp686O24166"]<-"NCR3LG1"
  df_syng[df_syng=="C10orf54"]<-"VSIR"
  colnames(df_syng)[1:7] <- c("Disease",
                              "Gene1","Gene2",
                              "IAP","Synergy_score",
                              "Gene1_expression","Gene2_expression")
  
  return(df_syng)
}