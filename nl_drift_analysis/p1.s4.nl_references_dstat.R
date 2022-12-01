###############################################################################

## rscript to modify the .ind file to label individuals based on Ireland, England, fs_{5..22}

ind_data <- "NL-YRI-Brit.admixtools.ind"
ind_data <- read.table(file = ind_data, header = F)
colnames(ind_data) <- c("IID","SEX","Q")

split_df <- data.frame(do.call("rbind", strsplit(as.character(ind_data$IID), ":", fixed = TRUE)))
colnames(split_df) <- c("DATA","LABEL","IID")

results <- data.frame(IID = split_df$IID, 
                      SEX = ind_data$SEX,
                      LABEL = split_df$LABEL)
  
results$NEW_LABEL <- ifelse(results$LABEL == "KGP_YRI", "YRI", "Ignore")
results$NEW_LABEL <- ifelse(results$LABEL == "IBD_1_1_1", "England", results$NEW_LABEL)
results$NEW_LABEL <- ifelse(results$LABEL %in% c("IBD_2_2_4","IBD_2_1_2"), "Ireland", results$NEW_LABEL)
results$NEW_LABEL <- ifelse(grepl("fs", results$LABEL), gsub(":","_",results$LABEL), results$NEW_LABEL)
results <- results[,c("IID","SEX","NEW_LABEL")]
  
outfile <- "NL-YRI-Brit.admixtools.dstat-test1.ind"
write.table(x = results, file = outfile, quote = F, col.names = F, row.names = F, sep = " ")

###############################################################################

## update the labels for NL-cluster specific tests of D(YRI, NL-cluster; Ireland, England)

ind_data <- "NL-YRI-Brit.admixtools.ind"
ind_data <- read.table(file = ind_data, header = F)
colnames(ind_data) <- c("IID","SEX","Q")

split_df <- data.frame(do.call("rbind", strsplit(as.character(ind_data$IID), ":", fixed = TRUE)))
colnames(split_df) <- c("DATA","LABEL","IID")

results <- data.frame(IID = split_df$IID, 
                      SEX = ind_data$SEX,
                      LABEL = split_df$LABEL)
  
results$NEW_LABEL <- ifelse(results$LABEL == "KGP_YRI", "YRI", "Ignore")
results$NEW_LABEL <- ifelse(results$LABEL == "IBD_1_1_1", "England", results$NEW_LABEL)
results$NEW_LABEL <- ifelse(results$LABEL %in% c("IBD_2_2_4","IBD_2_1_2"), "Ireland", results$NEW_LABEL)
results$NEW_LABEL <- ifelse(grepl("fs", results$LABEL), gsub(":","_",results$LABEL), results$NEW_LABEL)
results <- results[,c("IID","SEX","NEW_LABEL")]
  
clusters <- paste0("fs_", 5:22)
for (j in 1:length(clusters)) {
    
  results$CLST_LABEL <- ifelse(results$NEW_LABEL %in% c("YRI","England","Ireland"), results$NEW_LABEL, "Ignore")
  results$CLST_LABEL <- ifelse(results$NEW_LABEL == clusters[j], clusters[j], results$CLST_LABEL)
  
  tmp_ids <- results[results$NEW_LABEL != clusters[j] & grepl("fs",results$NEW_LABEL),"IID"]
  tmp_ids <- sample(tmp_ids, 500, replace = FALSE)
    
  results$CLST_LABEL <- ifelse(results$IID %in% tmp_ids, "NL500", results$CLST_LABEL)
    
  outfile <- paste0("NL-YRI-Brit.", clusters[j], ".admixtools.dstat-test2.ind")
  write.table(x = results[,c("IID","SEX","CLST_LABEL")], file = outfile, quote = F, col.names = F, row.names = F, sep = " ")
    
  pop_data <- data.frame(POP1 = c("YRI","YRI"), POP2 = rep(clusters[j], 2),
                         POP3 = c("England","Ireland"), POP4 = c("NL500","NL500"))

  outfile <- paste0("NL-YRI-Brit.", clusters[j], ".admixtools.dstat-test2.popfile")
  write.table(x = pop_data, file = outfile, quote = F, col.names = F, row.names = F, sep = " ")
}
