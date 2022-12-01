###############################################################################

# rscript to modify the .ind file for Ireland and England

ind_data <- "nl-gnomad-merged.ind"
ind_data <- read.table(file = ind_data, header = F)
colnames(ind_data) <- c("IID","SEX","Q")

split_df <- data.frame(do.call("rbind", strsplit(as.character(ind_data$IID), ":", fixed = TRUE)))
colnames(split_df) <- c("LABEL","IID")

results <- data.frame(IID = split_df$IID, 
                      SEX = ind_data$SEX,
                      LABEL = split_df$LABEL)
  
results$NEW_LABEL <- ifelse(results$LABEL == "KGP_YRI", "YRI", "Ignore")
results$NEW_LABEL <- ifelse(grepl("fs", results$LABEL), gsub(":","_",results$LABEL), results$NEW_LABEL)
results <- results[,c("IID","SEX","NEW_LABEL")]
  
outfile <- "nl-gnomad-merged.ascend-label.ind"
write.table(x = results, file = outfile, quote = F, col.names = F, row.names = F, sep = " ")
