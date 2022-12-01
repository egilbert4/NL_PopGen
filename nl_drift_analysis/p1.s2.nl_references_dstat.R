# update the labels of the combined dataset

nl_labels <- "NL-ancestry.K22-clust-membership.csv"
gnomad_labels <- "gnomad_idfile.csv"
irebrit_labels <- "NL-IreBrit.IBD-clusters-nf_ref.tsv"

# read in the fam data
fam_data <- "nl-gnomad-irebrit-merged.fam"
fam_data <- read.table(file = fam_data, header = F, sep = " ", stringsAsFactors = F)
colnames(fam_data) <- c("DATA","IID","PID","MID","SEX","PHE")

# read in the nfl labels
nl_labels <- read.table(file = nl_labels, header = T, sep = ",", stringsAsFactors = F)
colnames(nl_labels) <- c("IID","NEW_IID","CLST","LABEL")
nl_labels$LABEL <- paste0("fs_",nl_labels$LABEL)
nl_labels <- nl_labels[,c("IID","NEW_IID","LABEL")]
nl_labels$NEW_IID <- paste(nl_labels$LABEL, nl_labels$NEW_IID, sep = ":")

# read in the gnomad labels
gnomad_labels <- read.table(file = gnomad_labels, header = T, sep = ",", stringsAsFactors = F)
colnames(gnomad_labels) <- c("ID","FID","IID","KGP3","HGDP","POP","GROUP","LABEL","gnomad_release")

gnomad_labels[gnomad_labels$KGP3,"DATA"] <- c("KGP")
gnomad_labels[gnomad_labels$HGDP,"DATA"] <- c("HGDP")
gnomad_labels$LABEL <- paste0(gnomad_labels$DATA, "_", gnomad_labels$LABEL)
gnomad_labels <- gnomad_labels[,c("ID","DATA","LABEL")]
colnames(gnomad_labels) <- c("IID","DATA","LABEL")

gnomad_labels$NEW_IID <- paste(gnomad_labels$DATA, sprintf("%04d", 1:nrow(gnomad_labels) ), sep = "")
gnomad_labels$NEW_IID <- paste(gnomad_labels$LABEL, gnomad_labels$NEW_IID, sep = ":")

# read in the ire-brit labels
irebrit_labels <- read.table(file = irebrit_labels, header = T, sep = "\t", stringsAsFactors = F)

irebrit_labels$c1 <- ifelse(irebrit_labels$CLST_1 == 1, irebrit_labels$CLST_2, irebrit_labels$CLST_1 + 1)
irebrit_labels$c2 <- ifelse(irebrit_labels$CLST_1 == 1, irebrit_labels$CLST_3, irebrit_labels$CLST_2)
irebrit_labels$c3 <- ifelse(irebrit_labels$CLST_1 == 1, irebrit_labels$CLST_4, irebrit_labels$CLST_3)

irebrit_labels <- irebrit_labels[irebrit_labels$c1 != 3, c("IID","c1","c2","c3")]
irebrit_labels$LABEL <- paste("IBD", irebrit_labels$c1, irebrit_labels$c2, irebrit_labels$c3, sep = "_")
irebrit_labels$NEW_IID <- irebrit_labels$IID
irebrit_labels <- irebrit_labels[,c("IID","NEW_IID","LABEL")]
irebrit_labels$NEW_IID <- paste(irebrit_labels$LABEL, irebrit_labels$NEW_IID, sep = ":")


label_data <- rbind(nl_labels[,c("IID","NEW_IID")],
                   gnomad_labels[,c("IID","NEW_IID")],
                   irebrit_labels[,c("IID","NEW_IID")])

label_data <- merge(fam_data, label_data, by = "IID", all.x = TRUE)
label_data <- na.omit(label_data)

outfile <- paste0(workdir, "nl-gnomad-irebrit-merged.keep.idfile")
write.table(x = label_data[,c("DATA","IID")], file = outfile, quote = F, col.names = F, row.names = F, sep = " ")

outfile <- paste0(workdir, "nl-gnomad-irebrit-merged.update.idfile")
write.table(x = label_data[,c("DATA","IID","DATA","NEW_IID")], file = outfile, quote = F, col.names = F, row.names = F, sep = " ")

