###############################################################################
# merge the .map data with genetic map data and output

old_snp_data <- "temp.map"
old_snp_data <- read.table(old_snp_data, col.names = c("CHROM","SNPID","GEN","BP"))

chroms <- 1:22
gen_data <- list()
for (i in 1:length(chroms)) {
  gen_data[[i]] <- paste0("genetic_maps/build38/plink.chr", chroms[i], ".GRCh38.map")
  gen_data[[i]] <- read.table(gen_data[[i]], col.names = c("CHROM","BLANK","cM","BP"))
}

gen_data <- do.call(rbind, gen_data)

snp_data <- merge(old_snp_data[,c("SNPID","CHROM","BP")],
                  gen_data[,c("CHROM","BP","cM")],
                  by = c("CHROM","BP"), all.x = TRUE) 

snp_data <- snp_data[!is.na(snp_data$cM),]

outfile <- "temp.cm.map"
write.table(x = snp_data[,c("CHROM","SNPID","cM","BP")], file = outfile, sep = " ", row.names = F, col.names = F, quote = FALSE)

outfile <- "temp.cm.extract"
write.table(x = snp_data[,c("SNPID")], file = outfile, sep = " ", row.names = F, col.names = F, quote = FALSE)
