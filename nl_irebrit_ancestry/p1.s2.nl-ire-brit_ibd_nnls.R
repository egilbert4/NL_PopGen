## Script to process IBD segment data to estimate sharing proportions
## using nnls

###############################################################################
#### set up ####
# packages
library("tidyverse")
library("nnls")

###############################################################################
## General NNLS functions
## from authors of: https://www.nature.com/articles/nature14230

admix.nnls<-function(X,Y){
    ## Do an nnls admixture analysis (X ~ Y) where X is a vector
    if(class(X)!="numeric") stop("X must be numeric in admix.nnls")
    ourmix=getoverallfit(Y,X)$x
    ourmix=ourmix/sum(ourmix)
    ourmix
}

admix.nnls.all<-function(X,Y,verbose=TRUE){
    ## Do nnls admixture on each ROW of the MATRIX x
    #if(class(X)!="matrix") stop("X must be a matrix in admix.nnls.all")
    if (!(inherits(X, "matrix"))) stop("X must be numeric in admix.nnls")
    
    ret<-t(sapply(1:dim(X)[1],function(i){
        if(verbose) print(paste("Processing ind",i,"of",dim(X)[1]))
        admix.nnls(X[i,],Y)
    }))
    rownames(ret)<-rownames(X)
    ret
}

###############################################################################
# our the ibd nnls function
ibd_nnls <- function(popx, popsy, ibd_data) {
  # calculate X, a matrix of the per-target-individual contributions from the
  # 'source' populations, populations-y
  X <- ibd_data %>% 
    filter(CLST1 %in% popx & CLST2 %in% popy) %>%
    group_by(ID1) %>%
    mutate(cov = TotLen/sum(TotLen)) %>% ungroup() %>%
    group_by(ID1, CLST2) %>%
    summarise(sumcov = sum(cov), .groups = "drop") %>%
    mutate(CLST2 = paste0("C:",CLST2)) %>%
    pivot_wider(names_from = "CLST2", values_from = "sumcov", values_fill = 0) %>%
    column_to_rownames("ID1") %>%
    as.matrix()
  
  # calculate Y, the source pop contributions to themselves
  Y <- ibd_data %>% 
    filter(CLST1 %in% popy & CLST2 %in% popy) %>%
    group_by(ID1) %>%
    mutate(cov = TotLen/sum(TotLen)) %>% ungroup() %>%
    group_by(CLST1,ID1, CLST2) %>%
    summarise(sumcov = sum(cov), .groups = "drop") %>%
    group_by(CLST1, CLST2) %>%
    summarise(avcov = mean(sumcov), .groups = "drop") %>%
    mutate(CLST1 = paste0("C:",CLST1),
           CLST2 = paste0("C:",CLST2)) %>%
    pivot_wider(names_from = "CLST2", values_from = "avcov", values_fill = 0) %>%
    column_to_rownames("CLST1") %>%
    as.matrix()
  
  # perform nnls analysis
  results <- as_tibble(admix.nnls.all(X = as.matrix(X), Y = as.matrix(Y), verbose = FALSE)) %>%
    mutate(ID = rownames(X)) %>%
    pivot_longer(cols = starts_with("C:"),
                 names_to = "CLST2", values_to = "Prop") %>%
    mutate(CLST1 = popx) %>%
    select(CLST1,ID,CLST2,Prop) %>%
    rename(target = CLST1, source = CLST2)
  
  return(results)
}

###############################################################################
#### read in the cluster data ####
nl_clust <- "NL-ancestry.K22-clust-membership.csv"
nl_clust <- read_delim(file = nl_clust, delim = ",", col_names = T) %>%
  rename(CLST = K22) %>%
  arrange(CLST,IID) %>%
  select(IID,CLST) %>%
  mutate(CLST = paste0("fs:",CLST))

irebrit_clst <- "IreBrit.IBD-clusters-nf_ref.tsv"
irebrit_clst <- read_delim(file = irebrit_clst, delim = "\t", col_names = TRUE)

irebrit_clst <- irebrit_clst %>%
  mutate(CLST = paste0("IBD:",CLST))

clust_data <- bind_rows(nl_clust, irebrit_clst)
rm(nl_clust, irebrit_clst)

###############################################################################
#### process the ibd data to summarise per-pair sharing ####
# read in ibd_data
ibd_dir <- "ibd_data/"
ibd_data <- list()

chroms <- 1:22
ibd_minmax <- c(3,15)

for (i in 1:length(chroms)) {
  chrom <- chroms[i]
  writeLines(paste0("Processing chr", chrom, "..."))
  ibd_data[[i]] <- paste0(ibd_dir, "NL-IreBrit.chr",chrom,".mergeibd.ibd.gz")
  ibd_data[[i]] <- read_delim(file = ibd_data[[i]], delim = "\t",
                              col_names = c("ID1","HAP1","ID2","HAP2","CHROM","START","END","LOD","cM"))
  ibd_data[[i]] <- ibd_data[[i]] %>%
    filter(cM >= ibd_minmax[1] & cM <= ibd_minmax[2]) %>%
    select(ID1,ID2,cM)
}

(ibd_data <- do.call(rbind, ibd_data))

# summarise ibd_data
ibdsum_data <- ibd_data %>%
  group_by(ID1,ID2) %>%
  summarise(NSEG = n(), TotLen = sum(cM), .groups = "drop")

# annotate the summarised IBD data with ind-cluster
ibdsum_data <- ibdsum_data %>%
  left_join(., clust_data, by = c("ID1" = "IID")) %>% rename(CLST1 = CLST) %>%
  left_join(., clust_data, by = c("ID2" = "IID")) %>% rename(CLST2 = CLST)

#### save the summarised and annotated IBD data ####
outfile <- "NL-IreBrit.summarised-ibd.3-15.txt.gz"
write_delim(x = ibdsum_data, file = outfile, delim = "\t", col_names = T)

###############################################################################
#### read in the summarised and annotated IBD data ####
ibdsum_data <- "NL-IreBrit.summarised-ibd.3-15.txt.gz"
ibdsum_data <- read_delim(file = ibdsum_data, delim = "\t", col_names = T)

#### process summed IBD for the nnls method ####
# make the summarised data "symmetrical"
ibdsum_data <- bind_rows(ibdsum_data,
                         ibdsum_data %>% rename(ID1 = ID2, ID2 = ID1, CLST1 = CLST2, CLST2 = CLST1))

###############################################################################
#### perform the "all versus all" IBD-nnls analysis ####
# i.e., every NL/Ire/Brit cluster as a mixture of any else

# calculate nnls proportions across all clusters
 clusters <- clust_data %>%
   filter(!(CLST %in% paste0("fs:", 1:4))) %>%
   distinct(CLST) %>%
   arrange(CLST) %>%
   pull(CLST)
 
 tmp_nnls_data <- list()
 
 for (i in 1:length(clusters)) {
   popx <- clusters[i]
   popy <- clusters[clusters != popx]
   writeLines(paste0("Running ",popx,"..."))
   
   tmp_nnls_data[[i]] <- ibd_nnls(popx = popx, popsy = popy, ibd_data = ibdsum_data)
 }
 
 tmp_nnls_data <- do.call(rbind, tmp_nnls_data)
 
# summarise across the clusters
nnls_data <- tmp_nnls_data %>%
   group_by(target,source) %>%
   summarise(meanProp = mean(Prop),
             CIProp = (sd(Prop)/sqrt(n())*1.96),
             .groups = "drop") %>%
   mutate(source = gsub("^C:","",source))
rm(tmp_nnls_data)

# save this out
outfile <- "NL-IreBrit.ibd-nnls-all.results.txt"
write_delim(x = nnls_data, file = outfile, delim = "\t", col_names = T)

###############################################################################
#### now model NL as a mixture of Ire-Brit ####
# create inputs
nnls_irebrit_data <- list()
 
clusters <- clust_data %>%
  filter(!(CLST %in% paste0("fs:", 1:4))) %>%
  distinct(CLST) %>%
  arrange(CLST) %>%
  pull(CLST)
 
popsx <- clusters[grepl("fs:",clusters)]
popy <- clusters[!grepl("fs:",clusters)]
 
tmp_ibd_data <- ibdsum_data %>%
  filter(CLST1 %in% c(popsx, popy) & CLST2 %in% c(popsx, popy)) %>%
  group_by(ID1) %>%
  mutate(cov = TotLen/sum(TotLen))
 
# run the first pass of analyses
for (i in 1:length(popsx)) {
  popx <- popsx[i]
  writeLines(paste0("Running ",popx,"..."))
   
  nnls_irebrit_data[[i]] <- ibd_nnls(popx = popx, popsy = popy, ibd_data = tmp_ibd_data)
}
 
rm(tmp_ibd_data)
nnls_irebrit_data <- do.call(rbind, nnls_irebrit_data)
 
# filter out the sources that don't contribute >= 1%
popy <- nnls_irebrit_data %>%
  group_by(target, source) %>%
  summarise(avProp = mean(Prop), .groups = "drop") %>%
  group_by(source) %>%
  summarise(maxavProp = max(avProp), .groups = "drop") %>%
  filter(maxavProp >= 0.01) %>%
  mutate(source = gsub("C:","",source)) %>%
  pull(source)

# run again with the top contibutors to remove noise
nnls_irebrit_data <- list()
tmp_ibd_data <- ibdsum_data %>%
  filter(CLST1 %in% c(popsx, popy) & CLST2 %in% c(popsx, popy)) %>%
  group_by(ID1) %>%
  mutate(cov = TotLen/sum(TotLen))

# run the first pass of analyses
for (i in 1:length(popsx)) {
  popx <- popsx[i]
  writeLines(paste0("Running ",popx,"..."))
  
  nnls_irebrit_data[[i]] <- ibd_nnls(popx = popx, popsy = popy, ibd_data = tmp_ibd_data)
}
 
rm(tmp_ibd_data)
nnls_irebrit_data <- do.call(rbind, nnls_irebrit_data)
 
# now process for downstream
nnls_irebrit_data <- nnls_irebrit_data %>%
  group_by(target, source) %>%
  summarise(avProp = mean(Prop), .groups = "drop") %>%
  mutate(source = gsub("C:","",source))
 
# save this out
outfile <- "NL-IreBrit.ibd-nnls-irebritref.results.txt"
write_delim(x = nnls_irebrit_data, file = outfile, delim = "\t", col_names = T)
