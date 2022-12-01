## Script to calculate the average within-cluster IBD sharing
## and the ROH overall and BINNED

###############################################################################
## read in the Ire-Brit NL IBD data and summarise the per-individual-pair and 
## per-cluster-pair amount

# process the ibd segment data with (3,15) limits
ibd_data <- list()
chroms <- 1:22
ibd_minmax <- c(3,15)

for (i in 1:length(chroms)) {
  chrom <- chroms[i]
  writeLines(paste0("Processing chr", chrom, "..."))
  ibd_data[[i]] <- "NL-IreBrit.chr",chrom,".mergeibd.ibd.gz"
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

# annotate the summarised IBD data
ibdsum_data <- ibdsum_data %>%
  left_join(., clust_data, by = c("ID1" = "IID")) %>% rename(CLST1 = CLST) %>%
  left_join(., clust_data, by = c("ID2" = "IID")) %>% rename(CLST2 = CLST)

# make the summarised data "symmetrical"
ibdsum_data <- bind_rows(ibdsum_data,
                         ibdsum_data %>% rename(ID1 = ID2, ID2 = ID1, CLST1 = CLST2, CLST2 = CLST1))

###############################################################################
## part five calculate the within pop levels
# summarise the within-cluster sharing
ibdwithin_data <- filter(ibdsum_data, CLST1 == CLST2) %>%
  select(-CLST2) %>%
  rename(CLST = CLST1) %>%
  group_by(CLST,ID1) %>%
  summarise(IndavTotLen = mean(TotLen), IndavNSEG = mean(NSEG), .groups = "drop") %>%
  group_by(CLST) %>%
  summarise(avTotLen = mean(IndavTotLen), CITotLen = sd(IndavTotLen)/sqrt(n()),
            avNSEG = mean(IndavNSEG), CINSEG = sd(IndavNSEG)/sqrt(n()), .groups = "drop")

# annotate
ibdwithin_data <- left_join(ibdwithin_data, param_data, by ="CLST") %>%
  arrange(LABEL_1, LABEL_2, LABEL_3)
  
###############################################################################
## Process the ROH/HBD data detected by refinedIBD
chroms <- 1:22
hbd_data <- list()

for (i in 1:length(chroms)) {
  writeLines(paste0(chroms[i]))
  hbd_data[[i]] <- NL-IreBrit.chr",chroms[i],".refinedIBD.hbd.gz"
  hbd_data[[i]] <- read_table(hbd_data[[i]], show_col_types = F,
                              col_names = c("IID","HAP","ID2","HAP2","CHROM","START","END","LOD","cM"))

}

hbd_data <- bind_rows(hbd_data) %>%
  select(IID,cM) %>%
  left_join(clust_data, by = "IID")

## now we remove outlier individuals (more than 6 SD greater than cluster-mean
hbdind_data <- group_by(hbd_data,CLST,IID) %>%
  summarise(TotHBD = sum(cM), .groups = "drop") %>%
  left_join(select(param_data, CLST, LABEL_1,LABEL_3), by = "CLST")

tmp_data <- group_by(hbdind_data, CLST) %>%
  summarise(avTotHBD = mean(TotHBD), sdHBD = sd(TotHBD), .groups = "drop")

hbd_outliers <- left_join(hbdind_data, tmp_data, by = "CLST") %>%
  filter(TotHBD > (avTotHBD + (6*sdHBD))) %>%
  mutate(CLST_ThresholdHBD = avTotHBD + (5*sdHBD)) %>%
  select(IID, LABEL_1,LABEL_3, TotHBD,CLST_ThresholdHBD)

## summarise the HBD without outliers
hbdsum_data <- bind_rows(
  hbd_data %>%
    filter(!(IID %in% hbd_outliers$IID)) %>%
    mutate(CLASS = "> 1 cM") %>%
    group_by(CLST,IID,CLASS) %>%
    summarise(TotHBD = sum(cM), NHBD = n(),
              .groups = "drop") %>%
    group_by(CLST,CLASS) %>%
    summarise(avTotHBD = mean(TotHBD), avNHBD = mean(NHBD),
              ciTotHBD = 1.96 * (sd(TotHBD)/sqrt(n())), ciNHBD = 1.96 * (sd(NHBD)/sqrt(n())),
              .groups = "drop"),
  hbd_data %>%
    filter(!(IID %in% hbd_outliers$IID)) %>%
    mutate(CLASS = ifelse(cM >= 1 & cM < 4, "[1,4)", NA),
           CLASS = ifelse(cM >= 4 & cM < 8, "[4,8)", CLASS),
           CLASS = ifelse(cM >= 8 & cM < 12, "[8,12)", CLASS),
           CLASS = ifelse(cM >= 12 & cM < 20, "[12,20)", CLASS),
           CLASS = ifelse(cM >= 20 & cM < 300, "[20,300)", CLASS)) %>%
    group_by(CLST,IID,CLASS) %>%
    summarise(TotHBD = sum(cM), NHBD = n(),
              .groups = "drop") %>%
    group_by(CLST,CLASS) %>%
    summarise(avTotHBD = mean(TotHBD), avNHBD = mean(NHBD),
              ciTotHBD = 1.96 * (sd(TotHBD)/sqrt(n())), ciNHBD = 1.96 * (sd(NHBD)/sqrt(n())),
              .groups = "drop")
) %>%
  left_join(select(param_data, CLST,LABEL_1,LABEL_2,LABEL_3), by = "CLST") %>%
  mutate(CLASS = factor(CLASS, levels = c("> 1 cM", "[1,4)", "[4,8)","[8,12)","[12,20)","[20,300)"))) %>%
  arrange(CLST,CLASS)
