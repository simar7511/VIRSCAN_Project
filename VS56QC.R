# Generate VirScan QC report at run level ----
library(dplyr)
library(openxlsx)
library(readr)
library(tidyr)

options(stringsAsFactors = F)

# modify workingDir to reflect your environment
workingDir <- "/Users/rbasom/Bioinformatics/VirScan/2023_SummerQC/"
setwd(workingDir)

poolName <- "VS56"
zScoreThreshold <- "7"
indexFileName <- "VS56 Genomics Indexes.xlsx"
flowCell <- "kyahya/230227_VH00699_270_AACJKGJM5"
libraryVersion <- "Jan2023beads"

if(libraryVersion == "Oct2021beads"){
  load("mockIP/oct2021beads.RData")
  libraryName <- "Oct 2021 (V3) + CoV"
  beadsOnlyTable <- oct2021beads
}

if(libraryVersion == "Jan2023beads"){
  load("mockIP/jan2023beads.RData")
  libraryName <- "Jan 2023 (V3) + CoV"
  beadsOnlyTable <- jan2023beads
}


d <- read.xlsx(paste0("input/indices/", indexFileName)) %>%
  select(Pool, Sample.Name) %>%
  rename(Pool.ID = Pool,
         Sample.ID = Sample.Name) %>%
  mutate(Sample.ID = sub("_rep", ".rep", Sample.ID)) %>%
  mutate(path = flowCell)

qcTable <- d %>%
  mutate(VirScanID = sub(".rep[1-2]$", "", Sample.ID),
         rep = sub("^.+rep", "", Sample.ID)) %>%
  select(VirScanID, rep, Pool.ID, Sample.ID, path) %>%
  mutate(path = sub("^.+/", "", path),
         pearsonCor = NA) %>%
  rename(PoolID = Pool.ID,
         SampleID = Sample.ID,
         Flowcell = path)

pools <- data.frame(Set = poolName)



# Vir3 QC table ----
qcTableVir3 <- qcTable %>% head(0)
for(i in 1:nrow(pools)){
  resultsPath <- paste0(workingDir, "input/Vir3/", pools$Set[i], "/")
  sampleAnnotation <- read_csv(paste0(resultsPath, "wide_data/virscan_sample_annotation_table.csv")) %>%
    data.frame() %>%
    mutate(SampleID = sub("_S.+$", "", sample_name),
           SampleID = sub("_rep", ".rep", SampleID)) %>%
    select(SampleID, sample, raw_total_sequences, reads_mapped, error_rate, average_quality) %>%
    rename(VirScanID = sample)
  qcTab <- qcTable %>%
    inner_join(sampleAnnotation, by = c("SampleID", "VirScanID"))
  # number of peptides with at least 15 counts
  # Larman CDI Labs ANTYGEN doc
  # requires at least 15 counts/peptide in order for peptide to be included in hits file
  sampleAnnotationForCounts <- read_csv(paste0(resultsPath, "wide_data/virscan_sample_annotation_table.csv")) %>%
    filter(control_status == "empirical") %>%
    mutate(sample_id = paste0("X", sample_id)) %>%
    mutate(SampleID = sub("_S.+$", "", sample_name),
           SampleID = sub("_rep", ".rep", SampleID)) %>%
    select(sample_id, SampleID) %>%
    data.frame()
  counts <- read_csv(paste0(resultsPath, "wide_data/virscan_counts.csv")) %>%
    data.frame() %>%
    select(-...1)
  counts_min15 <- counts
  counts_min15[] <- apply(counts_min15, 2, function(x) ifelse(x < 15, 0, 1))
  counts_min15sum <- apply(counts_min15, 2, sum) %>%
    data.frame()
  counts_min15sum <- data.frame(sample_id = row.names(counts_min15sum),
                                peptidesGTE_15 = counts_min15sum$.) %>%
    inner_join(sampleAnnotationForCounts, by = "sample_id") %>%
    select(-sample_id)
  # peptides with at least 1 count
  counts_min1 <- counts
  counts_min1[] <- apply(counts_min1, 2, function(x) ifelse(x < 1, 0, 1))
  counts_min1sum <- apply(counts_min1, 2, sum) %>%
    data.frame()
  counts_min1sum <- data.frame(sample_id = row.names(counts_min1sum),
                               peptidesDetected = counts_min1sum$.) %>%
    inner_join(sampleAnnotationForCounts, by = "sample_id") %>%
    select(-sample_id)

  # CPM correlation between replicate samples
  countsPerMillion <- read_csv(paste0(resultsPath, "wide_data/virscan_cpm.csv")) %>%
    data.frame() %>%
    select(-...1) %>%
    unique() %>%
    t() %>%
    data.frame()
  countsPerMillion2 <- data.frame(row.names(countsPerMillion), countsPerMillion) %>%
    rename(sample_id = row.names.countsPerMillion.)
  countsPerMillion3 <- sampleAnnotationForCounts %>%
    inner_join(countsPerMillion2, by = "sample_id") %>%
    mutate(VirScanID = sub("\\.rep.+$", "", SampleID)) %>%
    mutate(rep = sub("^.+\\.", "", SampleID)) %>%
    select(VirScanID, everything())
  VirScanIDs <- unique(countsPerMillion3$VirScanID)
  for(j in 1:length(VirScanIDs)){
    z <- countsPerMillion3 %>%
      filter(VirScanID == VirScanIDs[j])
    if(nrow(z) == 2){
      z2 <- z %>%
        select(-VirScanID, -sample_id, -SampleID) %>%
        t() %>%
        data.frame()
      pearsonCor <- cor(x = as.numeric(z2$X1), y = as.numeric(z2$X2), method = "pearson", use = "pairwise.complete.obs") %>% round(2)
      qcTab$pearsonCor[qcTab$VirScanID == VirScanIDs[j]] <- pearsonCor
    }
  }

  # average read depth, number of reads/number of peptides
  # there are 3919 oligos represented more than once in the Vir3 annotation file
  # for read depth, consider number of alignments since a single ready can be associated with multiple peptides
  alignments <- apply(counts, 2, sum) %>%
    data.frame()
  alignments <- data.frame(sample_id = row.names(alignments),
                           alignments = alignments$.) %>%
    inner_join(sampleAnnotationForCounts, by = "sample_id") %>%
    select(-sample_id)
  qcTab <- qcTab %>%
    left_join(counts_min15sum, by = "SampleID") %>%
    left_join(counts_min1sum, by = "SampleID") %>%
    left_join(alignments, by = "SampleID") %>%
    mutate(readDepth = round(alignments/peptidesDetected, 2)) %>%
    mutate(pctReadsAligned = 100 * round(reads_mapped/raw_total_sequences, 3)) %>%
    mutate(pctDetected = 100 * round(peptidesDetected/115754, 3))

  # select organism hits
  organismSummary <- gzfile(paste0(resultsPath, "/aggregated_data/organism.summary.csv.gz")) %>%
    read_csv() %>%
    data.frame()
  organismTotalHits <- organismSummary %>%
    select(sample, n_hits_all) %>%
    group_by(sample) %>%
    mutate(sumAllSpecies = sum(n_hits_all)) %>%
    data.frame() %>%
    select(sample, sumAllSpecies) %>%
    unique() %>%
    rename(VirScanID = sample)
  organismSummaryHits <- organismSummary %>%
    filter(organism %in% c("Human herpesvirus 1",
                           "Human herpesvirus 2",
                           "Human herpesvirus 3",
                           "Human herpesvirus 4",
                           "Human herpesvirus 5",
                           "Human respiratory syncytial virus",
                           "Influenza A virus",
                           "Influenza B virus",
                           "Rhinovirus A",
                           "Rhinovirus B",
                           "Streptococcus pneumoniae",
                           "Measles virus",
                           "Simian foamy virus",
                           "Aravan virus")) %>%
    select(sample, organism, n_hits_all) %>%
    spread(organism, n_hits_all) %>%
    rename(VirScanID = sample)
  qcTab <- qcTab %>%
    left_join(organismTotalHits, by = "VirScanID") %>%
    left_join(organismSummaryHits, by = "VirScanID")

  qcTableVir3 <- bind_rows(qcTableVir3, qcTab)
}

qcTableVir3$Library <- libraryName


qcTableVir3 <- qcTableVir3 %>%
  select(VirScanID, rep, PoolID, SampleID, Flowcell, raw_total_sequences, reads_mapped, pctReadsAligned,
         error_rate, average_quality, peptidesGTE_15, peptidesDetected, pctDetected, alignments, readDepth,
         pearsonCor, sumAllSpecies,
         `Human herpesvirus 1`, `Human herpesvirus 2`, `Human herpesvirus 3`, `Human herpesvirus 4`, `Human herpesvirus 5`,
         `Human respiratory syncytial virus`, `Influenza A virus`, `Influenza B virus`, `Rhinovirus A`, `Rhinovirus B`,
         `Streptococcus pneumoniae`, `Measles virus`, `Simian foamy virus`, `Aravan virus`,
         Library) %>%
  mutate(numericPoolID = as.numeric(sub("VS", "", PoolID)),
         sampleNumber = as.numeric(sub("^.+_", "", VirScanID)),
         rep = as.numeric(rep)) %>%
  arrange(numericPoolID, sampleNumber, rep) %>%
  select(-numericPoolID, -sampleNumber)


# Vir3 QC thresholds ----
vir3thresh_readsMapped <- 200000
vir3thresh_averageQuality <- 32
vir3thresh_pctReadsAligned <- 60
vir3thresh_pctDetected <- 30
vir3thresh_readDepth <- 4
vir3thresh_pearsonCor <- .25

# tally Vir3 QC metrics below threshold ----
qcTableVir3 <- qcTableVir3 %>%
  mutate(`reads < 200000` = case_when(reads_mapped > vir3thresh_readsMapped ~ 0,
                                      reads_mapped <= vir3thresh_readsMapped ~ 1),
         `read depth < 4` = case_when(readDepth > vir3thresh_readDepth ~ 0,
                                      readDepth <= vir3thresh_readDepth ~ 1,),
         `detected < 30%` = case_when(pctDetected > vir3thresh_pctDetected ~ 0,
                                      pctDetected <= vir3thresh_pctDetected ~ 1),
         `cor < 25%` = case_when(pearsonCor > vir3thresh_pearsonCor ~ 0,
                                 pearsonCor <= vir3thresh_pearsonCor ~ 1)) %>%
  group_by(SampleID) %>%
  mutate(`QC Score Sum` = sum(c_across(`reads < 200000`:`cor < 25%`))) %>%
  data.frame(check.names = F)
qcTableVir3 <- qcTableVir3 %>%
  select(VirScanID, rep, PoolID, SampleID, Flowcell, raw_total_sequences, reads_mapped, pctReadsAligned,
         error_rate, average_quality, peptidesGTE_15, peptidesDetected, pctDetected, alignments, readDepth,
         pearsonCor,
         `reads < 200000`, `read depth < 4`, `detected < 30%`, `cor < 25%`, `QC Score Sum`,
         sumAllSpecies,
         `Human herpesvirus 1`, `Human herpesvirus 2`, `Human herpesvirus 3`, `Human herpesvirus 4`, `Human herpesvirus 5`,
         `Human respiratory syncytial virus`, `Influenza A virus`, `Influenza B virus`, `Rhinovirus A`, `Rhinovirus B`,
         `Streptococcus pneumoniae`, `Measles virus`, `Simian foamy virus`, `Aravan virus`,
         Library)



# CoV QC table ----
qcTableCoV <- qcTable %>% head(0)
for(i in 1:nrow(pools)){
  resultsPath <- paste0(workingDir, "input/CoV/", pools$Set[i], "/")
  sampleAnnotation <- read_csv(paste0(resultsPath, "wide_data/virscan_sample_annotation_table.csv")) %>%
    data.frame() %>%
    mutate(SampleID = sub("_S.+$", "", sample_name),
           SampleID = sub("_rep", ".rep", SampleID)) %>%
    select(SampleID, sample, raw_total_sequences, reads_mapped, error_rate, average_quality) %>%
    rename(VirScanID = sample)
  qcTab <- qcTable %>%
    inner_join(sampleAnnotation, by = c("SampleID", "VirScanID"))
  # number of peptides with at least 15 counts
  # Larman CDI Labs ANTYGEN doc
  # requires at least 15 counts/peptide in order for peptide to be included in hits file
  sampleAnnotationForCounts <- read_csv(paste0(resultsPath, "wide_data/virscan_sample_annotation_table.csv")) %>%
    filter(control_status == "empirical") %>%
    mutate(sample_id = paste0("X", sample_id)) %>%
    mutate(SampleID = sub("_S.+$", "", sample_name),
           SampleID = sub("_rep", ".rep", SampleID)) %>%
    select(sample_id, SampleID) %>%
    data.frame()
  counts <- read_csv(paste0(resultsPath, "wide_data/virscan_counts.csv")) %>%
    data.frame() %>%
    select(-...1)
  counts_min15 <- counts
  counts_min15[] <- apply(counts_min15, 2, function(x) ifelse(x < 15, 0, 1))
  counts_min15sum <- apply(counts_min15, 2, sum) %>%
    data.frame()
  counts_min15sum <- data.frame(sample_id = row.names(counts_min15sum),
                                peptidesGTE_15 = counts_min15sum$.) %>%
    inner_join(sampleAnnotationForCounts, by = "sample_id") %>%
    select(-sample_id)
  # peptides with at least 1 count
  counts_min1 <- counts
  counts_min1[] <- apply(counts_min1, 2, function(x) ifelse(x < 1, 0, 1))
  counts_min1sum <- apply(counts_min1, 2, sum) %>%
    data.frame()
  counts_min1sum <- data.frame(sample_id = row.names(counts_min1sum),
                               peptidesDetected = counts_min1sum$.) %>%
    inner_join(sampleAnnotationForCounts, by = "sample_id") %>%
    select(-sample_id)

  # CPM correlation between replicate samples
  countsPerMillion <- read_csv(paste0(resultsPath, "wide_data/virscan_cpm.csv")) %>%
    data.frame() %>%
    select(-...1) %>%
    unique() %>%
    t() %>%
    data.frame()
  countsPerMillion2 <- data.frame(row.names(countsPerMillion), countsPerMillion) %>%
    rename(sample_id = row.names.countsPerMillion.)
  countsPerMillion3 <- sampleAnnotationForCounts %>%
    inner_join(countsPerMillion2, by = "sample_id") %>%
    mutate(VirScanID = sub("\\.rep.+$", "", SampleID)) %>%
    mutate(rep = sub("^.+\\.", "", SampleID)) %>%
    select(VirScanID, everything())
  VirScanIDs <- unique(countsPerMillion3$VirScanID)
  for(j in 1:length(VirScanIDs)){
    z <- countsPerMillion3 %>%
      filter(VirScanID == VirScanIDs[j])
    if(nrow(z) == 2){
      z2 <- z %>%
        select(-VirScanID, -sample_id, -SampleID) %>%
        t() %>%
        data.frame()
      pearsonCor <- cor(x = as.numeric(z2$X1), y = as.numeric(z2$X2), method = "pearson", use = "pairwise.complete.obs") %>% round(2)
      qcTab$pearsonCor[qcTab$VirScanID == VirScanIDs[j]] <- pearsonCor
    }
  }

  alignments <- apply(counts, 2, sum) %>%
    data.frame()
  alignments <- data.frame(sample_id = row.names(alignments),
                           alignments = alignments$.) %>%
    inner_join(sampleAnnotationForCounts, by = "sample_id") %>%
    select(-sample_id)
  qcTab <- qcTab %>%
    left_join(counts_min15sum, by = "SampleID") %>%
    left_join(counts_min1sum, by = "SampleID") %>%
    left_join(alignments, by = "SampleID") %>%
    mutate(readDepth = round(alignments/peptidesDetected, 2)) %>%
    mutate(pctReadsAligned = 100 * round(reads_mapped/raw_total_sequences, 3)) %>%
    mutate(pctDetected = 100 * round(peptidesDetected/6932, 3))

  # select organism hits
  organismSummary <- gzfile(paste0(resultsPath, "/aggregated_data/organism.summary.csv.gz")) %>%
    read_csv() %>%
    data.frame()
  organismTotalHits <- organismSummary %>%
    select(sample, n_hits_all) %>%
    group_by(sample) %>%
    mutate(sumAllSpecies = sum(n_hits_all)) %>%
    data.frame() %>%
    select(sample, sumAllSpecies) %>%
    unique() %>%
    rename(VirScanID = sample)
  organismSummaryHits <- organismSummary %>%
    select(sample, organism, n_hits_all) %>%
    spread(organism, n_hits_all) %>%
    rename(VirScanID = sample)
  qcTab <- qcTab %>%
    left_join(organismTotalHits, by = "VirScanID") %>%
    left_join(organismSummaryHits, by = "VirScanID")



  qcTableCoV <- bind_rows(qcTableCoV, qcTab)
}
qcTableCoV$Library <- libraryName

qcTableCoV <- qcTableCoV %>%
  select(VirScanID, rep, PoolID, SampleID, Flowcell, raw_total_sequences, reads_mapped, pctReadsAligned,
         error_rate, average_quality, peptidesGTE_15, peptidesDetected, pctDetected, alignments, readDepth,
         pearsonCor, sumAllSpecies,
         `Bat coronavirus 279/2005 (BtCoV) (BtCoV/279/2005)`,
         `Bat coronavirus HKU3 (BtCoV) (SARS-like coronavirus HKU3)`,
         `Bat coronavirus Rp3/2004 (BtCoV/Rp3/2004) (SARS-like coronavirus Rp3)`,
         `Human coronavirus 229E (HCoV-229E)`,
         `Human coronavirus HKU1 (isolate N2) (HCoV-HKU1)`,
         `Human coronavirus NL63 (HCoV-NL63)`,
         `Human coronavirus OC43 (HCoV-OC43)`,
         `Human SARS coronavirus (SARS-CoV) (Severe acute respiratory syndrome coronavirus)`,
         `Middle East respiratory syndrome-related coronavirus`,
         `Severe acute respiratory syndrome coronavirus 2`,
         `Wuhan seafood market pneumonia virus`,
         Library) %>%
  mutate(numericPoolID = as.numeric(sub("VS", "", PoolID)),
         sampleNumber = as.numeric(sub("^.+_", "", VirScanID)),
         rep = as.numeric(rep)) %>%
  arrange(numericPoolID, sampleNumber, rep) %>%
  select(-numericPoolID, -sampleNumber)


# CoV QC thresholds
covthresh_readsMapped <- 10000
covthresh_averageQuality <- 32
covthresh_pctReadsAligned <- 2
covthresh_pctDetected <- 20
covthresh_readDepth <- 3
covthresh_pearsonCor <- .25

# tally CoV QC metrics below threshold
qcTableCoV.tally <- qcTableCoV

# tally CoV QC metrics below threshold ----
qcTableCoV <- qcTableCoV %>%
  mutate(`reads < 10000` = case_when(reads_mapped > covthresh_readsMapped ~ 0,
                                     reads_mapped <= covthresh_readsMapped ~ 1),
         `read depth < 4` = case_when(readDepth > covthresh_readDepth ~ 0,
                                      readDepth <= covthresh_readDepth ~ 1,),
         `detected < 20%` = case_when(pctDetected > covthresh_pctDetected ~ 0,
                                      pctDetected <= covthresh_pctDetected ~ 1),
         `cor < 25%` = case_when(pearsonCor > covthresh_pearsonCor ~ 0,
                                 pearsonCor <= covthresh_pearsonCor ~ 1)) %>%
  group_by(SampleID) %>%
  mutate(`QC Score Sum` = sum(c_across(`reads < 10000`:`cor < 25%`))) %>%
  data.frame(check.names = F)

qcTableCoV <- qcTableCoV %>%
  select(VirScanID, rep, PoolID, SampleID, Flowcell, raw_total_sequences, reads_mapped, pctReadsAligned,
         error_rate, average_quality, peptidesGTE_15, peptidesDetected, pctDetected, alignments, readDepth,
         pearsonCor,
         `reads < 10000`, `read depth < 4`, `detected < 20%`, `cor < 25%`, `QC Score Sum`,
         sumAllSpecies,
         `Bat coronavirus 279/2005 (BtCoV) (BtCoV/279/2005)`,
         `Bat coronavirus HKU3 (BtCoV) (SARS-like coronavirus HKU3)`,
         `Bat coronavirus Rp3/2004 (BtCoV/Rp3/2004) (SARS-like coronavirus Rp3)`,
         `Human coronavirus 229E (HCoV-229E)`,
         `Human coronavirus HKU1 (isolate N2) (HCoV-HKU1)`,
         `Human coronavirus NL63 (HCoV-NL63)`,
         `Human coronavirus OC43 (HCoV-OC43)`,
         `Human SARS coronavirus (SARS-CoV) (Severe acute respiratory syndrome coronavirus)`,
         `Middle East respiratory syndrome-related coronavirus`,
         `Severe acute respiratory syndrome coronavirus 2`,
         `Wuhan seafood market pneumonia virus`,
         Library)



# arrange by pool and rep
qcTableVir3 <- qcTableVir3 %>%
  mutate(PoolIndex = sub("VS","", VirScanID),
         PoolIndex = as.numeric(sub("_[0-9]+$", "", PoolIndex))) %>%
  mutate(SampleIndex = sub("VS[0-9]+_", "", SampleID),
         SampleIndex = as.numeric(sub(".rep[1-2]", "", SampleIndex))) %>%
  arrange(PoolIndex, SampleIndex, rep) %>%
  select(-PoolIndex, -SampleIndex)

qcTableCoV <- qcTableCoV %>%
  mutate(PoolIndex = sub("VS","", VirScanID),
         PoolIndex = as.numeric(sub("_[0-9]+$", "", PoolIndex))) %>%
  mutate(SampleIndex = sub("VS[0-9]+_", "", SampleID),
         SampleIndex = as.numeric(sub(".rep[1-2]", "", SampleIndex))) %>%
  arrange(PoolIndex, SampleIndex, rep) %>%
  select(-PoolIndex, -SampleIndex)


# write Vir3 and CoV QC tables to a single Excel file ----
# highlight problematic values
highlightStyle <- createStyle(bgFill = "#FFFF00")
tooLowStyle <- createStyle(fontSize = 12, fontColour = "red", textDecoration = "bold")
maxOverallStyle <- createStyle(fontSize = 12, fontColour = "skyblue", textDecoration = "bold")
bordersAllStyle <- createStyle(border = "TopBottomLeftRight")
bordersAllDateStyle <- createStyle(border = "TopBottomLeftRight", numFmt = "DATE")
wb <- createWorkbook()

# Vir3 worksheet ----
addWorksheet(wb, "Vir3")
writeData(wb, "Vir3", qcTableVir3)
freezePane(wb, "Vir3", firstRow = TRUE)
conditionalFormatting(wb,
                      "Vir3",
                      cols = grep("reads_mapped", names(qcTableVir3)),
                      rows = 1:(nrow(qcTableVir3) + 1),
                      rule = paste0("<= ", vir3thresh_readsMapped),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "Vir3",
                      cols = grep("average_quality", names(qcTableVir3)),
                      rows = 1:(nrow(qcTableVir3) + 1),
                      rule = paste0("<= ", vir3thresh_averageQuality),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "Vir3",
                      cols = grep("pctReadsAligned", names(qcTableVir3)),
                      rows = 1:(nrow(qcTableVir3) + 1),
                      rule = paste0("<= ", vir3thresh_pctReadsAligned),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "Vir3",
                      cols = grep("pctDetected", names(qcTableVir3)),
                      rows = 1:(nrow(qcTableVir3) + 1),
                      rule = paste0("<= ", vir3thresh_pctDetected),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "Vir3",
                      cols = grep("readDepth", names(qcTableVir3)),
                      rows = 1:(nrow(qcTableVir3) + 1),
                      rule = paste0("<= ", vir3thresh_readDepth),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "Vir3",
                      cols = grep("pearsonCor", names(qcTableVir3)),
                      rows = 1:(nrow(qcTableVir3) + 1),
                      rule = paste0("<= ", vir3thresh_pearsonCor),
                      style = tooLowStyle)

setColWidths(wb,
             "Vir3",
             cols = 1:ncol(qcTableVir3))
addStyle(wb,
         "Vir3",
         style = bordersAllStyle,
         cols = 1:(ncol(qcTableVir3)),
         rows = 1:(nrow(qcTableVir3) + 1),
         gridExpand = TRUE)


# CoV worksheet ----
addWorksheet(wb, "CoV")
writeData(wb, "CoV", qcTableCoV)
freezePane(wb, "CoV", firstRow = TRUE)
conditionalFormatting(wb,
                      "CoV",
                      cols = grep("reads_mapped", names(qcTableCoV)),
                      rows = 1:(nrow(qcTableCoV) + 1),
                      rule = paste0("<= ", covthresh_readsMapped),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "CoV",
                      cols = grep("average_quality", names(qcTableCoV)),
                      rows = 1:(nrow(qcTableCoV) + 1),
                      rule = paste0("<= ", covthresh_averageQuality),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "CoV",
                      cols = grep("pctReadsAligned", names(qcTableCoV)),
                      rows = 1:(nrow(qcTableCoV) + 1),
                      rule = paste0("<= ", covthresh_pctReadsAligned),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "CoV",
                      cols = grep("pctDetected", names(qcTableCoV)),
                      rows = 1:(nrow(qcTableCoV) + 1),
                      rule = paste0("<= ", covthresh_pctDetected),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "CoV",
                      cols = grep("readDepth", names(qcTableCoV)),
                      rows = 1:(nrow(qcTableCoV) + 1),
                      rule = paste0("<= ", covthresh_readDepth),
                      style = tooLowStyle)
conditionalFormatting(wb,
                      "CoV",
                      cols = grep("pearsonCor", names(qcTableCoV)),
                      rows = 1:(nrow(qcTableCoV) + 1),
                      rule = paste0("<= ", covthresh_pearsonCor),
                      style = tooLowStyle)

setColWidths(wb,
             "CoV",
             cols = 1:ncol(qcTableCoV))#,
addStyle(wb,
         "CoV",
         style = bordersAllStyle,
         cols = 1:(ncol(qcTableCoV)),
         rows = 1:(nrow(qcTableCoV) + 1),
         gridExpand = TRUE)


# beads only worksheet ----
addWorksheet(wb, libraryVersion)
writeData(wb, libraryVersion, beadsOnlyTable)
freezePane(wb, libraryVersion, firstRow = TRUE)


# write report to Excel file
todaysDate <- gsub("-|$", ".", Sys.Date())
saveWorkbook(wb, file = paste0("qcReports/", todaysDate, poolName, "_Z", zScoreThreshold, "_QC.xlsx"))
