rm(list=ls()); gc()

library(readxl)
library(stringr)

# =========================
# 1) Read metadata (Excel)
# =========================
META_XLSX <- "~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx"
meta <- data.frame(read_excel(META_XLSX))

req <- c("sample", "collection", "groupe","sti", "stald", "hold", "phase", "study.day")

##check for missing data##
missing <- setdiff(req, colnames(meta))
if (length(missing) > 0) {
  stop(paste("Missing columns in metadata:", paste(missing, collapse = ", ")))
}

meta$sample <- as.character(meta$sample)
meta$groupe <- as.character(meta$groupe)
meta$dummy <- runif(NROW(meta))



# =========================
# 2) CHeck format collection date (WARN, keep rows)
#    Expected format: dd.mm.yyyy e.g. 21.04.2022
# =========================
meta$collection_date <- as.Date(trimws(as.character(meta$collection)), format = "%d.%m.%Y")

bad_idx <- which(is.na(meta$collection_date))
if (length(bad_idx) > 0) {
  warning("Some collection dates are missing/unparsable. They will be filled with the control date (2024-03-15).")
  message("Rows affected (sample, collection, groupe):")
  print(meta[bad_idx, c("sample", "collection", "groupe")])
}

meta=meta[-which(is.na(meta$collection_date)==T),] 




# =========================
# 4) Build MiMARKS-like metadata table (for NCBI upload)
# =========================

forMimarks <- data.frame(
  sample_name          = meta$sample,
  sample_title         = "",
  bioproject_accession = "",
  organism             = "unclassified sequences; metagenomes; organismal metagenomes; pig gut metagenome",
  
  collection_date      = meta$collection_date,
  
  env_broad_scale      = "ENVO:01000998",
  env_local_scale      = "ENVO:01001029",
  env_medium           = "UBERON:0001988",
  geo_loc_name         = "Denmark",
  host                 = "Sus scrofa domesticus",
  lat_lon              = "not applicable",
  groupe                = meta$groupe,
  phase                = meta$phase,
  study.day            = meta$study.day,
  pen                  = meta$sti,
  stald                = meta$stald,
  hold              = meta$hold,
  dummy                = meta$dummy,
  
  stringsAsFactors     = FALSE
)


# =========================
# 5) Write output files
# =========================
write.table(
  x = forMimarks,
  file = "forMimarks.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Optional: semicolon CSV for Excel
write.table(
  x = forMimarks,
  file = "forMimarks.csv",
  sep = ";",
  row.names = FALSE,
  quote = TRUE
)

tmp <- forMimarks
tmp$sample_name <- NULL
tmp$sample_title <- NULL
tmp$description <- NULL  # if you have it

sum(duplicated(tmp))

# =========================
# Build SRA metadata sheet (matches NCBI columns in your screenshot)
# =========================

# Your raw data directory (folder that contains many sample subfolders)
FASTQ_DIR <- "C:/Users/pearoro/OneDrive - Danmarks Tekniske Universitet/Desktop/Novogen resutls/Result_X204SC24051149-Z01-F002_16SV34/result_X204SC24051149-Z01-F002/01.RawData"

# List sample folders (one folder per sample)
sample_folders <- list.dirs(FASTQ_DIR, recursive = FALSE, full.names = FALSE)
if (length(sample_folders) == 0) stop("No sample folders found under FASTQ_DIR")

# Helper normalizer (so S3.1 matches S3_1 etc.)
norm <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))

folder_map <- data.frame(
  folder = sample_folders,
  folder_norm = norm(sample_folders),
  stringsAsFactors = FALSE
)

meta$sample <- as.character(meta$sample)
meta$sample_norm <- norm(meta$sample)

# Find _1 and _2 FASTQs inside a given sample folder
get_pair <- function(sample_id, sample_norm) {
  folder <- folder_map$folder[match(sample_norm, folder_map$folder_norm)]
  if (is.na(folder)) return(c(f1 = "", f2 = ""))
  
  files <- list.files(file.path(FASTQ_DIR, folder), pattern="\\.fastq\\.gz$", full.names = FALSE)
  
  f1 <- files[grepl("_1\\.fastq\\.gz$", files, ignore.case = TRUE)]
  f2 <- files[grepl("_2\\.fastq\\.gz$", files, ignore.case = TRUE)]
  
  # Take first if multiple
  f1 <- if (length(f1) >= 1) f1[1] else ""
  f2 <- if (length(f2) >= 1) f2[1] else ""
  
  # IMPORTANT: NCBI wants filenames that match what you upload.
  # If you upload keeping folders, use file.path(folder, f1). If you upload flattened, use just f1.
  USE_FOLDERS_ON_UPLOAD <- FALSE
  
  if (USE_FOLDERS_ON_UPLOAD) {
    c(f1 = gsub("\\\\","/", file.path(folder, f1)),
      f2 = gsub("\\\\","/", file.path(folder, f2)))
  } else {
    c(f1 = f1, f2 = f2)
  }
}

pairs <- t(mapply(get_pair, meta$sample, meta$sample_norm))
meta$F1 <- pairs[, "f1"]
meta$F2 <- pairs[, "f2"]

# Report missing matches
unmatched <- meta$sample[meta$F1 == "" | meta$F2 == ""]
if (length(unmatched) > 0) {
  message("UNMATCHED samples (no _1 or _2 fastq found):")
  print(head(unmatched, 100))
  message("Total unmatched: ", length(unmatched))
}

# Build the SRA table with EXACT column names (spaces included)
SRA <- data.frame(
  "sample_name"         = meta$sample,
  "library_ID"          = meta$sample,
  "title"               = paste0(meta$sample, " 16S rRNA amplicon"),
  " library_strategy"    = "AMPLICON",
  "library_source"      = "GENOMIC",
  "library_selection"   = "PCR",
  "library_layout"      = "paired",
  "platform"            = "ILLUMINA",
  "instrument_model"    = "Illumina NovaSeq 6000",
  "design_description"  = "16S rRNA gene amplicon sequencing (V3-V4) NOVOGENE",
  " filetype"            = "fastq",
  "filename"            = meta$F1,
  "filename2"           = meta$F2,
  "dummy"               = meta$dummy,
  check.names           = FALSE
)

# Write outputs
write.table(SRA, "SRA_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(SRA, "SRA_metadata.csv", sep = ";", row.names = FALSE, quote = TRUE)

message("Wrote: SRA_metadata.tsv (and SRA_metadata.csv)")
