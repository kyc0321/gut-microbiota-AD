## 00_setup.R — Environment, credentials, dataset discovery
suppressPackageStartupMessages({
  library(TwoSampleMR); library(MRPRESSO); library(ieugwasr)
  library(MendelianRandomization); library(dplyr); library(data.table); library(ggplot2)
})

## Run from the project root directory (the folder containing data/, results/, scripts/).
for (d in c("data","results","figures","logs")) dir.create(d, showWarnings = FALSE)

source("scripts/_helpers.R")  # reads OPENGWAS_JWT from the environment; no token is stored in this code
stopifnot(nchar(ieugwasr::get_opengwas_jwt()) > 100)
u <- tryCatch(ieugwasr::user(), error = function(e) list(error = conditionMessage(e)))
cat("Auth check:\n"); str(u)

cat("\n=== Periodontitis candidates ===\n")
perio <- ieugwasr::gwasinfo() %>%
  filter(grepl("periodont|gingiv", trait, ignore.case = TRUE)) %>%
  select(id, trait, year, author, sample_size, ncase, ncontrol, population, consortium) %>%
  arrange(desc(sample_size))
print(perio, n = 30); fwrite(perio, "data/datasets_periodontitis.csv")

cat("\n=== Alzheimer candidates ===\n")
ad <- ieugwasr::gwasinfo() %>%
  filter(grepl("alzheim", trait, ignore.case = TRUE)) %>%
  select(id, trait, year, author, sample_size, ncase, ncontrol, population, consortium) %>%
  arrange(desc(sample_size))
print(ad, n = 30); fwrite(ad, "data/datasets_alzheimer.csv")

cat("\n=== MiBioGen gut microbiota count ===\n")
gut <- ieugwasr::gwasinfo() %>%
  filter(grepl("MiBioGen|Kurilshikov", paste(consortium, author), ignore.case = TRUE)) %>%
  select(id, trait, year, author, sample_size, population, consortium)
cat("MiBioGen taxa found:", nrow(gut), "\n")
fwrite(gut, "data/datasets_gut_microbiota.csv")
print(head(gut, 10))

cat("\n=== Oral microbiome candidates ===\n")
oral <- ieugwasr::gwasinfo() %>%
  filter(grepl("oral|saliva|tongue", trait, ignore.case = TRUE) &
         grepl("microb|bacter|abundan", trait, ignore.case = TRUE)) %>%
  select(id, trait, year, author, sample_size, population, consortium)
cat("Oral microbiome studies found:", nrow(oral), "\n")
print(oral); fwrite(oral, "data/datasets_oral_microbiome.csv")

saveRDS(list(perio=perio, ad=ad, gut=gut, oral=oral), "data/datasets_index.rds")
cat("\n[OK] 00_setup.R complete.\n")
