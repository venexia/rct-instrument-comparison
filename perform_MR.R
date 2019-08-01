# Set working directory =======================================================

setwd("E:/AntihypertensivesMR_letter/")
rm(list=ls())
graphics.off()

# Load packages ===============================================================

if (!require(pacman)) install.packages("pacman"); library(pacman)
p_load_gh("MRCIEU/TwoSampleMR")
p_load("readxl")
p_load("data.table")

# Load RCT data ===============================================================

rct <- fread("raw/rct_estimates.csv",data.table = FALSE)
colnames(rct) <- c("Drug","Outcome","Approach","SNPs/Trials","RR","LCI","UCI","cases","controls")

# Load instruments for drugs with RCT data ===================================

instruments <- read_excel("raw/ije-2018-12-1541-File004.xlsx",
                          sheet = "ST3", skip = 1)

instruments <- instruments[instruments$drug %in% rct$Drug,]

# Prepare stroke data =========================================================

stroke <- fread("raw/MEGASTROKE.1.AS.EUR.out",data.table = FALSE)
stroke <- stroke[stroke$MarkerName %in% instruments$SNP,]
fwrite(stroke,"data/stroke.csv")

# Loop through exposure-outcome combinations ==================================

df <- NULL

for (i in unique(instruments$drug)) {
  
  for (j in c("Coronary heart disease","Stroke")) {
    
    res <- NULL
    
    # Format exposure data ====================================================
    
    exposure <- format_data(instruments[instruments$drug==i,], type = "exposure")
    
    # Clump exposure data =====================================================
    
    exposure <- clump_data(exposure)
    
    # Extract outcome data ====================================================
    
    if (j=="Coronary heart disease") {
      outcome <- extract_outcome_data(exposure$SNP, outcomes = c(7))
    }
    
    if (j=="Stroke") {
      outcome <- read_outcome_data(
        snps = exposure$SNP,
        filename = "data/stroke.csv",
        sep = ",",
        snp_col = "MarkerName",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col = "Freq1",
        pval_col = "P-value")
    }
    
    # Harmonise data ==========================================================
    
    dat <- harmonise_data(exposure,outcome)
    
    # Perform MR ==============================================================
    
    res <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw"))
  
    # Add results to dataframe ================================================
    
    res$exposure <- i
    res$outcome <- j
  
    df <- rbind(df,res)
    
  }
}

# Add confidence intervals ====================================================

df$lci.mr <- df$b - qnorm(0.975)*df$se
df$uci.mr <- df$b + qnorm(0.975)*df$se 

# Scale results ===============================================================
# Divide by sd of SBP in UK Biobank (19.268) and multiply by scale factor =====

scale <- rbind(c("Angiotensin converting enzyme inhibitors",21.14),
               c("Beta-adrenoceptor blockers",9.51),
               c("Calcium channel blockers",8.90),
               c("Thiazides and related diuretics",12.56))

scale <- data.frame(scale,stringsAsFactors = FALSE)

colnames(scale) <- c("exposure","scale_factor")
scale$scale_factor <- as.numeric(scale$scale_factor)

df <- merge(df,scale,by=c("exposure"))

df$or <- 1/exp(df$b*(df$scale_factor/19.268))
df$uci <- 1/exp(df$lci.mr*(df$scale_factor/19.268))
df$lci <- 1/exp(df$uci.mr*(df$scale_factor/19.268))

# Save MR estinates ===========================================================

df <- df[,c("exposure","outcome","method","nsnp","or","lci","uci","pval","scale_factor")]

fwrite(df, "data/mr_estimates.csv", row.names = FALSE)

# Create combined file of MR and RCT estimates =================================

df$Approach <- "MR"
df <- df[,c("exposure","outcome","Approach","nsnp","or","lci","uci")]
colnames(df) <- c("Drug","Outcome","Approach","SNPs/Trials","RR","LCI","UCI")

rct <- rct[,c("Drug","Outcome","Approach","SNPs/Trials","RR","LCI","UCI")]

df <- rbind(df,rct)

# Add labels for figure =======================================================

df$label <- paste0(sprintf("%.2f",round(df$RR,2)),
                   " (",sprintf("%.2f",round(df$LCI,2)),
                   "-",sprintf("%.2f",round(df$UCI,2)),")")

# Save output =================================================================

fwrite(df, "output/figure_data.csv", row.names = FALSE)