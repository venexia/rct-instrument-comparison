# Set working directory =======================================================

setwd("E:/AntihypertensivesMR_letter/")
rm(list=ls())
graphics.off()

# Load packages ===============================================================

if (!require(pacman)) install.packages("pacman"); library(pacman)
p_load("readxl")
p_load("data.table")

# Source functions ============================================================

source("code/func_mRnd_binary.R")

# Define drug effects =========================================================

drugs <- rbind(c("Angiotensin converting enzyme inhibitors",21.14),
               c("Beta-adrenoceptor blockers",9.51),
               c("Calcium channel blockers",8.90),
               c("Thiazides and related diuretics",12.56))

drugs <- data.frame(drugs,stringsAsFactors = FALSE)

colnames(drugs) <- c("drug","drug_effect")
drugs$drug_effect <- as.numeric(drugs$drug_effect)

# Load data ===================================================================

df <- read_excel("raw/ije-2018-12-1541-File004.xlsx",
                 sheet = "ST3", skip = 1)

df <- df[df$drug %in% drugs$drug,c("drug","SNP","beta")]

# Restrict to SNPs used for analysis ==========================================

snps <- read_excel("raw/ije-2018-12-1541-File004.xlsx",
                   sheet = "ST4", skip = 1)

snps <- snps[snps$leveli %in% drugs$drug, c("leveli","snps.mr")]

s <- strsplit(snps$snps.mr, split = ";")

snps <- data.frame(drug = rep(snps$leveli, sapply(s, length)), 
                   SNP = unlist(s),
                   stringsAsFactors = FALSE)

df <- merge(snps,df,by=c("drug","SNP"),all.x = TRUE)

# Calculate absolute beta per 1mmHg ===========================================

df$beta <- abs(df$beta)/19.268

# Add drug effects ============================================================

df <- merge(df,drugs,by=c("drug"))

# Calculate approximate R^2 for each SNP ======================================

df$r2 <- df$beta/df$drug_effect

# Calculate total R^2 =========================================================

df <- aggregate(df$r2, by=list(drug = df$drug, drug_effect = df$drug_effect), FUN=sum)
colnames(df) <- c("drug","drug_effect","R2xz")

# Obtain case/control numbers =================================================

rct <- fread("raw/rct_estimates.csv",data.table = FALSE)
rct <- unique(rct[,c("Drug","Outcome","cases","controls")])
colnames(rct) <- c("drug","outcome","cases","controls")
df <- merge(rct,df,by = c("drug"),all.x = TRUE)

# Specify parameters for power calculations ===================================

df$alpha <- 0.05
df$N <- df$cases + df$controls
df$K <- df$cases/df$N

# Calculate power =============================================================

j <- 1e-4
df$min_power <- NA
df$max_power <- NA

for (i in 1:nrow(df)) {
  
  # Calculate minimum power =============================================================
  
  power <- 0
  or <- 1
  
  while(power < 0.8 & or > 0) {
    
    or <- or - j
    
    power <- mRnd_binary(alpha = df$alpha[i],
                         R2xz = df$R2xz[i],
                         N = df$N[i],
                         K = df$K[i],
                         OR = or)
    
  }
  
  df$min_power[i] <- ifelse(or<=j,NA,or)
  
  # Calculate maximum power =============================================================
  
  power <- 0
  or <- 1
  
  while(power < 0.8) {
    
    or <- or + j
    
    power <- mRnd_binary(alpha = df$alpha[i],
                         R2xz = df$R2xz[i],
                         N = df$N[i],
                         K = df$K[i],
                         OR = or)
    
  }
  
  df$max_power[i] <- or
  
}

# Scale estimates to drug effects =============================================

df$min_power_scaled <- ifelse(is.na(df$min_power),NA,sprintf("%.2f",round(exp((log(df$min_power)/19.27)*df$drug_effect),2)))
df$max_power_scaled <- ifelse(is.na(df$max_power),NA,sprintf("%.2f",round(exp((log(df$max_power)/19.27)*df$drug_effect),2)))

# Save ========================================================================

fwrite(df,"output/power.csv")