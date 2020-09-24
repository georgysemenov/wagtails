library(tidyverse)
data <- read.table("compound.genotypes1A.20.txt", header=T)
#this dataset includes only individuals sampled within alba and personata hybrid zone
# head.plumage.value is a measure of the amount of white coloration on the head and neck plumage (see Semenov et al. 2017 Mol. Ecol. for plumage quantification details)
# Input file contains genotypic information for two genomic regions on chromsome 1A (aa=alba allele, aP=heterozygote, PP=personata allele) and 20 (AA=alba allele, Ap=heterozygote,
# pp=personata alleles) encoded as biallelic loci. Alleles with partially dominant inheritance are capitilized.

# Model with addiditive effects in both 1A and 20
data$add1A <- gsub("PP", 1, data$genotype.1A)
data$add1A <- gsub("aP", 0.5, data$add1A)
data$add1A <- gsub("aa", 0, data$add1A)
data$add1A <- as.numeric(data$add1A)

data$add20 <- gsub("AA", 0, data$genotype.20)
data$add20 <- gsub("Ap", 0.5, data$add20)
data$add20 <- gsub("pp", 1, data$add20)
data$add20 <- as.numeric(data$add20)

additive <- lm(head.plumage.value ~ add1A + add20, data = data)

# Partially dominant model where:
# 1A is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model) and
# 20 is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model)
data$pers.dom.1A <- gsub("PP", 1, data$genotype.1A)
data$pers.dom.1A <- gsub("aP", 0.75, data$pers.dom.1A)
data$pers.dom.1A <- gsub("aa", 0, data$pers.dom.1A)
data$pers.dom.1A <- as.numeric(data$pers.dom.1A)

data$alba.dom20 <- gsub("AA", 1, data$genotype.20)
data$alba.dom20 <- gsub("Ap", 0.75, data$alba.dom20)
data$alba.dom20 <- gsub("pp", 0, data$alba.dom20)
data$alba.dom20 <- as.numeric(data$alba.dom20)

pers.dom.1A_alba.dom20 <- lm(head.plumage.value ~ pers.dom.1A + alba.dom20, data = data)

# Make a partially dominant model where 
# 1A is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model) and
# 20 is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model)
data$alba.dom1A <- gsub("PP", 0, data$genotype.1A)
data$alba.dom1A <- gsub("aP", 0.75, data$alba.dom1A)
data$alba.dom1A <- gsub("aa", 1, data$alba.dom1A)
data$alba.dom1A <- as.numeric(data$alba.dom1A)

data$pers.dom20 <- gsub("AA", 0, data$genotype.20)
data$pers.dom20 <- gsub("Ap", 0.75, data$pers.dom20)
data$pers.dom20 <- gsub("pp", 1, data$pers.dom20)
data$pers.dom20 <- as.numeric(data$pers.dom20)

alba.dom.1A_pers.dom20 <- lm(head.plumage.value ~ alba.dom1A + pers.dom20 , data = data)

# Make a partially dominant model where 
# 1A is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model)
# 20 is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model)

alba.dom.1A_alba.dom20 <- lm(head.plumage.value ~ alba.dom1A + alba.dom20 , data = data)

# Make a partially dominant model where 
# 1A is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model)
# 20 is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model)

pers.dom.1A_pers.dom20 <- lm(head.plumage.value ~ pers.dom.1A + pers.dom20 , data = data)

AIC(additive, pers.dom.1A_alba.dom20, alba.dom.1A_pers.dom20, alba.dom.1A_alba.dom20, pers.dom.1A_pers.dom20)

#                        df      AIC
# additive                4 442.3477
# pers.dom.1A_alba.dom20  4 427.3859
# alba.dom.1A_pers.dom20  4 478.7637
# alba.dom.1A_alba.dom20  4 441.5010
# pers.dom.1A_pers.dom20  4 475.7174

AICc(additive, pers.dom.1A_alba.dom20, alba.dom.1A_pers.dom20, alba.dom.1A_alba.dom20, pers.dom.1A_pers.dom20)

#                        df     AICc
# additive                4 443.0884
# pers.dom.1A_alba.dom20  4 428.1267
# alba.dom.1A_pers.dom20  4 479.5044
# alba.dom.1A_alba.dom20  4 442.2417
# pers.dom.1A_pers.dom20  4 476.4581

BIC(additive, pers.dom.1A_alba.dom20, alba.dom.1A_pers.dom20, alba.dom.1A_alba.dom20, pers.dom.1A_pers.dom20)

#                        df      BIC
# additive                4 450.6578
# pers.dom.1A_alba.dom20  4 435.6961
# alba.dom.1A_pers.dom20  4 487.0738
# alba.dom.1A_alba.dom20  4 449.8111
# pers.dom.1A_pers.dom20  4 484.0275

# pers.dom.1A_alba.dom20 (hypothesized mechnaism) is the best model across all tested

# Now adding the epistasis term 
pers.dom.1A_alba.dom20_epistasis <- lm(head.plumage.value ~ pers.dom.1A + alba.dom20 + pers.dom.1A:alba.dom20, data = data)

AICc(additive, pers.dom.1A_alba.dom20, alba.dom.1A_pers.dom20, alba.dom.1A_alba.dom20, pers.dom.1A_pers.dom20, pers.dom.1A_alba.dom20_epistasis)

#                                  df      AIC
# additive                          4 442.3477
# pers.dom.1A_alba.dom20            4 427.3859
# alba.dom.1A_pers.dom20            4 478.7637
# alba.dom.1A_alba.dom20            4 441.5010
# pers.dom.1A_pers.dom20            4 475.7174
# pers.dom.1A_alba.dom20_epistasis  5 423.3109

#Model with dominance and epistatsis is the best. AIC, BIC and AICc agree with that.

#Comparing the addiditve + epistatsis model to pers.dom.1A_alba.dom20_epistasis
additive_epistasis <- lm(head.plumage.value ~ add1A + add20 + add1A:add20, data = data)

AIC(additive, pers.dom.1A_alba.dom20, alba.dom.1A_pers.dom20, alba.dom.1A_alba.dom20, pers.dom.1A_pers.dom20, pers.dom.1A_alba.dom20_epistasis, additive_epistasis)

#                                  df      AIC
# additive                          4 442.3477
# pers.dom.1A_alba.dom20            4 427.3859
# alba.dom.1A_pers.dom20            4 478.7637
# alba.dom.1A_alba.dom20            4 441.5010
# pers.dom.1A_pers.dom20            4 475.7174
# pers.dom.1A_alba.dom20_epistasis  5 423.3109
# additive_epistasis                5 444.3477
