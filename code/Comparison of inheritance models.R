library(tidyverse)
data <- read.table("compound.genotypes1A.20.txt", header=T)
#this dataset includes only individuals sampled within the hyrbid zone

############### ADDITIVE MODEL ################

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
#adding the epistatic term
additive_epistasis <- lm(head.plumage.value ~ add1A + add20 + add1A:add20, data = data)

############### PARTIAL DOMINANCE MODELS ################

# Partially dominant model where:
# 1A is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model) and
# 20 is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model)
data$pers.dom.1A <- gsub("PP", 1, data$genotype.1A)
data$pers.dom.1A <- gsub("aP", 0.75, data$pers.dom.1A)
data$pers.dom.1A <- gsub("aa", 0, data$pers.dom.1A)
data$pers.dom.1A <- as.numeric(data$pers.dom.1A)

data$alba.dom20 <- gsub("AA", 0, data$genotype.20)
data$alba.dom20 <- gsub("Ap", 0.25, data$alba.dom20)
data$alba.dom20 <- gsub("pp", 1, data$alba.dom20)
data$alba.dom20 <- as.numeric(data$alba.dom20)

pers.dom.1A_alba.dom20 <- lm(head.plumage.value ~ pers.dom.1A + alba.dom20, data = data)
#adding the epistatic term
pers.dom.1A_alba.dom20_epistasis <- lm(head.plumage.value ~ pers.dom.1A + alba.dom20 + pers.dom.1A:alba.dom20, data = data)

# Partially dominant model where 
# 1A is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model) and
# 20 is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model)
data$alba.dom1A <- gsub("PP", 1, data$genotype.1A)
data$alba.dom1A <- gsub("aP", 0.25, data$alba.dom1A)
data$alba.dom1A <- gsub("aa", 0, data$alba.dom1A)
data$alba.dom1A <- as.numeric(data$alba.dom1A)

data$pers.dom20 <- gsub("AA", 0, data$genotype.20)
data$pers.dom20 <- gsub("Ap", 0.75, data$pers.dom20)
data$pers.dom20 <- gsub("pp", 1, data$pers.dom20)
data$pers.dom20 <- as.numeric(data$pers.dom20)

alba.dom.1A_pers.dom20 <- lm(head.plumage.value ~ alba.dom1A + pers.dom20 , data = data)

#adding the epistatic term
alba.dom.1A_pers.dom20_epistasis <- lm(head.plumage.value ~ alba.dom1A + pers.dom20 + alba.dom1A:pers.dom20, data = data)

# Partially dominant model where 
# 1A is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model)
# 20 is partially dominant in alba (heterozygote is 25% more similar to alba compared to addive model)

alba.dom.1A_alba.dom20 <- lm(head.plumage.value ~ alba.dom1A + alba.dom20 , data = data)

#adding the epistatic term
alba.dom.1A_alba.dom20_epistasis <- lm(head.plumage.value ~ alba.dom1A + alba.dom20 + alba.dom1A:alba.dom20, data = data)

# Partially dominant model where 
# 1A is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model)
# 20 is partially dominant in personata (heterozygote is 25% more similar to personata compared to addive model)

pers.dom.1A_pers.dom20 <- lm(head.plumage.value ~ pers.dom.1A + pers.dom20 , data = data)

#adding the epistatic term
pers.dom.1A_pers.dom20_epistasis <- lm(head.plumage.value ~ pers.dom.1A + pers.dom20 + pers.dom.1A:pers.dom20, data = data)

############### FULL DOMINANCE MODELS ################

# 1A is dominant in personata 
# 20 is dominant in alba 
data$pers.full.dom.1A <- gsub("PP", 1, data$genotype.1A)
data$pers.full.dom.1A <- gsub("aP", 1, data$pers.full.dom.1A)
data$pers.full.dom.1A <- gsub("aa", 0, data$pers.full.dom.1A)
data$pers.full.dom.1A <- as.numeric(data$pers.full.dom.1A)

data$alba.full.dom20 <- gsub("AA", 0, data$genotype.20)
data$alba.full.dom20 <- gsub("Ap", 0, data$alba.full.dom20)
data$alba.full.dom20 <- gsub("pp", 1, data$alba.full.dom20)
data$alba.full.dom20 <- as.numeric(data$alba.full.dom20)

pers.full.dom.1A_alba.full.dom20 <- lm(head.plumage.value ~ pers.full.dom.1A + alba.full.dom20, data = data)

#adding the epistatic term
pers.full.dom.1A_alba.full.dom20_epistasis <- lm(head.plumage.value ~ pers.full.dom.1A + alba.full.dom20 + pers.full.dom.1A:alba.full.dom20, data = data)
 
# 1A is dominant in alba 
# 20 is dominant in personata 
data$alba.full.dom1A <- gsub("PP", 1, data$genotype.1A)
data$alba.full.dom1A <- gsub("aP", 0, data$alba.full.dom1A)
data$alba.full.dom1A <- gsub("aa", 0, data$alba.full.dom1A)
data$alba.full.dom1A <- as.numeric(data$alba.full.dom1A)

data$pers.full.dom20 <- gsub("AA", 0, data$genotype.20)
data$pers.full.dom20 <- gsub("Ap", 1, data$pers.full.dom20)
data$pers.full.dom20 <- gsub("pp", 1, data$pers.full.dom20)
data$pers.full.dom20 <- as.numeric(data$pers.full.dom20)

alba.full.dom.1A_pers.full.dom20 <- lm(head.plumage.value ~ alba.full.dom1A + pers.full.dom20 , data = data)

#adding the epistatic term
alba.full.dom.1A_pers.full.dom20_epistasis <- lm(head.plumage.value ~ alba.full.dom1A + pers.full.dom20 + alba.full.dom1A:pers.full.dom20, data = data)


# 1A is dominant in alba 
# 20 is dominant in alba 

alba.full.dom.1A_alba.full.dom20 <- lm(head.plumage.value ~ alba.full.dom1A + alba.full.dom20 , data = data)

#adding the epistatic term
alba.full.dom.1A_alba.full.dom20_epsitasis <- lm(head.plumage.value ~ alba.full.dom1A + alba.full.dom20 + alba.full.dom1A:alba.full.dom20, data = data)


# 1A is dominant in personata
# 20 is dominant in personata 
pers.full.dom.1A_pers.full.dom20 <- lm(head.plumage.value ~ pers.full.dom.1A + pers.full.dom20 , data = data)

#adding the epistatic term
pers.full.dom.1A_pers.full.dom20_epistasis <- lm(head.plumage.value ~ pers.full.dom.1A + pers.full.dom20 + pers.full.dom.1A:pers.full.dom20, data = data)


############### MODEL COMPARISON ################

AIC(additive, 
    additive_epistasis,
    pers.dom.1A_alba.dom20,
    pers.dom.1A_alba.dom20_epistasis,
    alba.dom.1A_pers.dom20,
    alba.dom.1A_pers.dom20_epistasis,
    alba.dom.1A_alba.dom20,
    alba.dom.1A_alba.dom20_epistasis,
    pers.dom.1A_pers.dom20,
    pers.dom.1A_pers.dom20_epistasis,
    pers.full.dom.1A_alba.full.dom20,
    pers.full.dom.1A_alba.full.dom20_epistasis,
    alba.full.dom.1A_pers.full.dom20,
    alba.full.dom.1A_pers.full.dom20_epistasis,
    alba.full.dom.1A_alba.full.dom20,
    alba.full.dom.1A_alba.full.dom20_epsitasis,
    pers.full.dom.1A_pers.full.dom20,
    pers.full.dom.1A_pers.full.dom20_epistasis
    )

AICc(additive, 
    additive_epistasis,
    pers.dom.1A_alba.dom20,
    pers.dom.1A_alba.dom20_epistasis,
    alba.dom.1A_pers.dom20,
    alba.dom.1A_pers.dom20_epistasis,
    alba.dom.1A_alba.dom20,
    alba.dom.1A_alba.dom20_epistasis,
    pers.dom.1A_pers.dom20,
    pers.dom.1A_pers.dom20_epistasis,
    pers.full.dom.1A_alba.full.dom20,
    pers.full.dom.1A_alba.full.dom20_epistasis,
    alba.full.dom.1A_pers.full.dom20,
    alba.full.dom.1A_pers.full.dom20_epistasis,
    alba.full.dom.1A_alba.full.dom20,
    alba.full.dom.1A_alba.full.dom20_epsitasis,
    pers.full.dom.1A_pers.full.dom20,
    pers.full.dom.1A_pers.full.dom20_epistasis
)

BIC(additive, 
    additive_epistasis,
    pers.dom.1A_alba.dom20,
    pers.dom.1A_alba.dom20_epistasis,
    alba.dom.1A_pers.dom20,
    alba.dom.1A_pers.dom20_epistasis,
    alba.dom.1A_alba.dom20,
    alba.dom.1A_alba.dom20_epistasis,
    pers.dom.1A_pers.dom20,
    pers.dom.1A_pers.dom20_epistasis,
    pers.full.dom.1A_alba.full.dom20,
    pers.full.dom.1A_alba.full.dom20_epistasis,
    alba.full.dom.1A_pers.full.dom20,
    alba.full.dom.1A_pers.full.dom20_epistasis,
    alba.full.dom.1A_alba.full.dom20,
    alba.full.dom.1A_alba.full.dom20_epsitasis,
    pers.full.dom.1A_pers.full.dom20,
    pers.full.dom.1A_pers.full.dom20_epistasis
)    

# Summary of the best model  
summary(pers.dom.1A_alba.dom20_epistasis)
