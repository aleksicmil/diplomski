# Plots

library(VennDiagram)

draw.pairwise.venn(653, 951, 343, 
                   scaled = TRUE, 
                   col = c("lightsalmon", "lightblue2"),
                   fill = c("lightsalmon", "lightblue2"))


draw.pairwise.venn(653, 951, 343, 
                   scaled = TRUE, 
                   col = c("lightsalmon", "lightblue2"), 
                   fill = c("lightsalmon", "lightblue2"))

draw.triple.venn(653, 2747, 951, 
                 99, 147, 343, 55, 
                 scaled = T, 
                 col = c("lightsalmon", "greenyellow", "lightblue2"), 
                 fill = c("lightsalmon", "greenyellow", "lightblue2"))

# install.packages("tidyverse")
library(ggplot2)

pred <- as.data.frame(
  matrix(c(99, 2648, 554, 343, 608, 310),
         nrow = 6, ncol = 1, byrow = T))
pred[ ,2] <- as.factor(rep(c("FRENKI", "PANNZER"), each = 3))
pred[ ,3] <- as.factor(rep(c("TP", "FP", "FN"), 2))

rez <- as.data.frame(
  matrix(c(0.036, 0.151, 0.058, 0.36, 0.525, 0.427),
         nrow = 6, ncol = 1, byrow = T))
rez[ ,2] <- as.factor(rep(c("FRENKI", "PANNZER"), each = 3))
rez[ ,3] <- as.factor(rep(c("Preciznost", "Odziv", "F mera"), 2))

names(rez) <- c("score", "tool", "klasa")
names(pred) <- c("score", "tool", "klasa")

g1 <- ggplot(pred, aes(klasa, score))
g1 <- g1 + geom_col(aes(fill = tool),  
                    position = position_dodge(width = 1)) 

g2 <- ggplot(rez, aes(klasa, score))
g2 <- g2 +
  geom_col(aes(fill = tool), position = "dodge")

cafa <- read.table("Data/trening_cafa_C_ext.txt", header = F, sep = "\t")
table(cafa$V2 == "GO:0005739")
str(levels(cafa$V1))



