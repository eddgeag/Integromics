
library(ALDEx2)
library(dplyr)
library(ggplot2)
library(zCompositions)
library(lmPerm)
library(car)
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/otras_funciones_utiles.R")

datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")

genero <- as.data.frame(t(mean_aldex(datos$comunes$microbiota$genero)))
filo <- as.data.frame(t(mean_aldex(datos$comunes$microbiota$filo)))

grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE
sexo <- datos$comunes$variables_in_bacteria$SEX

funttest <- function(x,categotica) apply(x, 2, function(x) t.test(x ~ categotica)$p.value)
funwilcos <- function(x,categotica) apply(x, 2, function(x) wilcox.test(x ~ categotica)$p.value)
funaov <- function(x,categotica) apply(x, 2, function(x) summary(aov(x ~ categotica))[[1]]$`Pr(>F)`[1])
funpaov <- function(x) apply(x,2,function(x) summary(aovp(x ~ grupo*obesidad,contrasts = list(grupo="contr.sum",obesidad="contr.sum")))[[1]]$`Pr(Prob)`[1:3])
funaov <- function(x) apply(x,2,function(x) Anova(lm(x ~ grupo*obesidad,contrasts=list(grupo="contr.sum",obesidad="contr.sum")),white.adjust=T,type=3)$`Pr(>F)`[2:4])

p.values.perm <- t(funpaov(genero))
colnames(p.values.perm) <- c("grupo","obesidad","interaccion")
p.values.aov <- t(funaov(genero))
colnames(p.values.aov) <- c("grupo","obesidad","interaccion")
p.values.aov.adj <- apply(p.values.aov,2,function(x) p.adjust(x,"BH"))
p.values.paov.adj <- apply(p.values.perm,2,function(x) p.adjust(x,"bonferroni"))


