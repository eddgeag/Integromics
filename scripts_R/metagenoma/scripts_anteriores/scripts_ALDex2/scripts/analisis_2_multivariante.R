
library(ALDEx2)
library(ggplot2)
library(vegan)
library(psych)
library(car)
library(ggstatsplot)
library(reshape2)
library(ggsignif)
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/otras_funciones_utiles.R")
fichero.datos <- "../../datos/preprocesado_05_02_23/novoom.rds"
datos <- readRDS(fichero.datos)
grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE
sexo <- datos$comunes$variables_in_bacteria$SEX

### Empezamos con las diversidades

genero.abs <- datos$comunes$microbiota$genero.abs
filo.abs <- datos$comunes$microbiota$filo.abs


shannon.index <- data.frame(shannon=shannon(t(genero.abs)))
shannon.index$grupo <-grupo
shannon.index$obesidad <- obesidad

mod <- lm(shannon ~ grupo*obesidad,data=shannon.index,contrasts = list(grupo="contr.sum",obesidad="contr.sum"))

shapiro.test(residuals(mod))
Anova(mod=mod,type = 3)

pairwise.t.test(shannon.index$shannon,interaction(grupo,obesidad),p.adjust.method = "BH")

ggplot(shannon.index,aes(x=grupo,y=shannon,fill=obesidad))+geom_boxplot()+

