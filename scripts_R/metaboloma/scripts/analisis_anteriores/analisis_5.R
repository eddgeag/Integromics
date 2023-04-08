library(reshape2)
library(dplyr)
library(ggplot2)
library(limma)
library(car)
library(Hotelling)
library(factoextra)
library(FactoMineR)
library(ggpubr)
add_classificacion_metaboloma <- function(X, metaboloma) {
  aromatic <- colnames(metaboloma)[c(2, 3, 4, 1, 6, 35, 36, 34)]
  other <-
    colnames(metaboloma)[c(21, 9, 16, 26, 22, 18, 13, 29, 30, 27)]
  aa <- colnames(metaboloma)[c(33, 15, 20, 19, 23, 25)]
  carbo <-
    colnames(metaboloma)[!colnames(metaboloma) %in% c(aromatic, other, aa)]
  
  X$grupo_meta <-
    ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
  X$grupo_meta <-
    ifelse(colnames(metaboloma) %in% other,
           "OTHER",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(colnames(metaboloma) %in% aa,
           "DERIVADOS AA",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(
      colnames(metaboloma) %in% carbo,
      "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
      X$grupo_meta
    )
  
  return(X)
}

# library(MetaboAnalystR)
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/manova_vanvalen.R")

datos <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_none_.rds")

sexo <- datos$comunes$variables_in_bacteria$SEX
obesidad <- datos$comunes$variables_in_bacteria$OBESE
grupo <- datos$comunes$variables_in_bacteria$GROUP

metaboloma <- as.data.frame(t(datos$comunes$metaboloma))

metaboloma.log <- log(metaboloma)

sexo.control <- as.factor(as.character(grupo[grupo!="PCOS"]))
mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
herma <- as.factor(as.character(grupo[grupo!="Female"]))

pmanova <- adonis2(scale(metaboloma.log)~grupo*obesidad,method = "euclidean",permutations = 10000)$`Pr(>F)`[1:3]
females <- adonis2(scale(metaboloma.log)[grupo!="Male",] ~ mujeres,method = "euclidean",permutations = 10000)$`Pr(>F)`[1]
pcosMales<- adonis2(scale(metaboloma.log)[grupo!="Female",] ~ herma,method = "euclidean",permutations = 10000)$`Pr(>F)`[1]
sexo.c <- adonis2(scale(metaboloma.log)[grupo!="PCOS",] ~ sexo.control,method = "euclidean",permutations = 10000)$`Pr(>F)`[1]
p.multi <- p.adjust(c(females,pcosMales,sexo.c),"bonferroni")
names(p.multi) <- c("females","pcosmales","control")

mujeres<- scale(metaboloma.log)[grupo=="Female" ,]
hombres<- scale(metaboloma.log)[grupo=="Male" ,]
Pcos <- scale(metaboloma.log)[grupo=="PCOS" ,]

factor.mujeres <- as.factor(as.character(obesidad[grupo=="Female"]))
factor.males <- as.factor(as.character(obesidad[grupo=="Male"]))
factor.PCOS<- as.factor(as.character(obesidad[grupo=="PCOS"]))


p.mujeres<- adonis2(mujeres ~ factor.mujeres,method = "euclidean",permutations = 10000)$`Pr(>F)`[1]
p.pcos <- adonis2(Pcos~ factor.PCOS,method = "euclidean",permutations = 10000)$`Pr(>F)`[1]
p.males <- adonis2(hombres~ factor.males,method = "euclidean",permutations = 10000)$`Pr(>F)`[1]

p.inter <- p.adjust(c(p.mujeres,p.pcos,p.males),"bonferroni")
names(p.inter) <- c("mujeres","pcos","hombres")

pcx <- prcomp(scale(metaboloma.log))
varianzas <- round(100*(pcx$sdev^2/sum(pcx$sdev^2)),2)
ploteo <- data.frame(PC1=pcx$x[,1],PC2=pcx$x[,2],obesidad=obesidad,grupo=grupo,sexo=sexo)
p1 <- ggplot(ploteo,aes(PC1,PC2,color=obesidad))+geom_point()+xlab(paste("PC1 ",varianzas[1],"%"))+ylab(paste("PC2",varianzas[2],"%"))+facet_grid(~grupo)
p4 <- fviz_screeplot(pcx,top=10)
p5<-ggsave(p1,filename = "./12-2-22/scripts_analysis_reunion/scripts_metaboloma2/PCA.jpeg")


resultados <- data.frame(grupo=pmanova[1],obesidad=pmanova[2],interaccion=pmanova[3])
