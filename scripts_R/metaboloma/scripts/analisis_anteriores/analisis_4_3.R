library(reshape2)
library(dplyr)
library(ggplot2)
library(limma)
library(car)
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
metaboloma <- metaboloma/rowSums(metaboloma)
metaboloma.log <- log(metaboloma)
#### factores
sexo.control <- as.factor(as.character(grupo[grupo!="PCOS"]))
mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
herma <- as.factor(as.character(grupo[grupo!="Female"]))

mujeres.control <- as.factor(as.character(obesidad[grupo=="Female"]))
pcos.control <- as.factor(as.character(obesidad[grupo=="PCOS"]))
hombres.control <- as.factor(as.character(obesidad[grupo=="Male"]))

sexo.obesidad.general <- as.factor(as.character(sexo[obesidad=="Obese"]))
sexo.NoObesidad.general <- as.factor(as.character(sexo[obesidad=="No Obese"]))
sexo.obesidad.control <- as.factor(as.character(grupo[obesidad=="Obese" & grupo!="PCOS"]))
sexo.NoObesidad.control <- as.factor(as.character(grupo[obesidad=="No Obese" & grupo !="PCOS"]))
mujeres.obesidad <- as.factor(as.character(grupo[obesidad=="Obese" & grupo!="Male"]))
mujeres.NoObesidad <- as.factor(as.character(grupo[obesidad=="No Obese" & grupo!="Male"]))
herma.Obesidad <- as.factor(as.character(grupo[obesidad=="Obese" & grupo!="Female"]))
herma.NoObesidad <- as.factor(as.character(grupo[obesidad=="No Obese" & grupo!="Female"]))

### Metabolomas para el test

### grupo
metaboloma.grupo.control <- metaboloma.log[grupo!="PCOS",]
metaboloma.grupo.mujeres <- metaboloma.log[grupo!="Male",]
metaboloma.grupo.herma <- metaboloma.log[grupo!="Female",]
### obesos y no obesos para migrar variables grupo
metaboloma.grupo.control.obesos <- as.data.frame(metaboloma.log[grupo!="PCOS" & obesidad=="Obese",])
metaboloma.grupo.control.NoObesos <- as.data.frame(metaboloma.log[grupo!="PCOS" & obesidad=="No Obese",])
metaboloma.grupo.mujeres.obesos <- as.data.frame(metaboloma.log[grupo!="Male" & obesidad=="Obese",])
metaboloma.grupo.mujeres.NoObesos <- as.data.frame(metaboloma.log[grupo!="Male" & obesidad=="No Obese",])
metaboloma.grupo.herma.obesos <- as.data.frame(metaboloma.log[grupo!="Female" & obesidad=="Obese",])
metaboloma.grupo.herma.NoObesos <- as.data.frame(metaboloma.log[grupo!="Female" & obesidad=="No Obese",])
### varaibles grupo para mirar obesidad
metaboloma.obesidad.mujeres.control <- as.data.frame(metaboloma.log[grupo=="Female",])
metaboloma.obesidad.pcos.control <- as.data.frame(metaboloma.log[grupo=="PCOS",])
metaboloma.obesidad.hombres.control <- as.data.frame(metaboloma.log[grupo=="Male",])

### metabolomas para las medias

metaboloma.grupo.control.females.obesos <- metaboloma[grupo=="Female" & obesidad=="Obese",]
metaboloma.grupo.control.females.NoObesos <- metaboloma[grupo=="Female" & obesidad=="No Obese",]
metaboloma.grupo.control.males.obesos <- metaboloma[grupo=="Male" & obesidad=="Obese",]
metaboloma.grupo.control.males.NoObesos <- metaboloma[grupo=="Male" & obesidad=="No Obese",]
metaboloma.grupo.pcos.Obesos <- metaboloma[grupo=="PCOS" & obesidad=="Obese",]
metaboloma.grupo.pcos.NoObesos <- metaboloma[grupo=="PCOS" & obesidad=="No Obese",]

metaboloma.obesidad.pcos <- metaboloma[grupo=="PCOS",]
metaboloma.obesidad.mujeres <- metaboloma[grupo=="Female",]
metaboloma.obesidad.hombres <- metaboloma[grupo=="Male",]

metaboloma.obesos.control.hombres <- metaboloma[grupo=="Male" & obesidad=="Obese",]
metaboloma.NoObesos.control.hombres <- metaboloma[grupo=="Male" & obesidad=="No Obese",]
metaboloma.obesos.control.females <- metaboloma[grupo=="Female" & obesidad=="Obese",]
metaboloma.Nobesos.control.females <- metaboloma[grupo=="Female" & obesidad=="No Obese",]
metaboloma.obesos.control.pcos <- metaboloma[grupo=="PCOS" & obesidad=="Obese",]
metaboloma.NoObesos.control.pcos <- metaboloma[grupo=="PCOS" & obesidad=="No Obese",]



fun_categorica <- function(x,categorica) apply(x,2,function(y) t.test(y~categorica,var.equal=F)$p.value)

### grupo|obesidad para mirar por nivel de grupo la obesidad
### MUJERES
midf.obesidad.control.mujeres <- data.frame(p.value=fun_categorica(metaboloma.obesidad.mujeres.control,mujeres.control))
midf.obesidad.control.mujeres$p.adj <- p.adjust(midf.obesidad.control.mujeres$p.value,"BH")
midf.obesidad.control.mujeres <- add_classificacion_metaboloma(midf.obesidad.control.mujeres,metaboloma)
midf.obesidad.control.mujeres$Obese <- colMeans(metaboloma.grupo.control.females.obesos)
midf.obesidad.control.mujeres$NoObese <- colMeans(metaboloma.grupo.control.females.NoObesos)
midf.obesidad.control.mujeres$interpretacion <- ifelse(midf.obesidad.control.mujeres$Obese>midf.obesidad.control.mujeres$NoObese,"Aumenta","Disminuye")
midf.obesidad.control.mujeres$sig <- ifelse(midf.obesidad.control.mujeres$p.value<0.05,"*","")
midf.obesidad.control.mujeres$sig.adj <- ifelse(midf.obesidad.control.mujeres$p.adj<0.05,"**","")
### PCOS
midf.obesidad.control.PCOS <- data.frame(p.value=fun_categorica(metaboloma.obesidad.pcos.control,pcos.control))
midf.obesidad.control.PCOS$p.adj <- p.adjust(midf.obesidad.control.PCOS$p.value,"BH")
midf.obesidad.control.PCOS <- add_classificacion_metaboloma(midf.obesidad.control.PCOS,metaboloma)
midf.obesidad.control.PCOS$Obese <- colMeans(metaboloma.grupo.control.females.obesos)
midf.obesidad.control.PCOS$NoObese <- colMeans(metaboloma.grupo.control.females.NoObesos)
midf.obesidad.control.PCOS$interpretacion <- ifelse(midf.obesidad.control.PCOS$Obese>midf.obesidad.control.PCOS$NoObese,"Aumenta","Disminuye")
midf.obesidad.control.PCOS$sig <- ifelse(midf.obesidad.control.PCOS$p.value<0.05,"*","")
midf.obesidad.control.PCOS$sig.adj <- ifelse(midf.obesidad.control.PCOS$p.adj<0.05,"**","")

### HOMBRES
midf.obesidad.control.hombres <- data.frame(p.value=fun_categorica(metaboloma.obesidad.hombres.control,hombres.control))
midf.obesidad.control.hombres$p.adj <- p.adjust(midf.obesidad.control.hombres$p.value,"BH")
midf.obesidad.control.hombres <- add_classificacion_metaboloma(midf.obesidad.control.hombres,metaboloma)
midf.obesidad.control.hombres$Obese <- colMeans(metaboloma.grupo.control.males.obesos)
midf.obesidad.control.hombres$NoObese <- colMeans(metaboloma.grupo.control.males.NoObesos)
midf.obesidad.control.hombres$interpretacion <- ifelse(midf.obesidad.control.hombres$Obese>midf.obesidad.control.hombres$NoObese,"Aumenta","Disminuye")
midf.obesidad.control.hombres$sig <- ifelse(midf.obesidad.control.hombres$p.value<0.05,"*","")
midf.obesidad.control.hombres$sig.adj <- ifelse(midf.obesidad.control.hombres$p.adj<0.05,"**","")


### Ahora para ver por obesidad los niveles del grupo
## Sexo obesidad
midf.grupo.sexo.control.obesidad <- data.frame(p.value=fun_categorica(metaboloma.grupo.control.obesos,sexo.obesidad.control))
midf.grupo.sexo.control.obesidad$p.adj <- p.adjust(midf.grupo.sexo.control.obesidad$p.value,"BH")
midf.grupo.sexo.control.obesidad <- add_classificacion_metaboloma(midf.grupo.sexo.control.obesidad,metaboloma)
midf.grupo.sexo.control.obesidad$Male <- colMeans(metaboloma.obesos.control.hombres)
midf.grupo.sexo.control.obesidad$Female <- colMeans(metaboloma.obesos.control.females)
midf.grupo.sexo.control.obesidad$interpretacion <- ifelse(midf.grupo.sexo.control.obesidad$Male>midf.grupo.sexo.control.obesidad$Female,"Aumenta","Disminuye")
midf.grupo.sexo.control.obesidad$sig <- ifelse(midf.grupo.sexo.control.obesidad$p.value<0.05,"*","")
midf.grupo.sexo.control.obesidad$sig.adj <- ifelse(midf.grupo.sexo.control.obesidad$p.value<0.05,"**","")
## Seixo No Obesidad
midf.grupo.sexo.control.Noobesidad <- data.frame(p.value=fun_categorica(metaboloma.grupo.control.NoObesos,sexo.NoObesidad.control))
midf.grupo.sexo.control.Noobesidad$p.adj <- p.adjust(midf.grupo.sexo.control.Noobesidad$p.value,"BH")
midf.grupo.sexo.control.Noobesidad <- add_classificacion_metaboloma(midf.grupo.sexo.control.Noobesidad,metaboloma)
midf.grupo.sexo.control.Noobesidad$Male <- colMeans(metaboloma.NoObesos.control.hombres)
midf.grupo.sexo.control.Noobesidad$Female <- colMeans(metaboloma.Nobesos.control.females)
midf.grupo.sexo.control.Noobesidad$interpretacion <- ifelse(midf.grupo.sexo.control.Noobesidad$Male>midf.grupo.sexo.control.Noobesidad$Female,"Aumenta","Disminuye")
midf.grupo.sexo.control.Noobesidad$sig <- ifelse(midf.grupo.sexo.control.Noobesidad$p.value<0.05,"*","")
midf.grupo.sexo.control.Noobesidad$sig.adj <- ifelse(midf.grupo.sexo.control.Noobesidad$p.value<0.05,"**","")

### Mujeres Obesidad

midf.grupo.females.control.obesidad <- data.frame(p.value=fun_categorica(metaboloma.grupo.mujeres.obesos,mujeres.obesidad))
midf.grupo.females.control.obesidad$p.adj <- p.adjust(midf.grupo.females.control.obesidad$p.value,"BH")
midf.grupo.females.control.obesidad <- add_classificacion_metaboloma(midf.grupo.females.control.obesidad,metaboloma)
midf.grupo.females.control.obesidad$PCOS <- colMeans(metaboloma.grupo.pcos.Obesos)
midf.grupo.females.control.obesidad$Female <- colMeans(metaboloma.obesos.control.females)
midf.grupo.females.control.obesidad$interpretacion <- ifelse(midf.grupo.females.control.obesidad$PCOS>midf.grupo.females.control.obesidad$Female,"Aumenta","Disminuye")
midf.grupo.females.control.obesidad$sig <- ifelse(midf.grupo.females.control.obesidad$p.value,"*","")
midf.grupo.females.control.obesidad$sig.adj <- ifelse(midf.grupo.females.control.obesidad$p.adj,"**","")
### mujeres No obesidad 

midf.grupo.females.control.Noobesidad <- data.frame(p.value=fun_categorica(metaboloma.grupo.mujeres.NoObesos,mujeres.NoObesidad))
midf.grupo.females.control.Noobesidad$p.adj <- p.adjust(midf.grupo.females.control.Noobesidad$p.value,"BH")
midf.grupo.females.control.Noobesidad <- add_classificacion_metaboloma(midf.grupo.females.control.Noobesidad,metaboloma)
midf.grupo.females.control.Noobesidad$PCOS <- colMeans(metaboloma.grupo.pcos.NoObesos)
midf.grupo.females.control.Noobesidad$Female <- colMeans(metaboloma.Nobesos.control.females)
midf.grupo.females.control.Noobesidad$interpretacion <- ifelse(midf.grupo.females.control.Noobesidad$PCOS>midf.grupo.females.control.Noobesidad$Female,"Aumenta","Disminuye")
midf.grupo.females.control.Noobesidad$sig <- ifelse(midf.grupo.females.control.Noobesidad$p.value,"*","")
midf.grupo.females.control.Noobesidad$sig.adj <- ifelse(midf.grupo.females.control.Noobesidad$p.adj,"**","")

### hombres obesidad

midf.grupo.males.control.obesidad <- data.frame(p.value=fun_categorica(metaboloma.grupo.herma.obesos,herma.Obesidad))
midf.grupo.males.control.obesidad$p.adj <- p.adjust(midf.grupo.males.control.obesidad$p.value,"BH")
midf.grupo.males.control.obesidad <- add_classificacion_metaboloma(midf.grupo.males.control.obesidad,metaboloma)
midf.grupo.males.control.obesidad$PCOS <- colMeans(metaboloma.grupo.pcos.Obesos)
midf.grupo.males.control.obesidad$Male <- colMeans(metaboloma.grupo.control.males.obesos)
midf.grupo.males.control.obesidad$interpretacion <- ifelse(midf.grupo.males.control.obesidad$PCOS>midf.grupo.males.control.obesidad$Male,"Aumenta","Disminuye")
midf.grupo.males.control.obesidad$sig <- ifelse(midf.grupo.males.control.obesidad$p.value<0.05,"*","")
midf.grupo.males.control.obesidad$sig.adj <- ifelse(midf.grupo.males.control.obesidad$p.adj<0.05,"**","")

midf.grupo.males.control.Noobesidad <- data.frame(p.value=fun_categorica(metaboloma.grupo.herma.NoObesos,herma.NoObesidad))
midf.grupo.males.control.Noobesidad$p.adj <- p.adjust(midf.grupo.males.control.Noobesidad$p.value,"BH")
midf.grupo.males.control.Noobesidad <- add_classificacion_metaboloma(midf.grupo.males.control.Noobesidad,metaboloma)
midf.grupo.males.control.Noobesidad$PCOS <- colMeans(metaboloma.grupo.pcos.NoObesos)
midf.grupo.males.control.Noobesidad$Male <- colMeans(metaboloma.grupo.control.males.NoObesos)
midf.grupo.males.control.Noobesidad$interpretacion <- ifelse(midf.grupo.males.control.Noobesidad$PCOS>midf.grupo.males.control.Noobesidad$Male,"Aumenta","Disminuye")
midf.grupo.males.control.Noobesidad$sig <- ifelse(midf.grupo.males.control.Noobesidad$p.value<0.05,"*","")
midf.grupo.males.control.Noobesidad$sig.adj <- ifelse(midf.grupo.males.control.Noobesidad$p.adj<0.05,"**","")


### grupo - mujeres
midf.grupo.females <- data.frame(p.value=fun_categorica(metaboloma.grupo.mujeres,mujeres))
midf.grupo.females$p.adj <- p.adjust(midf.grupo.females$p.value,"BH")
midf.grupo.females <- add_classificacion_metaboloma(midf.grupo.females,metaboloma )
midf.grupo.females$PCOS <- colMeans(metaboloma[grupo=="PCOS",])
midf.grupo.females$Female <- colMeans(metaboloma[grupo=="Female",])
midf.grupo.females$interpretacion <- ifelse(midf.grupo.females$PCOS>midf.grupo.females$Female,"Aumenta","Disminuye")
midf.grupo.females$sig <- ifelse(midf.grupo.females$p.value<0.05,"*","")
midf.grupo.females$sig.adj <- ifelse(midf.grupo.females$p.adj<0.05,"**","")

### grupo -control
midf.grupo.control <- data.frame(p.value=fun_categorica(metaboloma.grupo.control,sexo.control))
midf.grupo.control$p.adj <- p.adjust(midf.grupo.control$p.value,"BH")
midf.grupo.control <- add_classificacion_metaboloma(midf.grupo.control,metaboloma )
midf.grupo.control$Male <- colMeans(metaboloma[grupo=="Male",])
midf.grupo.control$Female <- colMeans(metaboloma[grupo=="Female",])
midf.grupo.control$interpretacion <- ifelse(midf.grupo.control$Male>midf.grupo.control$Female,"Aumenta","Disminuye")
midf.grupo.control$sig <- ifelse(midf.grupo.control$p.value<0.05,"*","")
midf.grupo.control$sig.adj <- ifelse(midf.grupo.control$p.adj<0.05,"**","")

### gruoi -- herma

midf.grupo.herma <- data.frame(p.value=fun_categorica(metaboloma.grupo.herma,herma))
midf.grupo.herma$p.adj <- p.adjust(midf.grupo.herma$p.value,"BH")
midf.grupo.herma <- add_classificacion_metaboloma(midf.grupo.herma,metaboloma )
midf.grupo.herma$Male <- colMeans(metaboloma[grupo=="Male",])
midf.grupo.herma$PCOS <- colMeans(metaboloma[grupo=="PCOS",])
midf.grupo.herma$interpretacion <- ifelse(midf.grupo.herma$Male>midf.grupo.herma$PCOS,"Aumenta","Disminuye")
midf.grupo.herma$sig <- ifelse(midf.grupo.herma$p.value<0.05,"*","")
midf.grupo.herma$sig.adj <- ifelse(midf.grupo.herma$p.adj<0.05,"**","")

#### obesidad

obesidad.analisis <- data.frame(p.value=fun_categorica(metaboloma,obesidad))
obesidad.analisis$p.adj <- p.adjust(obesidad.analisis$p.value,"BH")
obesidad.analisis <- add_classificacion_metaboloma(obesidad.analisis,metaboloma)
obesidad.analisis$Obese <- colMeans(metaboloma[obesidad=="Obese",])
obesidad.analisis$NoObese <- colMeans(metaboloma[obesidad=="No Obese",])
obesidad.analisis$interpretacion <- ifelse(obesidad.analisis$Obese>obesidad.analisis$NoObese,"Aumenta","Disminuye")
obesidad.analisis$sig <- ifelse(obesidad.analisis$p.value<0.05,"*","")
obesidad.analisis$sig.adj <- ifelse(obesidad.analisis$p.adj<0.05,"**","")

### Sexo 
sexo.analisis <- data.frame(p.value=fun_categorica(metaboloma,sexo))
sexo.analisis$p.adj <- p.adjust(sexo.analisis$p.value,"BH")
sexo.analisis <- add_classificacion_metaboloma(sexo.analisis,metaboloma)
sexo.analisis$Obese <- colMeans(metaboloma[sexo=="Obese",])
sexo.analisis$NoObese <- colMeans(metaboloma[sexo=="No Obese",])
sexo.analisis$interpretacion <- ifelse(sexo.analisis$Obese>sexo.analisis$NoObese,"Aumenta","Disminuye")
sexo.analisis$sig <- ifelse(sexo.analisis$p.value<0.05,"*","")
sexo.analisis$sig.adj <- ifelse(sexo.analisis$p.adj<0.05,"**","")



### anova
mi.anova <- function(x) Anova(lm(x~grupo*obesidad),icontrasts=c("contr.sum","contr,sum"),white.adjust = T,type=3)$`Pr(>F)`[2:4]

res <- t(apply(metaboloma.log, 2, mi.anova))
colnames(res) <- c("grupo","obesidad","interaccion")
res <- as.data.frame(res)

tabla_pvalores <- data.frame(aov.obesidad=res$obesidad,
                             aov.grupo=res$grupo,
                             aov.interaccion=res$interaccion,
                             obesidad=obesidad.analisis$p.value,
                             Male.vs.Fmeale = midf.grupo.control$p.value,
                             Pcos.vs.Female = midf.grupo.females$p.value,
                             Pcos.vs.Male = midf.grupo.herma$p.value,
                             MalesObesity = midf.obesidad.control.hombres$p.value,
                             FemalesObesity = midf.obesidad.control.mujeres$p.value,
                             PCOS.Obesity = midf.obesidad.control.PCOS$p.value,
                             Obesity.Sex = midf.grupo.sexo.control.obesidad$p.value,
                             Obesity.Females = midf.grupo.females.control.obesidad$p.value,
                             Obesity.Herma = midf.grupo.males.control.obesidad$p.value,
                             NoObesity.Sex = midf.grupo.sexo.control.Noobesidad$p.value,
                             NoObesity.Females = midf.grupo.females.control.Noobesidad$p.value,
                             NoObesity.Herma = midf.grupo.males.control.Noobesidad$p.value)

tabla_p.adj <- tabla_pvalores
tabla_p.adj[,5:7] <- apply(tabla_p.adj[,5:7],1,function(x) p.adjust(x,"BH"))
tabla_p.adj[,8:16] <- apply(tabla_p.adj[,8:16],1,function(x) p.adjust(x,"BH"))

tabla_p.adj$obesidad<- ifelse(tabla_p.adj$aov.obesidad>0.05,NA,tabla_p.adj$obesidad)
tabla_p.adj$Male.vs.Fmeale <- ifelse(tabla_p.adj$aov.grupo>0.05,NA,tabla_p.adj$Male.vs.Fmeale)
tabla_p.adj$Pcos.vs.Female <- ifelse(tabla_p.adj$aov.grupo>0.05,NA,tabla_p.adj$Pcos.vs.Female)
tabla_p.adj$Pcos.vs.Male <- ifelse(tabla_p.adj$aov.grupo>0.05,NA,tabla_p.adj$Pcos.vs.Male)

tabla_p.adj$MalesObesity <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$MalesObesity)
tabla_p.adj$FemalesObesity <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$FemalesObesity)
tabla_p.adj$PCOS.Obesity <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$PCOS.Obesity)
tabla_p.adj$Obesity.Sex <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$Obesity.Sex)
tabla_p.adj$Obesity.Females <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$Obesity.Females)
tabla_p.adj$Obesity.Herma <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$Obesity.Herma)
tabla_p.adj$NoObesity.Sex<- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$NoObesity.Sex)
tabla_p.adj$NoObesity.Females <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$NoObesity.Females)
tabla_p.adj$NoObesity.Herma <- ifelse(tabla_p.adj$aov.interaccion>0.05,NA,tabla_p.adj$NoObesity.Herma)

tabla_p.adj$obesidad<- ifelse(tabla_p.adj$obesidad>0.05 & !is.na(tabla_p.adj$obesidad),NA,tabla_p.adj$obesidad)
tabla_p.adj$Male.vs.Fmeale <- ifelse(tabla_p.adj$Male.vs.Fmeale >0.05 & !is.na(tabla_p.adj$Male.vs.Fmeale),NA,tabla_p.adj$Male.vs.Fmeale)
tabla_p.adj$Pcos.vs.Female <- ifelse(tabla_p.adj$Pcos.vs.Female>0.05 & !is.na(tabla_p.adj$Pcos.vs.Female) ,NA,tabla_p.adj$Pcos.vs.Female)
tabla_p.adj$Pcos.vs.Male <- ifelse(tabla_p.adj$Pcos.vs.Male>0.05 & !is.na(tabla_p.adj$Pcos.vs.Male)  ,NA,tabla_p.adj$Pcos.vs.Male)

tabla_p.adj$MalesObesity <- ifelse(tabla_p.adj$MalesObesity >0.05 & !is.na(tabla_p.adj$MalesObesity),NA,tabla_p.adj$MalesObesity)
tabla_p.adj$FemalesObesity <- ifelse(tabla_p.adj$FemalesObesity>0.05 & !is.na(tabla_p.adj$FemalesObesity),NA,tabla_p.adj$FemalesObesity)
tabla_p.adj$PCOS.Obesity <- ifelse(tabla_p.adj$PCOS.Obesity>0.05 & !is.na(tabla_p.adj$PCOS.Obesity),NA,tabla_p.adj$PCOS.Obesity)
tabla_p.adj$Obesity.Sex <- ifelse(tabla_p.adj$Obesity.Sex>0.05 & !is.na(tabla_p.adj$Obesity.Sex),NA,tabla_p.adj$Obesity.Sex)
tabla_p.adj$Obesity.Females <- ifelse(tabla_p.adj$Obesity.Females>0.05 & !is.na(tabla_p.adj$Obesity.Females),NA,tabla_p.adj$Obesity.Females)
tabla_p.adj$Obesity.Herma <- ifelse(tabla_p.adj$Obesity.Herma>0.05 & !is.na(tabla_p.adj$Obesity.Herma),NA,tabla_p.adj$Obesity.Herma)
tabla_p.adj$NoObesity.Sex<- ifelse(tabla_p.adj$NoObesity.Sex>0.05 & !is.na(tabla_p.adj$NoObesity.Sex),NA,tabla_p.adj$NoObesity.Sex)
tabla_p.adj$NoObesity.Females <- ifelse(tabla_p.adj$NoObesity.Females>0.05 & !is.na(tabla_p.adj$NoObesity.Females),NA,tabla_p.adj$NoObesity.Females)
tabla_p.adj$NoObesity.Herma <- ifelse(tabla_p.adj$NoObesity.Herma>0.05 & !is.na(tabla_p.adj$NoObesity.Herma),NA,tabla_p.adj$NoObesity.Herma)

rownames(tabla_p.adj) <- colnames(metaboloma)
tabla_p.adj <- add_classificacion_metaboloma(tabla_p.adj,metaboloma)

interpretacion <- tabla_p.adj

interpretacion$obesidad <- ifelse(tabla_p.adj$aov.obesidad<0.05,obesidad.analisis$interpretacion,NA)
interpretacion$Male.vs.Fmeale <- ifelse(tabla_p.adj$Male.vs.Fmeale<0.05 & !is.na(tabla_p.adj$Male.vs.Fmeale),midf.grupo.control$interpretacion,NA)
interpretacion$Pcos.vs.Female<- ifelse(tabla_p.adj$Pcos.vs.Female<0.05 & !is.na(tabla_p.adj$Pcos.vs.Female),midf.grupo.females$interpretacion,NA)
interpretacion$Pcos.vs.Male<- ifelse(tabla_p.adj$Pcos.vs.Male<0.05 & !is.na(tabla_p.adj$Pcos.vs.Male),midf.grupo.herma$interpretacion,NA)
interpretacion$MalesObesity <- ifelse(tabla_p.adj$MalesObesity < 0.05 & !is.na(tabla_p.adj$MalesObesity),midf.obesidad.control.hombres$interpretacion,NA)
interpretacion$FemalesObesity <- ifelse(tabla_p.adj$FemalesObesity <0.05 & !is.na(tabla_p.adj$FemalesObesity),midf.obesidad.control.mujeres$interpretacion,NA)
interpretacion$PCOS.Obesity <- ifelse(tabla_p.adj$PCOS.Obesity <0.05 & !is.na(tabla_p.adj$PCOS.Obesity),midf.obesidad.control.PCOS$interpretacion,NA)
interpretacion$Obesity.Sex <- ifelse(tabla_p.adj$Obesity.Sex <0.05 & !is.na(tabla_p.adj$Obesity.Sex),midf.grupo.sexo.control.obesidad$interpretacion,NA)
interpretacion$Obesity.Females <- ifelse(tabla_p.adj$Obesity.Females <0.05 & !is.na(tabla_p.adj$Obesity.Females),midf.grupo.females.control.obesidad$interpretacion,NA)
interpretacion$Obesity.Herma <- ifelse(tabla_p.adj$Obesity.Herm< 0.05 & !is.na(tabla_p.adj$Obesity.Herma),midf.grupo.males.control.obesidad$interpretacion,NA)
interpretacion$NoObesity.Sex <- ifelse(tabla_p.adj$NoObesity.Sex < 0.05 & !is.na(tabla_p.adj$NoObesity.Sex),midf.grupo.sexo.control.Noobesidad$interpretacion,NA)
interpretacion$NoObesity.Females <- ifelse(tabla_p.adj$NoObesity.Females < 0.05 & !is.na(tabla_p.adj$NoObesity.Females),midf.grupo.sexo.control.Noobesidad$interpretacion,NA)
interpretacion$NoObesity.Herma <- ifelse(tabla_p.adj$NoObesity.Herma <0.05 & !is.na(tabla_p.adj$NoObesity.Herma),midf.grupo.males.control.Noobesidad$interpretacion,NA)

write.csv(interpretacion,"./12-2-22/scripts_analysis_reunion/scripts_metaboloma2/Resultados_univariante_metaboloma/interpretacion_comunes_rowsums.csv")
