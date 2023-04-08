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

sexo <- datos$totales$general_data$SEX
obesidad <- datos$totales$general_data$OBESE
grupo <- datos$totales$general_data$GROUP

comunes <- scale(t(datos$totales$metaboloma))

p.valor.sexo <- apply(metaboloma.total,2,function(x) t.test(x~sexo,var.equal=F)$p.value)
p.valor.sexo.adj <- p.adjust(p.valor.sexo,"BH")

res_sexo <- data.frame(p.valor.sexo=p.valor.sexo,p.valor.sexo.adj=p.valor.sexo.adj)


res_sexo<- add_classificacion_metaboloma(res_sexo,metaboloma.total)
res_sexo$Males <- apply((metaboloma.total[sexo=="Male",]),2,mean)
res_sexo$Females<- apply((metaboloma.total[sexo=="Female",]),2,mean)

res_sexo <- res_sexo[order(res_sexo$p.valor.sexo.adj),]

#### obesidad en hombres

p.valor.hombres <- apply(metaboloma.total[sexo=="Male",],2,function(x) t.test(x~as.factor(as.character(obesidad[sexo=="Male"])),var.equal=F)$p.value)
p.valor.hombres.adj <- p.adjust(p.valor.hombres,"BH")

res_hombres <- data.frame(p.valor.hombres=p.valor.hombres,p.valor.hombres.adj=p.valor.hombres.adj)


res_hombres<- add_classificacion_metaboloma(res_hombres,metaboloma.total)
res_hombres$medias_obesos <- apply((metaboloma.total[obesidad[sexo=="Male"]=="Obese",]),2,mean)
res_hombres$medias_NOobesos <- apply((metaboloma.total[obesidad[sexo=="Male"]=="No Obese",]),2,mean)

res_hombres <- res_hombres[order(res_hombres$p.valor.hombres.adj),]

#### obesidad en mujeres

p.valor.mujeres <- apply(metaboloma.total[sexo=="Female",],2,function(x) t.test(x~as.factor(as.character(obesidad[sexo=="Female"])),var.equal=F)$p.value)
p.valor.mujeres.adj <- p.adjust(p.valor.mujeres,"BH")

res_mujeres <- data.frame(p.valor.mujeres=p.valor.mujeres,p.valor.mujeres.adj=p.valor.mujeres.adj)


res_mujeres<- add_classificacion_metaboloma(res_mujeres,metaboloma.total)
res_mujeres$medias_obesos <- apply((metaboloma.total[obesidad[sexo=="Female"]=="Obese",]),2,mean)
res_mujeres$medias_NOobesos <- apply((metaboloma.total[obesidad[sexo=="Female"]=="No Obese",]),2,mean)

res_mujeres <- res_mujeres[order(res_mujeres$p.valor.mujeres.adj),]


#### obesidad en PCOS

p.valor.PCOS <- apply(metaboloma.total[grupo=="PCOS",],2,function(x) t.test(x~as.factor(as.character(obesidad[grupo=="PCOS"])),var.equal=F)$p.value)
p.valor.PCOS.adj <- p.adjust(p.valor.PCOS,"BH")

res_PCOS <- data.frame(p.valor.PCOS=p.valor.PCOS,p.valor.PCOS.adj=p.valor.PCOS.adj)


res_PCOS<- add_classificacion_metaboloma(res_PCOS,metaboloma.total)
res_PCOS$medias_obesos <- apply((metaboloma.total[obesidad[grupo=="PCOS"]=="Obese",]),2,mean)
res_PCOS$medias_NOobesos <- apply((metaboloma.total[obesidad[grupo=="PCOS"]=="No Obese",]),2,mean)

res_PCOS <- res_PCOS[order(res_PCOS$p.valor.PCOS.adj),]


#### obesidad en mujeres control

p.valor.mujeres.control <- apply(metaboloma.total[grupo=="Female",],2,function(x) t.test(x~as.factor(as.character(obesidad[grupo=="Female"])),var.equal=F)$p.value)
p.valor.mujeres.control.adj <- p.adjust(p.valor.mujeres.control,"BH")

res_mujeres.control <- data.frame(p.valor.mujeres.control=p.valor.mujeres.control,p.valor.mujeres.control.adj=p.valor.mujeres.control.adj)


res_mujeres.control<- add_classificacion_metaboloma(res_mujeres.control,metaboloma.total)
res_mujeres.control$medias_obesos <- apply((metaboloma.total[obesidad[grupo=="Female"]=="Obese",]),2,mean)
res_mujeres.control$medias_NOobesos <- apply((metaboloma.total[obesidad[grupo=="Female"]=="No Obese",]),2,mean)

res_mujeres.control <- res_mujeres.control[order(res_mujeres.control$p.valor.mujeres.control.adj),]


#### VEMOS OBESIDAD EN EL GRUPO control

p.valor.sujetos.obesos.control <- apply(metaboloma.total[obesidad=="Obese" & grupo!="PCOS",],2,function(x) t.test(x~as.factor(as.character(grupo[obesidad=="Obese" & grupo!="PCOS"])),var.equal=F)$p.value)
p.valor.sujetos.obesos.control.adj <- p.adjust(p.valor.sujetos.obesos.control,"BH")

res_sujetos.obesos.control <- data.frame(p.valor.sujetos.obesos.control=p.valor.sujetos.obesos.control,p.valor.sujetos.obesos.control.adj=p.valor.sujetos.obesos.control.adj)


res_sujetos.obesos.control<- add_classificacion_metaboloma(res_sujetos.obesos.control,metaboloma.total)
res_sujetos.obesos.control$medias_Females <- apply((metaboloma.total[grupo[obesidad=="Obese" & grupo!="PCOS"]=="Female",]),2,mean)
res_sujetos.obesos.control$medias_Males <- apply((metaboloma.total[grupo[obesidad=="Obese" & grupo!="PCOS"]=="Male",]),2,mean)

res_sujetos.obesos.control <- res_sujetos.obesos.control[order(res_sujetos.obesos.control$p.valor.sujetos.obesos.control.adj),]

### VEMOS obesidad en las mujeres en genral

p.valor.mujeres.general <- apply(metaboloma.total[obesidad=="Obese" & grupo!="Male",],2,function(x) t.test(x~as.factor(as.character(grupo[obesidad=="Obese" & grupo!="Male"])),var.equal=F)$p.value)
p.valor.mujeres.general.adj <- p.adjust(p.valor.mujeres.general,"BH")

res_mujeres.general <- data.frame(p.valor.mujeres.general=p.valor.mujeres.general,p.valor.mujeres.general.adj=p.valor.mujeres.general.adj)


res_mujeres.general<- add_classificacion_metaboloma(res_mujeres.general,metaboloma.total)
res_mujeres.general$medias_PCOS <- apply((metaboloma.total[grupo[obesidad=="Obese" & grupo!="Male"]=="PCOS",]),2,mean)
res_mujeres.general$medias_Females <- apply((metaboloma.total[grupo[obesidad=="Obese" & grupo!="Male"]=="Female",]),2,mean)

res_mujeres.general <- res_mujeres.general[order(res_mujeres.general$p.valor.mujeres.general.adj),]


### VEMOS obesidad en las hombres vs PCOS en genral

p.valor.hombres.pcos <- apply(metaboloma.total[obesidad=="Obese" & grupo!="Female",],2,function(x) t.test(x~as.factor(as.character(grupo[obesidad=="Obese" & grupo!="Female"])),var.equal=F)$p.value)
p.valor.hombres.pcos.adj <- p.adjust(p.valor.hombres.pcos,"BH")

res_hombres.pcos <- data.frame(p.valor.hombres.pcos=p.valor.hombres.pcos,p.valor.hombres.pcos.adj=p.valor.hombres.pcos.adj)


res_hombres.pcos<- add_classificacion_metaboloma(res_hombres.pcos,metaboloma.total)
res_hombres.pcos$medias_PCOS <- apply((metaboloma.total[grupo[obesidad=="Obese" & grupo!="Female"]=="PCOS",]),2,mean)
res_hombres.pcos$medias_Females <- apply((metaboloma.total[grupo[obesidad=="Obese" & grupo!="Female"]=="Male",]),2,mean)

res_hombres.pcos <- res_hombres.pcos[order(res_hombres.pcos$p.valor.hombres.pcos.adj),]



#### VEMOS NO OBESIDAD EN EL GRUPO control

p.valor.sujetos.obesos.control <- apply(metaboloma.total[obesidad=="No Obese" & grupo!="PCOS",],2,function(x) t.test(x~as.factor(as.character(grupo[obesidad=="No Obese" & grupo!="PCOS"])),var.equal=F)$p.value)
p.valor.sujetos.obesos.control.adj <- p.adjust(p.valor.sujetos.obesos.control,"BH")

res_sujetos.obesos.control <- data.frame(p.valor.sujetos.obesos.control=p.valor.sujetos.obesos.control,p.valor.sujetos.obesos.control.adj=p.valor.sujetos.obesos.control.adj)


res_sujetos.obesos.control<- add_classificacion_metaboloma(res_sujetos.obesos.control,metaboloma.total)
res_sujetos.obesos.control$medias_Females <- apply((metaboloma.total[grupo[obesidad=="No Obese" & grupo!="PCOS"]=="Female",]),2,mean)
res_sujetos.obesos.control$medias_Males <- apply((metaboloma.total[grupo[obesidad=="No Obese" & grupo!="PCOS"]=="Male",]),2,mean)

res_sujetos.Noobesos.control <- res_sujetos.obesos.control[order(res_sujetos.obesos.control$p.valor.sujetos.obesos.control.adj),]

### VEMOS obesidad en las mujeres en genral

p.valor.mujeres.general <- apply(metaboloma.total[obesidad=="No Obese" & grupo!="Male",],2,function(x) t.test(x~as.factor(as.character(grupo[obesidad=="No Obese" & grupo!="Male"])),var.equal=F)$p.value)
p.valor.mujeres.general.adj <- p.adjust(p.valor.mujeres.general,"BH")

res_mujeres.general <- data.frame(p.valor.mujeres.general=p.valor.mujeres.general,p.valor.mujeres.general.adj=p.valor.mujeres.general.adj)


res_mujeres.general<- add_classificacion_metaboloma(res_mujeres.general,metaboloma.total)
res_mujeres.general$medias_PCOS <- apply((metaboloma.total[grupo[obesidad=="No Obese" & grupo!="Male"]=="PCOS",]),2,mean)
res_mujeres.general$medias_Females <- apply((metaboloma.total[grupo[obesidad=="No Obese" & grupo!="Male"]=="Female",]),2,mean)

res_mujeresNO.obesas.general <- res_mujeres.general[order(res_mujeres.general$p.valor.mujeres.general.adj),]


### VEMOS obesidad en las hombres vs PCOS en genral

p.valor.hombres.pcos <- apply(metaboloma.total[obesidad=="No Obese" & grupo!="Female",],2,function(x) t.test(x~as.factor(as.character(grupo[obesidad=="No Obese" & grupo!="Female"])),var.equal=F)$p.value)
p.valor.hombres.pcos.adj <- p.adjust(p.valor.hombres.pcos,"BH")

res_hombres.pcos <- data.frame(p.valor.hombres.pcos=p.valor.hombres.pcos,p.valor.hombres.pcos.adj=p.valor.hombres.pcos.adj)


res_hombres.pcos<- add_classificacion_metaboloma(res_hombres.pcos,metaboloma.total)
res_hombres.pcos$medias_PCOS <- apply((metaboloma.total[grupo[obesidad=="No Obese" & grupo!="Female"]=="PCOS",]),2,mean)
res_hombres.pcos$medias_males <- apply((metaboloma.total[grupo[obesidad=="No Obese" & grupo!="Female"]=="Male",]),2,mean)

res_hombres.pcos.noObesos <- res_hombres.pcos[order(res_hombres.pcos$p.valor.hombres.pcos.adj),]

fun_anova <- function(x){
  mod <- lm(x ~ grupo*obesidad,contrasts=list(grupo=contr.sum,obesidad=contr.sum,white.adjust=T))
  s <- Anova(mod,type = 3)$`Pr(>F)`
  ret <- c(grupo=s[2],obesidad=s[3],interaccion=s[4])
  return(ret)
}


metaboloma.total <- as.data.frame(t(datos$totales$metaboloma))
p.valor.obesidad <- apply(metaboloma.total,2,function(x) t.test(x~obesidad,var.equal=T)$p.value)
p.valor.obesidad.adj <- p.adjust(p.valor.obesidad,"BH")

res_obesidad <- data.frame(p.valor.obesidad=p.valor.obesidad,p.valor.obesidad.adj=p.valor.obesidad.adj)


res_obesidad<- add_classificacion_metaboloma(res_obesidad,metaboloma.total)
res_obesidad$medias_obesos <- apply((metaboloma.total[obesidad=="Obese",]),2,mean)
res_obesidad$medias_NOobesos <- apply((metaboloma.total[obesidad=="No Obese",]),2,mean)

res_obesidad <- res_obesidad[order(res_obesidad$p.valor.obesidad.adj),]



glmgeneral <- as.data.frame(t(apply((metaboloma.total),2, fun_anova)))

glmgeneral.grupo <- rownames(glmgeneral[glmgeneral$grupo<0.05,])
glmgeneral.obesidad <- rownames(glmgeneral[glmgeneral$obesidad<0.05,])
glmgeneral.interaccion <- rownames(glmgeneral[glmgeneral$interaccion<0.05,])


res_obesidad$glm_res <- glmgeneral$obesidad

res_obesidad$interpretacion <- ifelse(res_obesidad$medias_obesos>res_obesidad$medias_NOobesos,"AUMENTA","DISMINUYE")

res_obesidad$interpretacion.p <- ifelse(res_obesidad$p.valor.obesidad<0.05,"*","")
res_obesidad$interpretacion.padj <- ifelse(res_obesidad$p.valor.obesidad.adj<0.05,"**","")
res_obesidad$interpretacion.p.glm <- ifelse(res_obesidad$glm_res<0.05,"***","")
res_obesidad$AA <- rownames(res_obesidad)



res_sexo$interpretacion <- ifelse(res_sexo$Males>res_sexo$Females,"AUMENTA","DISMINUYE")
res_sexo$metabolitos <- rownames(res_sexo)
res_sexo$signif <- ifelse(res_sexo$p.valor.sexo<0.05,"*","")
res_sexo$signif.adj <- ifelse(res_sexo$p.valor.sexo.adj<0.05,"**","")


res_mujeres$interpretacion <- ifelse(res_mujeres.general$medias_PCOS>res_mujeres.general$medias_Females,"AUMENTA","DISMINUYE")
res_mujeres$metabolitos <- rownames(res_mujeres.general)
res_mujeres$signif <- ifelse(res_mujeres.general$p.valor.mujeres.general<0.05,"*","")
res_mujeres$signif.adj <- ifelse(res_mujeres.general$p.valor.mujeres.general.adj<0.05,"**","")


mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
p.valor.grupo.mujeres <- apply(metaboloma.total[grupo!="Male",],2,function(x) t.test(x~mujeres,var.equal=F)$p.value)
p.valor.grupo.mujeres.adj <- p.adjust(p.valor.grupo.mujeres,"BH")

res_grupo.mujeres <- data.frame(p.valor.grupo.mujeres=p.valor.grupo.mujeres,p.valor.grupo.mujeres.adj=p.valor.grupo.mujeres.adj)


res_grupo.mujeres<- add_classificacion_metaboloma(res_grupo.mujeres,metaboloma.total)
res_grupo.mujeres$PCOS <- apply((metaboloma.total[mujeres=="PCOS",]),2,mean)
res_grupo.mujeres$Females<- apply((metaboloma.total[mujeres=="Female",]),2,mean)

res_grupo.mujeres <- res_grupo.mujeres[order(res_grupo.mujeres$p.valor.grupo.mujeres.adj),]


res_grupo.mujeres$interpretacion <- ifelse(res_grupo.mujeres$PCOS>res_grupo.mujeres$Females,"AUMENTA","DISMINUYE")
res_grupo.mujeres$metabolitos <- rownames(res_grupo.mujeres)
res_grupo.mujeres$signif <- ifelse(res_grupo.mujeres$p.valor.grupo.mujeres<0.05,"*","")
res_grupo.mujeres$signif.adj <- ifelse(res_grupo.mujeres$p.valor.grupo.mujeres.adj<0.05,"**","")

