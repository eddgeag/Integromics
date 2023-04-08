
library(ALDEx2)
library(limma)
library(Hotelling)
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/manova_vanvalen.R")
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/otras_funciones_utiles.R")

datos.obesidad <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_obesidad_.rds")

datos.grupo <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_grupo_.rds")

datos.sexo <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_sex_.rds")

datos.interaccion<- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_interaccion_.rds")

datos.none<- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_none_.rds")

datos <- list(none=datos.none,sexo=datos.sexo,grupo=datos.grupo,obesidad=datos.obesidad,interaccion=datos.interaccion)

fun_summary_genero <- function(X){
  obesidad <- X$comunes$variables_in_bacteria$OBESE
  grupo <- X$comunes$variables_in_bacteria$GROUP
  sexo <- X$comunes$variables_in_bacteria$SEX
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  escalados <- scale(t(mean_aldex(X$comunes$microbiota$aldex$genero.aldex)))
  
  pvalores <- vector("numeric",length=4)
  pvalores[1] <-adonis2(escalados~obesidad,method = "euclidean")$`Pr(>F)`[1]
  pvalores[2] <-adonis2(escalados~grupo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[3] <-adonis2(escalados~sexo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[4] <-adonis2(escalados~grupo*obesidad,method = "euclidean")$`Pr(>F)`[3]
  pvalores[5] <-adonis2(escalados[grupo!="Male",]~mujeres,method = "euclidean")$`Pr(>F)`[1]
  names(pvalores) <- c("obesidad","grupo","sexo","interaccion","mujeres")
  return(pvalores)
}

fun_summary_filo <- function(X){
  obesidad <- X$comunes$variables_in_bacteria$OBESE
  grupo <- X$comunes$variables_in_bacteria$GROUP
  sexo <- X$comunes$variables_in_bacteria$SEX
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  escalados <- scale(t(mean_aldex(X$comunes$microbiota$aldex$filo.aldex)))
  
  pvalores <- vector("numeric",length=4)
  pvalores[1] <-adonis2(escalados~obesidad,method = "euclidean")$`Pr(>F)`[1]
  pvalores[2] <-adonis2(escalados~grupo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[3] <-adonis2(escalados~sexo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[4] <-adonis2(escalados~grupo*obesidad,method = "euclidean")$`Pr(>F)`[3]
  pvalores[5] <-adonis2(escalados[grupo!="Male",]~mujeres,method = "euclidean")$`Pr(>F)`[1]
  names(pvalores) <- c("obesidad","grupo","sexo","interaccion","mujeres")
  return(pvalores)
}

fun_summary_filo.filt <- function(X){
  obesidad <- X$comunes$variables_in_bacteria$OBESE
  grupo <- X$comunes$variables_in_bacteria$GROUP
  sexo <- X$comunes$variables_in_bacteria$SEX
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  escalados <- scale(t(mean_aldex(X$comunes$microbiota$aldex.filt$filo.aldex)))
  
  pvalores <- vector("numeric",length=4)
  pvalores[1] <-adonis2(escalados~obesidad,method = "euclidean")$`Pr(>F)`[1]
  pvalores[2] <-adonis2(escalados~grupo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[3] <-adonis2(escalados~sexo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[4] <-adonis2(escalados~grupo*obesidad,method = "euclidean")$`Pr(>F)`[3]
  pvalores[5] <-adonis2(escalados[grupo!="Male",]~mujeres,method = "euclidean")$`Pr(>F)`[1]
  names(pvalores) <- c("obesidad","grupo","sexo","interaccion","mujeres")
  return(pvalores)
}

fun_summary_genero.filt <- function(X){
  obesidad <- X$comunes$variables_in_bacteria$OBESE
  grupo <- X$comunes$variables_in_bacteria$GROUP
  sexo <- X$comunes$variables_in_bacteria$SEX
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  escalados <- scale(t(mean_aldex(X$comunes$microbiota$aldex.filt$genero.aldex)))
  
  pvalores <- vector("numeric",length=4)
  pvalores[1] <-adonis2(escalados~obesidad,method = "euclidean")$`Pr(>F)`[1]
  pvalores[2] <-adonis2(escalados~grupo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[3] <-adonis2(escalados~sexo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[4] <-adonis2(escalados~grupo*obesidad,method = "euclidean")$`Pr(>F)`[3]
  pvalores[5] <-adonis2(escalados[grupo!="Male",]~mujeres,method = "euclidean")$`Pr(>F)`[1]
  names(pvalores) <- c("obesidad","grupo","sexo","interaccion","mujeres")
  return(pvalores)
}


fun_summary_genero.bayes <- function(X){
  obesidad <- X$comunes$variables_in_bacteria$OBESE
  grupo <- X$comunes$variables_in_bacteria$GROUP
  sexo <- X$comunes$variables_in_bacteria$SEX
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  escalados <- scale(((X$comunes$microbiota$bayes.clr$genero.clr)))
  
  pvalores <- vector("numeric",length=4)
  pvalores[1] <-adonis2(escalados~obesidad,method = "euclidean")$`Pr(>F)`[1]
  pvalores[2] <-adonis2(escalados~grupo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[3] <-adonis2(escalados~sexo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[4] <-adonis2(escalados~grupo*obesidad,method = "euclidean")$`Pr(>F)`[3]
  pvalores[5] <-adonis2(escalados[grupo!="Male",]~mujeres,method = "euclidean")$`Pr(>F)`[1]
  names(pvalores) <- c("obesidad","grupo","sexo","interaccion","mujeres")
  return(pvalores)
}

fun_summary_filo.bayes <- function(X){
  obesidad <- X$comunes$variables_in_bacteria$OBESE
  grupo <- X$comunes$variables_in_bacteria$GROUP
  sexo <- X$comunes$variables_in_bacteria$SEX
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  escalados <- scale(((X$comunes$microbiota$bayes.clr$filo.clr)))
  
  pvalores <- vector("numeric",length=4)
  pvalores[1] <-adonis2(escalados~obesidad,method = "euclidean")$`Pr(>F)`[1]
  pvalores[2] <-adonis2(escalados~grupo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[3] <-adonis2(escalados~sexo,method = "euclidean")$`Pr(>F)`[1]
  pvalores[4] <-adonis2(escalados~grupo*obesidad,method = "euclidean")$`Pr(>F)`[3]
  pvalores[5] <-adonis2(escalados[grupo!="Male",]~mujeres,method = "euclidean")$`Pr(>F)`[1]
  names(pvalores) <- c("obesidad","grupo","sexo","interaccion","mujeres")
  return(pvalores)
}

resultados.genero <- as.data.frame(lapply(datos, fun_summary_genero))

resultados.filo <- as.data.frame(lapply(datos, fun_summary_filo))

resultados.filo.filt <- as.data.frame(lapply(datos, fun_summary_filo.filt))

resultados.genero.filt <- as.data.frame(lapply(datos, fun_summary_genero.filt))

resultados.genero.bayes <- as.data.frame(lapply(datos, fun_summary_genero.bayes))

resultados.filo.bayes <- as.data.frame(lapply(datos, fun_summary_filo.bayes))


none <- scale(t(mean_aldex(datos$none$comunes$microbiota$aldex$genero.aldex)))
none.catego <- as.data.frame(none)
sexo <- datos$none$comunes$variables_in_bacteria$SEX
none.catego$sexo <- sexo

none.melt <- melt(none.catego)

ggplot(none.melt,aes(value,color=sexo))+geom_density()

pcx <- prcomp(t(none))
varianza <-100* pcx$sdev^2/sum(pcx$sdev^2)

dataplot <- bind_cols(pcx$rotation[,1:2],obesidad = sexo)

ggplot(dataplot,aes(PC1,PC2,color=obesidad))+geom_point()


mds <- cmdscale(dist(none))
dataplot2 <- bind_cols(PC1=mds[,1],PC2=mds[,2],obesidad = sexo)
ggplot(dataplot2,aes(PC1,PC2,color=obesidad))+geom_point()







