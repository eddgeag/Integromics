library(ggplot2)
library(vegan)
library(reshape2)
library(ALDEx2)
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/calculo_medianas.R")
datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")

set.seed(126581)
sexo <- datos$comunes$variables_in_bacteria$SEX
grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE
genero <- t(mean_aldex(datos$comunes$microbiota$genero))
filo <- t(mean_aldex(datos$comunes$microbiota$filo))
resultado.genero <- analyze_data(genero,grupo,obesidad,correccion=1)
medianas.genero <- compute_medians(genero,grupo,obesidad)
resultado.filo <- analyze_data(filo,grupo,obesidad,correccion=1)
medianas.filo<- compute_medians(filo,grupo,obesidad)

directorio.filo <- "./scripts_R/metagenoma/resultados_univariante_filo"
if(!dir.exists(directorio.filo)){
  dir.create(directorio.filo)
}
directorio.genero <- "./scripts_R/metagenoma/resultados_univariante_genero"
if(!dir.exists(directorio.genero)){
  dir.create(directorio.genero)
}

write.csv(resultado.filo$p_valores,file.path(directorio.filo,"p_valores.csv"))
write.csv(resultado.filo$interpretacion,file.path(directorio.filo,"interpretacion.csv"))
write.csv(medianas.filo,file.path(directorio.filo,"medianas.csv"))



write.csv(resultado.genero$p_valores,file.path(directorio.genero,"p_valores.csv"))
write.csv(resultado.genero$interpretacion,file.path(directorio.genero,"interpretacion.csv"))
write.csv(medianas.genero,file.path(directorio.genero,"medianas.csv"))

