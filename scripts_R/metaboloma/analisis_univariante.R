library(ggplot2)

library(reshape2)
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/calculo_medianas.R")
datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")
set.seed(126581)


X <- (datos$totales$metaboloma)
X <- t(log2(X/colSums(X)))

grupo <- datos$totales$general_data$GROUP
obesidad <- datos$totales$general_data$OBESE

resultado <- analyze_data(as.matrix(X),grupo,obesidad,correccion = 3)

resultado$interpretacion<-add_classificacion_metaboloma(resultado$interpretacion)
resultado$interpretacion <- resultado$interpretacion[order(resultado$interpretacion$grupo_meta),]

directorio.totales <- "./scripts_R/metaboloma/resultado.total"
if(!dir.exists(directorio.totales)){
  dir.create(directorio.totales)
}

write.csv(resultado$p.values,file.path(directorio.totales,"p_values.csv"))
write.csv(resultado$interpretacion  ,file.path(directorio.totales,"interpretacion.csv"))
write.csv(compute_medians(X,grupo,obesidad),file.path(directorio.totales,"medianas.csv"))

X <- (datos$comunes$metaboloma)
X <- t(log2(X/colSums(X)))


grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE

resultado <- analyze_data(as.matrix(X),grupo,obesidad,correccion = 3)

resultado$interpretacion<-add_classificacion_metaboloma(resultado$interpretacion)


resultado$interpretacion <- resultado$interpretacion[order(resultado$interpretacion$grupo_meta),]

directorio.comunes <- "./scripts_R/metaboloma/resultado.comun"
if(!dir.exists(directorio.comunes)){
  dir.create(directorio.comunes)
}

write.csv(resultado$p.values,file.path(directorio.comunes,"p_values.csv"))
write.csv(resultado$interpretacion  ,file.path(directorio.comunes,"interpretacion.csv"))
write.csv(compute_medians(X,grupo,obesidad),file.path(directorio.comunes,"medianas.csv"))
