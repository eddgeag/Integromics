
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")
set.seed(126581)

directorio <- "./scripts_R/IP/resultados_univariantes_totales"

if(!dir.exists(directorio)){
  dir.create(directorio)
}

X <- scale(as.matrix(datos$totales$ip))
grupo <- datos$totales$general_data$GROUP
obesidad <- datos$totales$general_data$OBESE

resultado <- analyze_data(as.matrix(X),grupo,obesidad,correccion = 3)
write.csv(resultado$p_valores,file.path(directorio,"p_valores.csv"))
write.csv(resultado$interpretacion,file.path(directorio,"interpretacion.csv"))
write.csv(compute_medians(X,grupo,obesidad),file.path(directorio,"medianas.csv"))


directorio <- "./scripts_R/IP/resultados_univariantes_comunes"

if(!dir.exists(directorio)){
  dir.create(directorio)
}
X <- scale(as.matrix(datos$comunes$ip))
grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE

resultado <- analyze_data(as.matrix(X),grupo,obesidad,correccion = 3)

write.csv(resultado$p_valores,file.path(directorio,"p_valores.csv"))
write.csv(resultado$interpretacion,file.path(directorio,"interpretacion.csv"))
write.csv(compute_medians(X,grupo,obesidad),file.path(directorio,"medianas.csv"))



