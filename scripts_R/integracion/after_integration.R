library(ALDEx2)
library(MOFA2)
library(psych)
library(ggheatmap)


source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
transformacion <- function(X) {
  return((t(log2(X / colSums(
    X
  )))))
}
datos  <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")


metaboloma <- (transformacion(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX

## Primero extraemos el mejor numero de componentes
modelos <- lapply(list.files("./scripts_R/integracion/quasi/",
                             full.names = T), load_model)

modelos <- lapply(list.files(path = "./scripts_R/integracion/modelo_optimo_comp",
                             full.names = T), load_model)

w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

### hasta aqui hemos obtenido el mejor modelo que nos dira el numero de componetes

### ahora cargamos los datos permutados
modelos <- lapply(list.files(path = "./scripts_R/integracion/quasi/",
                             full.names = T), load_model)
### obtenemos la distribucion nula
pesos.null <- lapply(modelos, function(x) get_weights(x,scale = T))
pesos.null <- lapply(pesos.null, function(x) Reduce(rbind,x))
### ahora de los del modelooriginal obtenelos los pesos
pesos.orig <- Reduce(rbind,get_weights(modelo,scale = T))
### ahora tomamos cuales de los valores permutacionales son mayores que la variable a testear

colMedian <- function(x) apply(x,1,median)
medianas <- lapply(pesos.null, colMedian)
medianas.red <- Reduce(cbind,medianas)
diferencias_mean <- Reduce("+",diferencias)/1000

diferencias.tmp.stats <- lapply(pesos.null, function (x) (a))
diferencias.tmp.stats <- Reduce("+",diferencias.tmp.stats)/1000

diferencias.tmp <- lapply(pesos.null, function (x) (x>pesos.orig))


diferencias_mean.tmp <- Reduce("+",diferencias.tmp)/1000

diferencias_sig_log <- (diferencias_mean.tmp<0.05)
rownames(diferencias_sig_log) <- rownames(diferencias_mean.tmp)
nombres <- names(which(apply(diferencias_sig_log,1,function(x) any(x==T))))


medianas_fac <- analyze_data(as.matrix(factores),grupo,obesidad)$medianas
poi <- pesos.orig[nombres,]


## tendrima una multiplicacion de las varsofi con factores con los pesos mas importatnes

reconstruccion <- factores %*% t(poi)
csa <- t(as.matrix(poi) %*% as.matrix(medianas_fac))
aux <- corr.test(rbind(reconstruccion,csa),method = "spearman")$r

csa2 <- (reconstruccion) %*% t(csa)










 