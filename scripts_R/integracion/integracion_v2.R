library(glmnet)
library(caret)
library(ALDEx2)
library(MOFA2)
library(psych)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")

transformacion <- function(X) {
  return((t(log2(X / colSums(
    X
  )))))
}

funcion_pca <- function(X){
  pca <- prcomp(X,scale. = T)
  
  pca.plot <- data.frame(PC1=pca$x[,1],
                         PC2=pca$x[,2],
                         obesidad=obesidad,
                         grupo = grupo,
                         sexo=sexo)
  p1 <- ggplot(pca.plot,aes(PC1,PC2,color=grupo))+geom_point()
  p2 <- ggplot(pca.plot,aes(PC1,PC2,color=grupo))+geom_point()+facet_grid(~obesidad)
  p3 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+geom_point()+facet_grid(~obesidad)
  p4 <- ggarrange(p1,p2,p3)
  p1 <- ggplot(pca.plot,aes(PC1,PC2,color=obesidad))+geom_point()
  p2 <- ggplot(pca.plot,aes(PC1,PC2,color=obesidad))+geom_point()+facet_grid(~grupo)
  p3 <- ggplot(pca.plot,aes(PC1,PC2,color=obesidad))+geom_point()+facet_grid(~sexo)
  
  
  p5 <- ggarrange(p1,p2,p3)
  p1 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+geom_point()
  p2 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+geom_point()+facet_grid(~grupo)
  p3 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+geom_point()+facet_grid(~sexo)
  
  
  p6<- ggarrange(p1,p2,p3)
  
  return(list(p4,p5,p6))
}
datos <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")
obesidad <- datos$comunes$obesidad
grupo <- datos$comunes$grupo
sexo <- datos$comunes$variables_in_bacteria$SEX

metaboloma <- (transformacion(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))

### OBTENEMOS EL MODELO MOFA Y PESOS SIGNIFICATIVOS

modelos <- lapply(list.files(path = "./scripts_R/integracion/modelo_optimo_comp",
                             full.names = T), load_model)

w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

factores <- length(colnames(get_factors(modelo)[[1]]))


pesos <- Reduce("rbind",get_weights(modelo,scale=T))

estadistico <- apply(abs(pesos), 1, sum)

etiquetas <- sample(1:ncol(pesos))

distribucion_nula <- replicate(10000, {
  etiquetas_aleatorias <- sample(1:ncol(pesos))
  apply(abs(pesos[, etiquetas_aleatorias]), 1, sum)
})

valores_p <- apply(distribucion_nula, 1, function(x) sum(x >= estadistico) / length(x))

variables_significativas <- which(valores_p < 0.05)

rownames(pesos)[variables_significativas]

factores <- get_factors(modelo)[[1]]
reconstruccion <- (factores) %*% t(pesos[variables_significativas,])

modelo.splsda <- block.splsda(X=list(metaboloma=factores %*% t(pesos),
                                     metagenoma=factores %*% t(pesos)),
                                     as.numeric(grupo),ncomp = 5,scale = T)

modelo.perf.splsda <- perf(modelo.splsda,auc=T,folds=2,
                           signif.threshold = 0.05,nrepeat = 10)

ncomps <- modelo.perf.splsda$choice.ncomp

modelo.tune <- tune.splsda(X = (factores %*% t(pesos)),
                           Y=grupo,ncomp = ncomps,
                           test.keepX  = which(valores_p<0.05),
                           validation = "Mfold",folds=10,progressBar = T)
modelo.tune <- splsda(X = (factores %*% t(pesos)),
                           Y=grupo,ncomp = ncomps,
                            keepX = modelo.tune$choice.keepX)
csa <- unique(rownames(pesos)[c(which(modelo.tune$loadings$X[,1]>0),
                         which(modelo.tune$loadings$X[,2]>0))])

datos.plot <- data.frame(PC1=modelo.tune$variates$X[,1],
                         PC2=modelo.tune$variates$X[,2],
                         grupo=grupo,
                         obesidad=obesidad,
                         sexo=sexo)
ggplot(datos.plot,aes(PC1,PC2,color=grupo))+geom_point()

diablos <- readRDS("./scripts_R/integracion/diablos.rds")


diablo.grupo <- diablos$grupo
diablo.obesidad <- diablos$obesidad

datos.plot <- data.frame(PC1=diablo.grupo$variates$metaboloma[,1],
                         PC2=diablo.grupo$variates$metagenoma[,2],
                         grupo=grupo,
                         obesidad=obesidad,
                         sexo=sexo)
ggplot(datos.plot,aes(PC1,PC2,color=obesidad))+geom_point(size=2)
pesos.diablo <- unique(c(rownames(diablo.grupo$loadings$metaboloma)[which(apply(diablo.grupo$loadings$metaboloma,
                                                                                1,
                                                                                function(x) any(x>0)))
                                                                    ],
                             rownames(diablo.grupo$loadings$metagenoma)[which(apply(diablo.grupo$loadings$metagenoma,
                                                                                    1,
                                                                                    function(x) any(x>0)))
                             ],
                             rownames(diablo.obesidad$loadings$metaboloma)[which(apply(diablo.obesidad$loadings$metaboloma,
                                                                                       1,
                                                                                       function(x) any(x>0)))
                             ],
                             rownames(diablo.obesidad$loadings$metagenoma)[which(apply(diablo.grupo$loadings$metagenoma,
                                                                                1,
                                                                                function(x) any(x>0)))
                                                                    ]))
pesos.metaboloma <- pesos.diablo[grep("mean",x=pesos.diablo)]
pesos.metagenoma <- pesos.diablo[-grep("mean",x=pesos.diablo)]

diablo.entero <- block.splsda(X=list(metaboloma=metaboloma,
                                   metagenoma=metagenoma
                                   ),
                              grupo,
                              ncomp=2,
                              keepX = list(metaboloma=c(22,22),
                                           metagenoma=c(106,106)))



datos.plot <- data.frame(PC1=diablo.entero$variates$metagenoma[,2],
                         PC2=diablo.entero$variates$metaboloma[,1],
                         grupo=grupo,
                         obesidad=obesidad,
                         sexo=sexo)
ggplot(datos.plot,aes(PC1,PC2,color=grupo,shape=obesidad))+geom_point(size=2)






