


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
datos  <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")


metaboloma <- (transformacion(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))

mofa.obj <- create_mofa_from_matrix(list(
  metaboloma = t(metaboloma),
  metagenoma = t(metagenoma)
))

data_opts <- get_default_data_options(mofa.obj)
data_opts$scale_views <- T

model_opts <- get_default_model_options(mofa.obj)

model_opts$num_factors <- 8

train_opts <- get_default_training_options(mofa.obj)

### probar estocastico

MOFAobject <- prepare_mofa(
  object = mofa.obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

directorio_modelo <- "./scripts_R/integracion/modelo_prueba"
if (!dir.exists(directorio_modelo)) {
  dir.create(directorio_modelo)
}



## Primero vemos cual es el numero optimo de componentes

# semillas <- seq(40,80,10)
# 
# directorio_modelo <- "./scripts_R/integracion/modelo_optimo_comp"
# if (!dir.exists(directorio_modelo)) {
#   dir.create(directorio_modelo)
# }
# mofa_componentes <- function(ncomp,semilla){
#   metaboloma <- (transformacion(datos$comunes$metaboloma))
#   metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
#   mofa.obj <- create_mofa_from_matrix(list(
#     metaboloma = t(metaboloma),
#     metagenoma = t(metagenoma)
#   ))
#   
#   data_opts <- get_default_data_options(mofa.obj)
#   data_opts$scale_views <- T
#   
#   model_opts <- get_default_model_options(mofa.obj)
#   
#   model_opts$num_factors <- ncomp
#   
#   train_opts <- get_default_training_options(mofa.obj)
#   train_opts$seed <- semilla
#   
#   ### probar estocastico
#   
#   MOFAobject <- prepare_mofa(
#     object = mofa.obj,
#     data_options = data_opts,
#     model_options = model_opts,
#     training_options = train_opts
#   )
#   
# 
#   outfile <- file.path(directorio_modelo, paste0(semilla,"_",ncomp,"model.hdf5"))
#   
#   MOFAobject.trained <- run_mofa(MOFAobject, outfile)
#   
# }
# ncomps <- 2:11
# for(s in semillas){
#   
#   lapply(ncomps, function(x) mofa_componentes(ncomp = x,semilla = s))
#   
# }

modelos <- lapply(list.files(path = "./scripts_R/integracion/modelo_optimo_comp",
                             full.names = T), load_model)

w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

factores <- length(colnames(get_factors(modelo)[[1]]))


pesos <- Reduce("rbind",get_weights(modelo,scale=T))

estadistico <- apply(abs(pesos), 1, sum)

etiquetas <- sample(1:ncol(pesos))

distribucion_nula <- replicate(1000, {
  etiquetas_aleatorias <- sample(1:ncol(pesos))
  apply(abs(pesos[, etiquetas_aleatorias]), 1, sum)
})

valores_p <- apply(distribucion_nula, 1, function(x) sum(x >= estadistico) / length(x))

variables_significativas <- which(valores_p < 0.05)

rownames(pesos)[variables_significativas]

factores <- get_factors(modelo)[[1]]
reconstruccion <- (factores) %*% t(pesos[variables_significativas,])

grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo<- datos$comunes$variables_in_bacteria$SEX
medias_latentes <- analyze_data(factores,grupo,obesidad)$medianas

prueba.cor <- corr.test((medias_latentes),t(pesos),method ="spearman",adjust = "none")$p

# prueba.cor <- apply(prueba.cor,2,function(x) p.adjust(x,"BH"))

feat_Extract <- which(apply(prueba.cor,2,function(x) any(x<0.05)))
# w <- which(prueba.cor<0.05)
nombres1 <- rownames(pesos)[variables_significativas]
nombres2 <- rownames(pesos)[feat_Extract]

nombres <- c(nombres1,nombres2)

reconstruccion <- factores %*% t(pesos[colnames(reconstruccion.scale)[coef(model)@i],])

pca <- prcomp(reconstruccion,scale. = F)

pca.plot <- data.frame(PC1=pca$x[,1],
                       PC2=pca$x[,2],
                       obesidad=obesidad,
                       grupo = grupo)

ggplot(pca.plot,aes(PC1,PC2,color=grupo))+geom_point()

library(glmnet)
reconstruccion.scale <- as.data.frame(cbind(metaboloma,metagenoma))
reconstruccion.scale$clase <- grupo
levels(reconstruccion.scale$clase) <- make.names(levels(reconstruccion.scale$clase))
idx <- sample(1:nrow(reconstruccion),nrow(reconstruccion)*0.6)
train.data <- reconstruccion.scale[idx,]
test.data <- reconstruccion.scale[-idx,]
library(caret)
x <- as.matrix(train.data[,-ncol(train.data)])
y <- train.data[,ncol(train.data)]
cv.lasso <- cv.glmnet(x,y ,alpha = 1, family = "multinomial",nfolds=10,grouped=F)
# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, family = "multinomial",
                lambda = cv.lasso$lambda.min)

# Display regression coefficients
coef(model)
# Make predictions on the test data
x.test <- model.matrix(test.data$clase ~., test.data)[,-1]
predicted.classes <-stats::predict(model,x.test,type = "class")
# Model accuracy
confusionMatrix(grupo[-idx],as.factor(predicted.classes),)

pcos.col <- colnames(reconstruccion.scale)[coef(model)$PCOS@i]
fem.col <-colnames(reconstruccion.scale)[coef(model)$Female@i]
male.col <- colnames(reconstruccion.scale)[coef(model)$Male@i]

nomnres2 <- c(pcos.col,fem.col,male.col)

reconstruccion <-  factores %*% t(pesos[male.col,])
pca <- prcomp(reconstruccion,scale. = T)

pca.plot <- data.frame(PC1=pca$x[,1],
                       PC2=pca$x[,2],
                       obesidad=obesidad,
                       grupo = grupo,
                       sexo=sexo)

ggplot(pca.plot,aes(PC1,PC2,color=sexo))+geom_point()

