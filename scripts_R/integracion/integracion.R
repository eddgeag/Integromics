library(glmnet)
library(caret)
library(ALDEx2)
library(MOFA2)
library(psych)
library(vegan)
library(mixOmics)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
transformacion <- function(X) {
  return((t(log2(X / colSums(
    X
  )))))
}

mofa_componentes <- function(ncomp,semilla){
  metaboloma <- (transformacion(datos$comunes$metaboloma))
  metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
  mofa.obj <- create_mofa_from_matrix(list(
    metaboloma = t(metaboloma),
    metagenoma = t(metagenoma)
  ))
  directorio_modelo <- "./scripts_R/integracion/modelo_optimo_comp"
  if (!dir.exists(directorio_modelo)) {
    dir.create(directorio_modelo)
  }
  data_opts <- get_default_data_options(mofa.obj)
  data_opts$scale_views <- T
  
  model_opts <- get_default_model_options(mofa.obj)
  
  model_opts$num_factors <- ncomp
  
  train_opts <- get_default_training_options(mofa.obj)
  train_opts$seed <- semilla
  
  ### probar estocastico
  
  MOFAobject <- prepare_mofa(
    object = mofa.obj,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  
  outfile <- file.path(directorio_modelo, paste0(semilla,"_",ncomp,"model.hdf5"))
  
  MOFAobject.trained <- run_mofa(MOFAobject, outfile)
  
}

datos  <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")


metaboloma <- (transformacion(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX

## Primero vemos cual es el numero optimo de componentes

semillas <- seq(40,80,10)

directorio_modelo <- "./scripts_R/integracion/modelo_optimo_comp"
if (!dir.exists(directorio_modelo)) {
  dir.create(directorio_modelo)
}

ncomps <- 2:11
if(length(list.files(directorio_modelo))==0){
  for(s in semillas){
    
    lapply(ncomps, function(x) mofa_componentes(ncomp = x,semilla = s))
    
  }
}


modelos <- lapply(list.files(path = directorio_modelo,
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

#### Entrenamos al DIABLO
### obesidad
set.seed(123456789)
idx <- sample(1:nrow(metaboloma),size = nrow(metaboloma)*0.7)

train.met <- metaboloma[idx,]
train.mic <- metagenoma[idx,]
test.met <- metaboloma[-idx,]
test.mic <- metagenoma[-idx,]
grupo.train <- grupo[idx]
obesidad.train <- obesidad[idx]
grupo.test <- grupo[-idx]
obesidad.test <- obesidad[-idx]

X.list <- list(metaboloma=as.matrix(train.met),metagenoma=as.matrix(train.mic))
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))


model.preperf <- block.splsda(X=list(metaboloma=train.met,
                                     metaganoma=train.mic),
                              Y=obesidad.train,ncomp = 10,design = diseno)
model.perf <- perf(model.preperf,dist = "all",validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05)

ncomps <- 3
tune.diablo <- tune.block.splsda(list(metaboloma=train.met,
                                      metaganoma=train.mic),
                                 Y = obesidad.train,design = diseno,
                                 ncomp = ncomps,
                                 test.keepX = list(metaboloma=which(colnames(metaboloma) %in% names(variables_significativas)),
                                                   metagenoma=which(colnames(metagenoma) %in% names(variables_significativas))),
                                 validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05,progressBar = T)


keep.list <- tune.diablo$choice.keepX
modelo.diablo.obesidad <- block.splsda(list(metaboloma=train.met,
                                   metaganoma=train.mic),
                              Y=obesidad.train,
                              design = diseno,ncomp = ncomps,keepX = keep.list)


predict.diablo.obesidad <-  stats::predict(modelo.diablo.obesidad,list(metaboloma=test.met,
                                                     metaganoma=test.mic),type="class")
(confusion.mat.obesidad = caret::confusionMatrix(data = as.factor(predict.diablo.obesidad$class$max.dist$metaboloma[,2]),
                                       reference=obesidad.test))$overall[1]

### sexo
set.seed(123456789)
idx <- sample(1:nrow(metaboloma),size = nrow(metaboloma)*0.7)

train.met <- metaboloma[idx,]
train.mic <- metagenoma[idx,]
test.met <- metaboloma[-idx,]
test.mic <- metagenoma[-idx,]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]

X.list <- list(metaboloma=as.matrix(train.met),metagenoma=as.matrix(train.mic))
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))


model.preperf <- block.splsda(X=list(metaboloma=train.met,
                                     metaganoma=train.mic),
                              Y=sexo.train,ncomp = 1,design = diseno)
model.perf <- perf(model.preperf,dist = "all",validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05)

ncomps <- 3

tune.diablo <- tune.block.splsda(list(metaboloma=train.met,
                                      metaganoma=train.mic),
                                 Y = sexo.train,design = diseno,
                                 ncomp = ncomps,
                                 test.keepX = list(metaboloma=which(colnames(metaboloma) %in% names(variables_significativas)),
                                                   metagenoma=which(colnames(metagenoma) %in% names(variables_significativas))),
                                 validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05,progressBar = T)


keep.list <- tune.diablo$choice.keepX
modelo.diablo.sexo <- block.splsda(list(metaboloma=train.met,
                                         metaganoma=train.mic),
                                    Y=sexo.train,
                                    design = diseno,ncomp = ncomps,keepX = keep.list)


predict.diablo.sexo <-  stats::predict(modelo.diablo.sexo,list(metaboloma=test.met,
                                                     metaganoma=test.mic),type="class")
(confusion.mat.sexo = caret::confusionMatrix(data = as.factor(predict.diablo.sexo$class$max.dist$metaboloma[,2]),
                                                 reference=sexo.test))$ovarall[1]

### Tenemos que ver las contribuciones

csa <- plotLoadings(modelo.diablo.obesidad,plot = F,contrib = "max",method = "median",block = 2)

library(Mik)





