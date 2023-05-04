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

modelos <- lapply(list.files(path = "../integromics_code/scripts_R/integracion/modelo_optimo_comp/",
                             full.names = T), load_model)

null.modelos <- lapply(list.files(path="../integromics_code/scripts_R/integracion/quasi/",full.names = T),load_model)
w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

factores <- get_factors(modelo,scale = T)[[1]]
pesos <- Reduce(rbind,get_weights(modelo))
pesos.metaboloma <- get_weights(modelo)[[1]]

pesos.null <- lapply(null.modelos, function(x) Reduce(rbind,get_weights(x)))

mean.diff <- lapply(pesos.null, function(x) abs(x-pesos))


tensor <- (simplify2array(pesos.null))

diferencias <- lapply(pesos.null, function(x) abs(x-pesos)==0)

diferencias <- Reduce("+",diferencias)/length(pesos.null)

diferencias <- diferencias[complete.cases(diferencias),]

### Ahora obetenmos de estas diferencias, las significantes en factores

### partimos en entrenamiento y test

set.seed(123456789)
idx <- sample(1:nrow(metaboloma),nrow(metaboloma)*0.7,replace = F)
factores.train <- factores[idx,]
obesidad.train <- obesidad[idx]
factores.test <-factores[-idx,]
obesidad.test <-obesidad[-idx]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]
### regresion logistica binaria lasso para obesidad en factores
lambda <- 10^seq(10, -2, length = 100)
cv.lasso <- glmnet::cv.glmnet(factores.train,obesidad.train ,
                              alpha = 1, 
                              family = "binomial",
                              nfolds=10,grouped=F,lambda = lambda)
model.obesidad.factores <- glmnet(factores.train,
                                  obesidad.train, 
                                  alpha = 1, 
                                  family = "binomial",
                                  lambda = cv.lasso$lambda.min)


factores.clase.test <- bind_cols(factores.test,obesidad=obesidad.test)
x.test <- model.matrix(obesidad ~., factores.clase.test)[,-1]
predicted.classes <-as.factor(stats::predict(model.obesidad.factores,
                                             newx = x.test,type="class"))

(confusion <- confusionMatrix(predicted.classes,obesidad.test)$overall[1])
                                  
colnames(factores)[coef(model.obesidad.factores)@i]

### regresion logistica binaria lasso para sexo en factores
lambda <- 10^seq(10, -2, length = 100)
cv.lasso <- glmnet::cv.glmnet(factores.train,sexo.train ,
                              alpha = 1, 
                              family = "binomial",
                              nfolds=10,grouped=F,
                              lambda = lambda)
model.sexo.factores <- glmnet(factores.train,
                                  sexo.train, 
                                  alpha = 1, 
                                  family = "binomial",
                                  lambda = cv.lasso$lambda.min)


factores.clase.test <- bind_cols(factores.test,sexo=sexo.test)
x.test <- model.matrix(sexo ~., factores.clase.test)[,-1]
predicted.classes <-as.factor(stats::predict(model.sexo.factores,
                                             newx = x.test,type="class"))

(confusion <- confusionMatrix(predicted.classes,sexo.test)$overall[1])

colnames(factores)[coef(model.sexo.factores)@i]

factores.nombres <- unique(c(colnames(factores)[coef(model.sexo.factores)@i],
                             colnames(factores)[coef(model.obesidad.factores)@i]))
### Hacemos la reconstruccion

reconstruccion  <- factores[,factores.nombres] %*% t(diferencias[,factores.nombres])


### Realizamos seleccion de varaibles

### partimos en entrenamiento y test

set.seed(123456789)
idx <- sample(1:nrow(reconstruccion),nrow(reconstruccion)*0.7,replace = F)
reconstruccion.train <- reconstruccion[idx,]
obesidad.train <- obesidad[idx]
reconstruccion.test <-reconstruccion[-idx,]
obesidad.test <-obesidad[-idx]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]
### regresion logistica binaria lasso para obesidad en factores
lambda <- 10^seq(10, -2, length = 100)
cv.lasso <- glmnet::cv.glmnet(reconstruccion.train,obesidad.train ,
                              alpha = 1, 
                              family = "binomial",
                              nfolds=10,grouped=F,lambda = lambda)
model.obesidad.factores <- glmnet(reconstruccion.train,
                                  obesidad.train, 
                                  alpha = 1, 
                                  family = "binomial",
                                  lambda = cv.lasso$lambda.min)


factores.clase.test <- bind_cols(reconstruccion.test,obesidad=obesidad.test)
x.test <- model.matrix(obesidad ~., factores.clase.test)[,-1]
predicted.classes <-as.factor(stats::predict(model.obesidad.factores,
                                             newx = x.test,type="class"))

(confusion <- confusionMatrix(predicted.classes,obesidad.test)$overall[1])

colnames(reconstruccion)[coef(model.obesidad.factores)@i]

### regresion logistica binaria lasso para sexo en factores
lambda <- 10^seq(10, -2, length = 100)
cv.lasso <- glmnet::cv.glmnet(reconstruccion.train,sexo.train ,
                              alpha = 1, 
                              family = "binomial",
                              nfolds=10,grouped=F,
                              lambda = lambda)
model.sexo.factores <- glmnet(reconstruccion.train,
                              sexo.train, 
                              alpha = 1, 
                              family = "binomial",
                              lambda = cv.lasso$lambda.min)


factores.clase.test <- bind_cols(reconstruccion.test,sexo=sexo.test)
x.test <- model.matrix(sexo ~., factores.clase.test)[,-1]
predicted.classes <-as.factor(stats::predict(model.sexo.factores,
                                             newx = x.test,type="class"))

(confusion <- confusionMatrix(predicted.classes,sexo.test)$overall[1])

colnames(reconstruccion)[coef(model.sexo.factores)@i]

pesos.nombres <- unique(c(colnames(reconstruccion)[coef(model.sexo.factores)@i],
                             colnames(reconstruccion)[coef(model.obesidad.factores)@i]))
### Hacemos la reconstruccion

reconstruccion_2  <- factores[,factores.nombres] %*% t(diferencias[pesos.nombres,factores.nombres])

reconstruccion_total <- factores %*% t(pesos)
### Realizamos seleccion de varaibles
funcion_pca(reconstruccion)
funcion_pca(reconstruccion_total)
funcion_pca(reconstruccion_2)

### Ahora viene la obtencion de coeficientes sin mofa2


set.seed(123456789)
X <- scale(cbind(metaboloma,metagenoma))
idx <- sample(1:nrow(X),nrow(X)*0.7,replace = F)
X.train <- X[idx,]
obesidad.train <- obesidad[idx]
X.test <-X[-idx,]
obesidad.test <-obesidad[-idx]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]
### regresion logistica binaria lasso para obesidad en factores
lambda <- 10^seq(10, -2, length = 100)
cv.lasso <- glmnet::cv.glmnet(X.train,obesidad.train ,
                              alpha = 1, 
                              family = "binomial",
                              nfolds=10,grouped=F,lambda = lambda)
model.obesidad.factores <- glmnet(X.train,
                                  obesidad.train, 
                                  alpha = 1, 
                                  family = "binomial",
                                  lambda = cv.lasso$lambda.min)


factores.clase.test <- bind_cols(X.test,obesidad=obesidad.test)
x.test <- model.matrix(obesidad ~., factores.clase.test)[,-1]
predicted.classes <-as.factor(stats::predict(model.obesidad.factores,
                                             newx = x.test,type="class"))

(confusion <- confusionMatrix(predicted.classes,obesidad.test)$overall[1])

pesos_obesidad <- colnames(X)[coef(model.obesidad.factores)@i]
funcion_pca(X[,pesos_obesidad])
### regresion logistica binaria lasso para sexo en factores
lambda <- 10^seq(10, -2, length = 100)
cv.lasso <- glmnet::cv.glmnet(X.train,sexo.train ,
                              alpha = 1, 
                              family = "binomial",
                              nfolds=10,grouped=F,
                              lambda = lambda)
model.sexo.factores <- glmnet(X.train,
                              sexo.train, 
                              alpha = 1, 
                              family = "binomial",
                              lambda = cv.lasso$lambda.min)


factores.clase.test <- bind_cols(X.test,sexo=sexo.test)
x.test <- model.matrix(sexo ~., factores.clase.test)[,-1]
predicted.classes <-as.factor(stats::predict(model.sexo.factores,
                                             newx = x.test,type="class"))

(confusion <- confusionMatrix(predicted.classes,sexo.test)$overall[1])

pesos.sexo <- colnames(X)[coef(model.sexo.factores)@i]

pesos.nombres <- unique(c(colnames(X)[coef(model.sexo.factores)@i],
                          colnames(X)[coef(model.obesidad.factores)@i]))

funcion_pca(X[,pesos.nombres])

reconstruccion_5 <- factores[,factores.nombres] %*% t(pesos[pesos.nombres,factores.nombres])

funcion_pca(reconstruccion_5)

#### Ahora a hacer lo mismo pero con PLS-DA multibloque

library(mixOmics) # import the mixOmics library

X.list <- list(metaboloma=as.matrix(metaboloma),metagenoma=as.matrix(metagenoma))
Y <- obesidad
# names(Y) <- rownames(metaboloma)
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))
colnames(diseno) <- c("metaboloma","metagenoma")
rownames(diseno) <- colnames(diseno)
diablo.grupo <- mixOmics::block.plsda(X.list, as.numeric(obesidad),design = diseno,ncomp = 10) # run the method


perf.diablo <- perf(diablo.grupo, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values


ncomp <- perf.diablo$choice.ncomp$AveragedPredict[2,] # what is the optimal value of components according to perf()

pesos.totales <- unique(c(pesos_obesidad,pesos.sexo))
list.keepX <- list(metaboloma=which(colnames(metaboloma) %in% pesos.totales),
                   metagenoma=which(colnames(metagenoma) %in% pesos.totales))




tune.splsda.srbct <- tune.block.splsda(X.list, Y, ncomp = ncomp, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,progressBar = T,
                                 design = diseno) #'

list.keepX <- tune.splsda.srbct$choice.keepX

final.diablo.model = block.splsda(X = X.list, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = diseno)



datos.plot <- as.data.frame((bind_cols(final.diablo.model$variates$metaboloma[,1:2],final.diablo.model$variates$metagenoma[,1:2],
          grupo=grupo,obesidad=obesidad,sexo=sexo)))

colnames(datos.plot) <-c("mPC1","mPC2","gPC1","gPC2","grupo","obesidad","sexo")

ggplot(datos.plot,aes(mPC1,gPC1,color=obesidad))+geom_point()

diablo.obesidad.metaboloma <- final.diablo.model$loadings$metaboloma
diablo.obesidad.metagenoma <- final.diablo.model$loadings$metaboloma


##### Grupo


Y <- grupo
# names(Y) <- rownames(metaboloma)
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))
colnames(diseno) <- c("metaboloma","metagenoma")
rownames(diseno) <- colnames(diseno)
diablo.grupo <- mixOmics::block.plsda(X.list, as.factor(grupo),design = diseno,ncomp = 10) # run the method


perf.diablo <- perf(diablo.grupo, validation = "Mfold", 
                    folds = 5, nrepeat = 10, # use repeated cross-validation
                    progressBar = FALSE, auc = TRUE) # include AUC values


ncomp <- perf.diablo$choice.ncomp$AveragedPredict[2,] # what is the optimal value of components according to perf()

pesos.totales <- unique(c(pesos_obesidad,pesos.sexo))
list.keepX <- list(metaboloma=which(colnames(metaboloma) %in% pesos.totales),
                   metagenoma=which(colnames(metagenoma) %in% pesos.totales))




tune.splsda.srbct <- tune.block.splsda(X.list, Y, ncomp = ncomp, # calculate for first 4 components
                                       validation = 'Mfold',
                                       folds = 5, nrepeat = 10, # use repeated cross-validation
                                       dist = 'max.dist', # use max.dist measure
                                       measure = "BER", # use balanced error rate of dist measure
                                       test.keepX = list.keepX,progressBar = T,
                                       design = diseno) #

list.keepX <- tune.splsda.srbct$choice.keepX

final.diablo.model.grupo = block.splsda(X = X.list, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = diseno)



datos.plot <- as.data.frame((bind_cols(final.diablo.model.grupo$variates$metaboloma[,1:2],final.diablo.model.grupo$variates$metagenoma[,1:2],
                                       grupo=grupo,obesidad=obesidad,sexo=sexo)))

colnames(datos.plot) <-c("mPC1","mPC2","gPC1","gPC2","grupo","obesidad","sexo")

ggplot(datos.plot,aes(mPC1,gPC1,color=obesidad))+geom_point()

diablo.grupo.metaboloma <- final.diablo.model.grupo$loadings$metaboloma
diablo.grupo.metagenoma <- final.diablo.model.grupo$loadings$metaboloma

##### sexo


Y <- sexo
# names(Y) <- rownames(metaboloma)
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))
colnames(diseno) <- c("metaboloma","metagenoma")
rownames(diseno) <- colnames(diseno)
diablo.sexo <- mixOmics::block.plsda(X.list, as.factor(sexo),design = diseno,ncomp = 10) # run the method


perf.diablo <- perf(diablo.sexo, validation = "Mfold", 
                    folds = 5, nrepeat = 10, # use repeated cross-validation
                    progressBar = FALSE, auc = TRUE) # include AUC values


ncomp <- perf.diablo$choice.ncomp$AveragedPredict[2,] # what is the optimal value of components according to perf()

pesos.totales <- unique(c(pesos_obesidad,pesos.sexo))
list.keepX <- list(metaboloma=which(colnames(metaboloma) %in% pesos.totales),
                   metagenoma=which(colnames(metagenoma) %in% pesos.totales))




tune.splsda.srbct <- tune.block.splsda(X.list, Y, ncomp = ncomp, # calculate for first 4 components
                                       validation = 'Mfold',
                                       folds = 5, nrepeat = 10, # use repeated cross-validation
                                       dist = 'max.dist', # use max.dist measure
                                       measure = "BER", # use balanced error rate of dist measure
                                       test.keepX = list.keepX,progressBar = T,
                                       design = diseno) #

list.keepX <- tune.splsda.srbct$choice.keepX

final.diablo.model.sexo = block.splsda(X = X.list, Y = Y, ncomp = ncomp, 
                                        keepX = list.keepX, design = diseno)



datos.plot <- as.data.frame((bind_cols(final.diablo.model.grupo$variates$metaboloma[,1:2],final.diablo.model.grupo$variates$metagenoma[,1:2],
                                       grupo=grupo,obesidad=obesidad,sexo=sexo)))

colnames(datos.plot) <-c("mPC1","mPC2","gPC1","gPC2","grupo","obesidad","sexo")

ggplot(datos.plot,aes(mPC1,gPC2,color=grupo))+geom_point()+facet_grid(~obesidad)

diablo.sexo.metaboloma <- final.diablo.model.sexo$loadings$metaboloma
diablo.sexo.metagenoma <- final.diablo.model.sexo$loadings$metaboloma

diablos <- list(obesidad=final.diablo.model,
                grupo=final.diablo.model.grupo,
                sexo=final.diablo.model.sexo)

saveRDS(diablos,"./scripts_R/integracion/diablos.rds")
