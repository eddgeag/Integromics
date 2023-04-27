


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

semillas <- seq(40,80,10)

directorio_modelo <- "./scripts_R/integracion/modelo_optimo_comp"
if (!dir.exists(directorio_modelo)) {
  dir.create(directorio_modelo)
}
mofa_componentes <- function(ncomp,semilla){
  metaboloma <- (transformacion(datos$comunes$metaboloma))
  metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
  mofa.obj <- create_mofa_from_matrix(list(
    metaboloma = t(metaboloma),
    metagenoma = t(metagenoma)
  ))
  
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
ncomps <- 2:11
for(s in semillas){
  
  lapply(ncomps, function(x) mofa_componentes(ncomp = x,semilla = s))
  
}

modelos <- lapply(list.files(path = "./scripts_R/integracion/modelo_optimo_comp",
                             full.names = T), load_model)

w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

factores <- length(colnames(get_factors(modelo)[[1]]))

directorio_modelo_quasi <- "./scripts_R/integracion/quasi"

mofa_perm <- function(ncomp,semilla,dummy){
  metaboloma <- (transformacion(datos$comunes$metaboloma))
  idx <- sample(nrow(metaboloma))
  metaboloma <- metaboloma[idx,]
  metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
  metagenoma <- metagenoma[idx,]
  print(dim(metagenoma))
  print(dim(metaboloma))
  mofa.obj <- create_mofa_from_matrix(list(
    metaboloma = t(metaboloma),
    metagenoma = t(metagenoma)
  ))
  
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
  
  
  outfile <- file.path(directorio_modelo_quasi, paste0("model_",dummy,".hdf5"))
  
  MOFAobject.trained <- run_mofa(MOFAobject, outfile)
  
}

metaboloma <- (transformacion(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))

perms <- 1:1000
metaboloma.perm <- lapply(perms, function(x) mofa_perm(factores,modelo@training_options$seed,x))


modelos <- lapply(list.files("./scripts_R/integracion/quasi/",
                             full.names = T), load_model)

pesos.null <- lapply(modelos, function(x) get_weights(x,scale = T))

pesos.orig <- get_weights(modelo,scale = T)

p.valores <- lapply(pesos.null, function(x) apply(x,2,function(y) mean(abs()))

# ### Primero sera la comparacion tipo ggcorrplot
# 
# covariables <-
#   as.data.frame(datos$comunes$clinicos)
# 
# covariables <- covariables[, c("EDAD",
#                                "BMI",
#                                "WC",
#                                "WHR",
#                                "SHBG",
#                                "Free_ESTRA",
#                                "Total_ESTR",
#                                "hsCRP",
#                                "HOMAIRmean")]
# 
# 
# factores <- as.data.frame(get_factors(MOFAobject.trained)[[1]])
# 
# 
# 
# matrix_cov <-
#   matrix(NA, nrow = ncol(factores), ncol = ncol(covariables))
# for (j in 1:ncol(covariables)) {
#   variable <- covariables[,j]
# 
#     matrix_cov[, j] <-
#       apply(factores, 2, function(x)
#         corr.test(x, variable, method = "spearman")$p)
#     
#   
# }
# 
# colnames(matrix_cov) <- colnames(covariables)
# rownames(matrix_cov) <- colnames(factores)
# grupo <- datos$comunes$variables_in_bacteria$GROUP
# obesidad <- datos$comunes$variables_in_bacteria$OBESE
# resultado <- analyze_data(as.matrix((factores)),
#                           grupo=grupo,
#                           obesidad=obesidad,correccion = 2)$p_valores
# covariables_cor <- bind_cols(resultado,matrix_cov)
# covariables_cor.t <- as.data.frame(t(covariables_cor))
# 
# 
# 
# view2.factor1 <- plot_weights(MOFAobject.trained,return_data = T,view = 2,factors = 1)
