


## Preaparamos librerias CRAN
list.of.packages <-
  c(
    "xlsx",
    "kableExtra",
    "dplyr",
    "ggplot2",
    "egg",
    "cowplot",
    "patchwork",
    "gridExtra",
    "UsingR",
    "car",
    "lattice",
    "ggpubr",
    "ggbreak",
    "GGally",
    "reshape2",
    "ggcorrplot",
    "corrplot",
    "ggrepel",
    "factoextra",
    "chemometrics",
    "sparsepca",
    "vegan",
    "ca",
    "ade4",
    "pls",
    "OmicsPLS",
    "psych"
  )

### Preparamos librerias BIoconductor
list.of.bioc.packages <-
  c("phyloseq",
    "pcaMethods",
    "ropls",
    "mixOmics",
    "limma",
    "DESeq2")

### chequeamos que todo esté bien
res <- unlist(lapply(list.of.packages, require, character.only = T))
### chqueamos que todo este bien en bioconducot
res.bioc <-
  unlist(lapply(list.of.bioc.packages, require, character.only = T))

### Procesamos los datos de las bacterias

proces_bacteria <- function(datos, tipo) {
  Order <- bacteria_phylum.abs$Order
  variables_in_bacteria <-
    general_data[general_data$Orden %in% Order ,]
  Sample <- variables_in_bacteria$Paciente
  
  
  if (tipo == "genera") {
    Order <- Order
    
    X <- datos[-1, -1]
    X <- apply(X, 2, as.numeric)
    colnames(X) <- datos[1, 2:ncol(datos)]
    rownames(X) <- Sample
    X.rel <- 100 * X / rowSums(X)
    
    
  } else if (tipo == "phylum") {
    X <- datos[, -c(1:2)]
    colnames(X) <- colnames(datos)[3:ncol(datos)]
    rownames(X) <- Sample
    X.rel <- 100 * X / rowSums(X)
  }
  
  ###
  GROUP <- variables_in_bacteria$GROUP
  SEX <- variables_in_bacteria$SEX
  OBESE <- variables_in_bacteria$OBESE
  
  X <- as.data.frame(X)
  X.rel <- as.data.frame(X.rel)
  
  X$GROUP <- GROUP
  X$SEX <- SEX
  X$OBESE <- OBESE
  
  
  X.rel$GROUP <- GROUP
  X.rel$SEX <- SEX
  X.rel$OBESE <- OBESE
  
  return(list(abs = X, rel = X.rel))
  
}

### para quitar los bajos contajes
low.count.removal = function(data, # OTU count data frame of size n (sample) x p (OTU)
                             percent = 0.01 ) {
                             keep.otu = which(colSums(data) * 100 / (sum(colSums(data))) > percent)
                             data.filter = data[, keep.otu]
                             return(list(data.filter = data.filter, keep.otu = keep.otu))
                             }

low.count.removal.counts = function(data, # OTU count data frame of size n (sample) x p (OTU)
                                    abundance = 10) {
                                    keep.otu <- which(colSums(data) > abundance)
                                    data.filter <- data[, keep.otu]
                                    return(list(data.filter = data.filter, keep.otu = keep.otu))
                                    }

relative.fn <- function(data) {
  return(100 * data / rowSums(data))
  
}

### para extraer las bacterias comunes

common.bacteria <- function(datos_i, datos.comunes) {
  ## retorna los indices con las variables comunes
  
  interaccion <-
    interaction(datos.comunes$GROUP, datos.comunes$OBESE)
  datos.comunes$interaccion <- interaccion
  grupo <- levels(interaccion)
  w <- vector(mode = "list", length = 6)
  names(w) <- levels(interaccion)
  for (i in 1:6) {
    w[[names(w)[i]]] <- which(apply(as.matrix(datos_i[rownames(datos_i) %in%
                                                        datos.comunes[datos.comunes$interaccion == grupo[i], "Paciente"], ]), 2, function(x)
                                                          sum(x) == 0))
  }
  
  w2 <- vector(mode = "list", length = 6)
  names(w2) <- levels(interaccion)
  for (i in 1:6) {
    w2[[names(w2)[i]]] <- which(apply(as.matrix(datos_i[rownames(datos_i) %in%
                                                          datos.comunes[datos.comunes$interaccion == grupo[i], "Paciente"], ]), 2, function(x)
                                                            sum(x) > 0))
  }
  
  w3 <- lapply(w2, names)
  common_variables <- Reduce(intersect, x = w3)
  
  return(common_variables)
}


## leemos los datos clinicos

general_data <-
  read.xlsx("../../datos/Integromics_1.xlsx", sheetIndex = 1)
general_data$SEX <-
  factor(general_data$SEX,
         levels = c(0, 2),
         labels = c("Female", "Male"))
general_data$OBESE <-
  factor(general_data$OBESE,
         levels = c(0, 1),
         labels = c("No Obese", "Obese"))
general_data$GROUP <-
  factor(
    general_data$GROUP,
    levels = c(0, 1, 2),
    labels = c("Female",  "PCOS", "Male")
  )

##  leemos los datos IP

IP_file <-
  file.path("../../datos/Integromics_IPmarkers.xlsx") ## file path of the data
IP.tmp <- read.xlsx(IP_file, sheetIndex = 1)
IP <- IP.tmp[, -7] # remove the missin values column
IP$SEX <- factor(IP$SEX,
                 levels = c(0, 2),
                 labels = c("Female", "Male"))
IP$OBESE <-
  factor(IP$OBESE,
         levels = c(0, 1),
         labels = c("No Obese", "Obese"))
IP$GROUP <-
  factor(IP$GROUP,
         levels = c(0, 1, 2),
         labels = c("Female", "PCOS", "Male"))

### leemos el metaboloma
metabolome_file <-
  file.path("../../datos/Integromics_Metabolome.xlsx") ## file path of the data
metabolome.tmp <-
  read.xlsx(metabolome_file, sheetIndex = 1) # read the data
metabolites_name <-
  read.xlsx(metabolome_file, sheetIndex = 2)$Metabolitos # read the names of the metabolites
any(is.na(metabolome.tmp)) # there is a row of missing values.
metabolome <-
  metabolome.tmp[-nrow(metabolome.tmp), ] ## remove the missing value row
metabolome$GROUP <- general_data$GROUP ## add GROUP variable
metabolome$OBESE <- general_data$OBESE ## add OBESE variables.


# Preparando los datos para la integracion
bacteria_phylum.abs <- as.data.frame(read.xlsx(file.path("../../datos/integromics_microbiota.xlsx"), sheetIndex = 1))
bacteria_genera.abs <- as.data.frame(t(read.xlsx(file.path("../../datos/integromics_microbiota.xlsx"), sheetIndex = 3)))
phylum.list <- proces_bacteria(bacteria_phylum.abs,tipo="phylum")
genera.list <- proces_bacteria(bacteria_genera.abs,tipo="genera")

phylum.abs <- phylum.list$abs;phylum.rel<-phylum.list$rel
genera.abs <- genera.list$abs;genera.rel<-genera.list$rel


variables_in_bacteria <-
  general_data[general_data$Orden %in% bacteria_phylum.abs$Order,]



## Datos Clinicos
# **NOTE we are going to include also on both levels of integration body measures**
# We exclude ISI because it is measured after the load of glucose, the mean values, to exclude redundancy of the data, ant the postprandial levels.


sujetos <- variables_in_bacteria$Paciente
Clin.tmp <-
  variables_in_bacteria[, -c(1, 2, 3, 4, 5)] # get rid off orden, paciente, sex, group and obese
auc_clin <- colnames(Clin.tmp)[grep("AUC", colnames(Clin.tmp))]
mean_basal_clin <- colnames(Clin.tmp)[grep("_B$", colnames(Clin.tmp))]
sog <- colnames(Clin.tmp)[grep("^SO[GLP]", colnames(Clin.tmp))]
variables_clin <- colnames(Clin.tmp)
testosterona <- colnames(Clin.tmp)[grep("TEST", colnames(Clin.tmp))]
interr <- colnames(Clin.tmp)[grep("interaccion", colnames(Clin.tmp))]
ratio <- colnames(Clin.tmp)[grep("Ratio", colnames(Clin.tmp))]
w <-
  c(
    which(variables_clin %in% auc_clin),
    which(variables_clin %in% mean_basal_clin),
    which(variables_clin %in% sog),
    which(variables_clin %in% testosterona),
    which(variables_clin %in% interr),
    which(variables_clin %in% ratio)
  )
Clin <- Clin.tmp[, -w]
rownames(Clin) <- sujetos
colnames(Clin)



## Datos IP

IP.tmp <- IP[IP$Paciente %in% variables_in_bacteria$Paciente, ]
redundant_variables <-
  c("Orden", "Paciente", "SEX", "GROUP", "OBESE", "BMI")
variables_ip <- colnames(IP.tmp)
auc_ip <- variables_ip[grep("AUC", variables_ip)]
mean_ip <- variables_ip[grep("_0$", variables_ip)]
w <-
  c(
    which(variables_ip %in% redundant_variables),
    which(variables_ip %in% auc_ip),
    which(variables_ip %in% mean_ip)
  )
IP.basal <- IP.tmp[, -w]
rownames(IP.basal) <- sujetos



## Metabolome data (37)


metabolome.basal.tmp <-
  metabolome[metabolome$Order %in% variables_in_bacteria$Orden, ]
redundant_data_meta <- c("Order", "Sample", "GROUP", "OBESE")
variables_meta <- colnames(metabolome.basal.tmp)
mean_meta <- variables_meta[grep("_[GLP]0$", variables_meta)]
auc_meta <- variables_meta[grep("AU[CG]", variables_meta)]
w <-
  c(
    which(variables_meta %in% redundant_data_meta),
    which(variables_meta %in% mean_meta),
    which(variables_meta %in% auc_meta)
  )
metabolome.basal_mean <- metabolome.basal.tmp[, -w]
rownames(metabolome.basal_mean) <- sujetos
colnames(metabolome.basal_mean)

GRUPO <- variables_in_bacteria$GROUP
OBESIDAD <- variables_in_bacteria$OBESE 

## procesado a 94 genero
genera.abs.tmp <- genera.abs
redundant_data_ph <- c("GROUP","SEX","OBESE")
genus.tmp.tmp <- genera.abs.tmp[,-c(ncol(genera.abs.tmp)-2,ncol(genera.abs.tmp)-1,ncol(genera.abs.tmp))]
genus <- genus.tmp.tmp[,colSums(genus.tmp.tmp)>0]
rownames(genus) <- sujetos
## removal of genus unidentified

genus<-genus[,-grep("^ud-",colnames(genus))]

common_genus <-common.bacteria(genus,variables_in_bacteria)

genus<-genus[,common_genus]
# genus.abs <- low.count.removal.counts(genus,10)$data.filter


met37<-t(metabolome.basal_mean)
met37_norm <- met37/rowSums(met37)
met37_log <- t(apply(met37_norm,1,log))
design <- model.matrix(~OBESIDAD*GRUPO)
METABOLOMA<-t(voom(met37_norm,normalize.method = "quantile",design = design)$E)
library(ALDEx2)

GENERO <-aldex.clr(t(genus.abs))
library(vegan)
mc.instances <- numMCInstances(GENERO)
mc.all <- getMonteCarloInstances(GENERO)
obes.adonis <- vector("numeric",length = mc.instances)
grup.adonis <- vector("numeric",length = mc.instances)
inter.adonis <- vector("numeric",length = mc.instances)
  for(mc.i in 1:mc.instances) {
    t.input <- (sapply(mc.all, function(y) {
      y[, mc.i]
    }))
    res<-adonis2(t(t.input)~OBESIDAD*GRUPO,method = "euclidean")
    obes.adonis[mc.i] <- res$`Pr(>F)`[1]
    grup.adonis[mc.i] <- res$`Pr(>F)`[2]
    inter.adonis[mc.i] <-res$`Pr(>F)`[3]
  }

res.adonis <- cbind(obesidad=obes.adonis,grupo=grup.adonis,interaccion=inter.adonis)
colMeans(res.adonis)



# INTEGRACION
# 
# library(MOFA2)
# set.seed(123)
# 
# seed <- sample(1000,10)
# mc.instances <- numMCInstances(GENERO)
# mc.all <- getMonteCarloInstances(GENERO)
# i<-0
# directorioMOFA <-"./12-2-22/scripts_analysis_reunion/MOFA_analisis_1_94_mettrans"
# for(s in seed){
#   for(mc.i in 1:mc.instances) {
#     i<-i+1
#     t.input <- (sapply(mc.all, function(y) {
#       y[, mc.i]
#     }))
# 
#     genero <- t.input
#     datos <- list(metaboloma=(X),microbioma=(genero))
#     MOFAobject <-create_mofa_from_matrix(datos,groups=NULL)
#     samples_metadata(MOFAobject)<- data.frame(sample=pacientes,obesidad=(obesidad),grupo=(grupo),edad=Clin$EDAD,BMI=Clin$BMI,WC=Clin$WC,WHR=Clin$WHR,SHGB=Clin$SHBG,hsCRP=Clin$hsCRP,TE2=Clin$Total_ESTR,E2=Clin$Free_ESTRA,obesidad=as.numeric(obesidad),grupo=as.numeric(grupo))
#     rownames(samples_metadata(MOFAobject))<-pacientes
#     data_opts <- get_default_data_options(MOFAobject)
#     data_opts$scale_views <-T
#     model_opts <- get_default_model_options(MOFAobject)
#     model_opts$likelihoods=c("gaussian","gaussian")
#     model_opts$num_factors=10
#     head(model_opts)
#     model_opts$spikeslab_factors<-F
#     model_opts$ard_factors <-T
#     train_opts <- get_default_training_options(MOFAobject)
#     train_opts$seed <- s
#     # outfile = file.path(getwd(),paste0(seq(1:10),"_model.hdf5"))
#     MOFAobject <- prepare_mofa(
#       object = MOFAobject,
#       data_options = data_opts,
#       model_options = model_opts,
#       training_options = train_opts)
#     if(!dir.exists())
#     outfile <- file.path(getwd(),"modelsMOFA/",paste0("model_","seed_",s,"_instace_",mc.i,".hdf5"))
#     reticulate::use_python("/usr/bin/python3.10")
# 
#     MOFAobject <- run_mofa(MOFAobject, outfile,use_basilisk = F)
#     MOFAobject@covariates<-data.frame(obesidad=as.numeric(obesidad),grupo=as.numeric(grupo),edad=Clin$EDAD,BMI=Clin$BMI,WC=Clin$WC,WHR=Clin$WHR,SHGB=Clin$SHBG,hsCRP=Clin$hsCRP,TE2=Clin$Total_ESTR,E2=Clin$Free_ESTRA)
# 
#   }
# }
# 
# 
# 
