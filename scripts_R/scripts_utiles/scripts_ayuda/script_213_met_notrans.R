





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
    "psych",
    "lmPerm",
    "car"
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
    general_data[general_data$Orden %in% Order , ]
  Sample <- variables_in_bacteria$Paciente
  
  
  if (tipo == "genera") {
    Order <- Order
    
    X <- datos[-1,-1]
    X <- apply(X, 2, as.numeric)
    colnames(X) <- datos[1, 2:ncol(datos)]
    rownames(X) <- Sample
    X.rel <- 100 * X / rowSums(X)
    
    
  } else if (tipo == "phylum") {
    X <- datos[,-c(1:2)]
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
                             percent = 0.01) {
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
    w[[names(w)[i]]] <-
      which(apply(as.matrix(datos_i[rownames(datos_i) %in%
                                      datos.comunes[datos.comunes$interaccion == grupo[i], "Paciente"],]), 2, function(x)
                                        sum(x) == 0))
  }
  
  w2 <- vector(mode = "list", length = 6)
  names(w2) <- levels(interaccion)
  for (i in 1:6) {
    w2[[names(w2)[i]]] <-
      which(apply(as.matrix(datos_i[rownames(datos_i) %in%
                                      datos.comunes[datos.comunes$interaccion == grupo[i], "Paciente"],]), 2, function(x)
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
IP <- IP.tmp[,-7] # remove the missin values column
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
  metabolome.tmp[-nrow(metabolome.tmp),] ## remove the missing value row
metabolome$GROUP <- general_data$GROUP ## add GROUP variable
metabolome$OBESE <- general_data$OBESE ## add OBESE variables.


# Preparando los datos para la integracion
bacteria_phylum.abs <-
  as.data.frame(read.xlsx(
    file.path("../../datos/integromics_microbiota.xlsx"),
    sheetIndex = 1
  ))
bacteria_genera.abs <-
  as.data.frame(t(read.xlsx(
    file.path("../../datos/integromics_microbiota.xlsx"),
    sheetIndex = 3
  )))
phylum.list <- proces_bacteria(bacteria_phylum.abs, tipo = "phylum")
genera.list <- proces_bacteria(bacteria_genera.abs, tipo = "genera")

phylum.abs <- phylum.list$abs
phylum.rel <- phylum.list$rel
genera.abs <- genera.list$abs
genera.rel <- genera.list$rel


variables_in_bacteria <-
  general_data[general_data$Orden %in% bacteria_phylum.abs$Order, ]



## Datos Clinicos
# **NOTE we are going to include also on both levels of integration body measures**
# We exclude ISI because it is measured after the load of glucose, the mean values, to exclude redundancy of the data, ant the postprandial levels.


sujetos <- variables_in_bacteria$Paciente
Clin.tmp <-
  variables_in_bacteria[,-c(1, 2, 3, 4, 5)] # get rid off orden, paciente, sex, group and obese
auc_clin <- colnames(Clin.tmp)[grep("AUC", colnames(Clin.tmp))]
mean_basal_clin <-
  colnames(Clin.tmp)[grep("_B$", colnames(Clin.tmp))]
sog <- colnames(Clin.tmp)[grep("^SO[GLP]", colnames(Clin.tmp))]
variables_clin <- colnames(Clin.tmp)
testosterona <- colnames(Clin.tmp)[grep("TEST", colnames(Clin.tmp))]
interr <-
  colnames(Clin.tmp)[grep("interaccion", colnames(Clin.tmp))]
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
Clin <- Clin.tmp[,-w]
rownames(Clin) <- sujetos
colnames(Clin)



## Datos IP

IP.tmp <- IP[IP$Paciente %in% variables_in_bacteria$Paciente,]
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
IP.basal <- IP.tmp[,-w]
rownames(IP.basal) <- sujetos



###########3
###########
#############
############
#####METABOLOMA#############


metabolome.basal.tmp <- metabolome
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
metabolome.basal_mean <- metabolome.basal.tmp[,-w]
rownames(metabolome.basal_mean) <- general_data$Paciente

metabolito <- metabolome.basal_mean$Leu_mean
grupo <- general_data$GROUP
obesidad <- general_data$OBESE

library(car)
library(lmPerm)

## obesidad
metaboloma <- metabolome.basal_mean

if (ncol(metaboloma) != 37) {
  stop(errorCondition("No estan correctas las dimensiones"))
}
vector_dif_obesidad <- vector("numeric", length = ncol(metaboloma))
vector_pvalue_obesidad <-
  vector("numeric", length = ncol(metaboloma))
vector_prueba_obesidad <-
  vector("character", length = ncol(metaboloma))

for (i in 1:ncol(metaboloma)) {
  metabolito <- metaboloma[, i]
  dif_media <-
    mean(metabolito[obesidad == "No Obese"]) - mean(metabolito[obesidad == "Obese"])
  vector_dif_obesidad[i] <- dif_media
  mod <- lm(metabolito ~ obesidad)
  p.value_F <- leveneTest(mod)$`Pr(>F)`[1]
  p.value_norm <- shapiro.test(residuals(mod))$p.value
  if (p.value_F < 0.05 && p.value_norm > 0.05) {
    vector_pvalue_obesidad[i] <-
      t.test(metabolito ~ obesidad,
             var.equal = F,
             paired = F)$p.value
    vector_prueba_obesidad[i] <- "ttest"
  } else if (p.value_F > 0.05 && p.value_norm > 0.05) {
    vector_pvalue_obesidad[i] <-
      t.test(metabolito ~ obesidad,
             var.equal = T,
             paired = F)$p.value
    vector_prueba_obesidad[i] <- "ttest"
    
    
  } else if (p.value_F > 0.05 && p.value_norm < 0.05) {
    vector_pvalue_obesidad[i] <-
      wilcox.test(
        metabolito ~ obesidad,
        var.equal = T,
        paired = F,
        exact = F
      )$p.value
    vector_prueba_obesidad[i] <- "wilcox"
    
    
  } else if (p.value_F < 0.05 | p.value_norm < 0.05) {
    vector_pvalue_obesidad[i] <-
      wilcox.test(
        metabolito ~ obesidad,
        var.equal = T,
        paired = F,
        exact = F
      )$p.value
    vector_prueba_obesidad[i] <- "wilcox"
    
    
  }
  
  
}

### datos obesidad

res_obesidad <-
  data.frame(diferencia = vector_dif_obesidad, p_valor = vector_pvalue_obesidad)
res_obesidad$p.adj.value <- p.adjust(res_obesidad$p_valor, "BH")
rownames(res_obesidad) <- colnames(metaboloma)
res_obesidad$significante <-
  ifelse(res_obesidad$p.adj.value < 0.05, rownames(res_obesidad), NA)
res_obesidad$tipo_prueba <- vector_prueba_obesidad

### datos grupo
combinacion <- t(combn(levels(grupo), 2))
matriz_pvalores <- matrix(NA, nrow = ncol(metaboloma), ncol = 4)
matriz_location <- matrix(NA, nrow = ncol(metaboloma), ncol = 3)
vector_prueba_grupo<-
  vector("character", length = ncol(metaboloma))

for (i in 1:ncol(metaboloma)) {
  metabolito <- metaboloma[, i]
  mod_1 <- lm(metabolito ~ grupo)
  if (shapiro.test(residuals(mod_1))$p.value < 0.05) {
    matriz_pvalores[i, 4] <- kruskal.test(metabolito, grupo)$p.value
    vector_prueba_obesidad[i] <- "np"
  } else{
    matriz_pvalores[i, 4] <- Anova(mod_1)$`Pr(>F)`[1]
    vector_prueba_obesidad[i] <- "param"
    
  }
  
  for (j in 2:(nrow(combinacion) + 1)) {
    j <- j - 1
    variables_interes <- combinacion[j, ]
    mod <- lm(metabolito[grupo == variables_interes[1] |
                           grupo == variables_interes[2]] ~ grupo[grupo == variables_interes[1] |
                                                                    grupo == variables_interes[2]])
    p.value_F <- leveneTest(mod)$`Pr(>F)`[1]
    p.value_norm <- shapiro.test(residuals(mod))$p.value
    
    
    if (p.value_F < 0.05 && p.value_norm > 0.05) {
      matriz_pvalores[i, j] <-
        t.test(metabolito[grupo == variables_interes[1] |
                            grupo == variables_interes[2]] ~ grupo[grupo == variables_interes[1] |
                                                                     grupo == variables_interes[2]], var.equal =
                 F)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[grupo == variables_interes[1]]), mean(metabolito[grupo ==
                                                                                  variables_interes[2]])))
 
      
    } else if (p.value_F < 0.05 && p.value_norm < 0.05) {
      matriz_pvalores[i, j] <-
        wilcox.test(metabolito[grupo == variables_interes[1] |
                                 grupo == variables_interes[2]] ~ grupo[grupo == variables_interes[1] |
                                                                          grupo == variables_interes[2]], exact =
                      F)$p.value
      matriz_location[i, j] <-
        diff(c(median(metabolito[grupo == variables_interes[1]]), median(metabolito[grupo ==
                                                                                      variables_interes[2]])))
      
    } else if (p.value_F > 0.05 && p.value_norm > 0.05) {
      matriz_pvalores[i, j] <-
        t.test(metabolito[grupo == variables_interes[1] |
                            grupo == variables_interes[2]] ~ grupo[grupo == variables_interes[1] |
                                                                     grupo == variables_interes[2]], var.equal =
                 T)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[grupo == variables_interes[1]]), mean(metabolito[grupo ==
                                                                                  variables_interes[2]])))
      
    } else if (p.value_F > 0.05 && p.value_norm < 0.05) {
      matriz_pvalores[i, j] <-
        wilcox.test(metabolito[grupo == variables_interes[1] |
                                 grupo == variables_interes[2]] ~ grupo[grupo == variables_interes[1] |
                                                                          grupo == variables_interes[2]], exact =
                      F)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[grupo == variables_interes[1]]), mean(metabolito[grupo ==
                                                                                  variables_interes[2]])))
      
    }
    
    
  }
  
}

colnames(matriz_pvalores) <-
  c("Female-PCOS", "Female-Male", "PCOS-Male", "interaccion")
rownames(matriz_pvalores) <- colnames(metaboloma)
matriz_pvalores.adj <-
  as.data.frame(t(apply(matriz_pvalores, 1, function(x)
    c(p.adjust(x[1:3], "BH"), x[4]))))
colnames(matriz_location) <-
  c("mu Female-PCOS", "mu Female-Male", "mu PCOS-Male")
matriz_grupo <- cbind(matriz_location, matriz_pvalores.adj,vector_prueba_grupo)

### interaccion
vector_interaccion <-
  vector(mode = "numeric", length = ncol(metaboloma))
vector_obesidad <-
  vector(mode = "numeric", length = ncol(metaboloma))
vector_grupo <- vector(mode = "numeric", length = ncol(metaboloma))
aux.pairwise <- function(x, y) {
  #x variable dependiente el metabolito
  #y obesidad o grupo
  mod <- lm(x ~ y)
  p.homo <- leveneTest(mod)$`Pr(>F)`[1]
  p.norm <- shapiro.test(residuals(mod))$p.value
  if (p.homo > 0.05 && p.norm > 0.05) {
    media <-
      diff(c(mean(x[which(obesidad == levels(y)[1])], na.rm = T),
             mean(x[which(obesidad == levels(y)[2])], na.rm = T)))
    p.valor <- t.test(x ~ y, var.equal = T)$p.value
    param <- "ttest"
    
  } else if (p.homo < 0.05 && p.norm > 0.05) {
    media <-
      diff(c(mean(x[which(obesidad == levels(y)[1])], na.rm = T),
             mean(x[which(obesidad == levels(y)[2])], na.rm = T)))
    p.valor <- t.test(x ~ y, var.equal = F)$p.value
    param <- "ttest"
  } else if (p.norm < 0.05) {
    media <-
      diff(c(median(x[which(obesidad == levels(y)[1])], na.rm = T),
             median(x[which(obesidad == levels(y)[2])], na.rm = T)))
    p.valor <- wilcox.test(x ~ y, var.equal = F, exact = F)$p.value
    param <- "wilcox"
  }
  return(c(
    media = media,
    p.valor = p.valor,
    param = param
  ))
}

vector_obesidad <- vector("numeric", length = ncol(metaboloma))
vector_grupo <- vector("numeric", length = ncol(metaboloma))
vector_interaccion <- vector("numeric", length = ncol(metaboloma))
vector_pcos <- vector("numeric", length = ncol(metaboloma))
vector_female <- vector("numeric", length = ncol(metaboloma))
vector_male <- vector("numeric", length = ncol(metaboloma))
vector_pcos_media <- vector("numeric", length = ncol(metaboloma))
vector_female_media <- vector("numeric", length = ncol(metaboloma))
vector_male_media <- vector("numeric", length = ncol(metaboloma))
vector_male_prueba <- vector("numeric", length = ncol(metaboloma))
vector_female_prueba <- vector("numeric", length = ncol(metaboloma))
vector_pcos_prueba <- vector("numeric", length = ncol(metaboloma))


for (i in 1:ncol(metaboloma)) {
  metabolito <- metaboloma[, i]
  mod_1 <- lm(metabolito ~ grupo * obesidad)
  
  if (shapiro.test(residuals(mod_1))$p.value < 0.05) {
    s <- summary(aovp(metabolito ~ grupo * obesidad))[[1]]$`Pr(Prob)`
    vector_interaccion_i <- s[3]
    vector_obesidad_i <- s[2]
    vector_grupo_i <- s[1]
  } else{
    s <- Anova(mod_1)
    vector_interaccion_i <- s[3, 4]
    vector_obesidad_i <- s[2, 4]
    vector_grupo_i <- s[1, 4]
  }
  
  
  ## grupo|obesos
  PCOS <- metabolito[grupo == "PCOS"]
  PCOS.obesidad <- obesidad[grupo == "PCOS"]
  FEMALE <- metabolito[grupo == "Female"]
  FEMALE.obesidad <- obesidad[grupo == "Female"]
  MALE <- metabolito[grupo == "Male"]
  MALE.obesidad <- obesidad[grupo == "Male"]
  
  res.pcos <- aux.pairwise(PCOS, PCOS.obesidad)
  res.female <- aux.pairwise(FEMALE, FEMALE.obesidad)
  res.male <- aux.pairwise(MALE, MALE.obesidad)
  
  ## preparacion resultados
  vector_obesidad[i] <- vector_obesidad_i
  vector_grupo[i] <- vector_grupo_i
  vector_interaccion[i] <- vector_interaccion_i
  vector_pcos[i] <- res.pcos[2]
  vector_female[i] <- res.female[2]
  vector_male[i] <- res.male[2]
  vector_pcos_media[i] <- res.pcos[1]
  vector_female_media[i] <- res.female[1]
  vector_male_media[i] <- res.male[1]
  vector_pcos_prueba[i] <- res.pcos[3]
  vector_female_prueba[i] <- res.female[3]
  vector_male_prueba[i] <- res.male[3]
  
  
}

res.interaccion_grupo_obesidad <-
  data.frame(
    glm_obesidad = vector_obesidad,
    glm_grupo = vector_grupo,
    glm_interaccion = vector_interaccion,
    pcos = vector_pcos,
    pcos_media = vector_pcos_media,
    pcos_prueba = vector_pcos_prueba,
    female = vector_female,
    female_media = vector_female_media,
    female_prueba = vector_female_prueba,
    male = vector_male,
    male_media = vector_male_media,
    male_prueba = vector_male_prueba
  )


rownames(res.interaccion_grupo_obesidad) <- colnames(metaboloma)

res.interaccion_grupo_obesidad[, c("pcos", "female", "male")] <-
  apply(res.interaccion_grupo_obesidad[, c("pcos", "female", "male")], 2, function(x)
    p.adjust(x, "BH"))

### datos grupo NO OBESOS
combinacion <- t(combn(levels(grupo), 2))
matriz_pvalores <- matrix(NA, nrow = ncol(metaboloma), ncol = 4)
matriz_location <- matrix(NA, nrow = ncol(metaboloma), ncol = 3)

for (i in 1:ncol(metaboloma)) {
  metabolito <- metaboloma[obesidad == "No Obese", i]
  Y <- grupo[obesidad == "No Obese"]
  Y <- droplevels(Y)
  mod_1 <- lm(metabolito ~ Y)
  if (shapiro.test(residuals(mod_1))$p.value < 0.05) {
    matriz_pvalores[i, 4] <- kruskal.test(metabolito, Y)$p.value
  } else{
    matriz_pvalores[i, 4] <- Anova(mod_1)$`Pr(>F)`[1]
  }
  for (j in 2:(nrow(combinacion) + 1)) {
    j <- j - 1
    variables_interes <- combinacion[j, ]
    mod <- lm(metabolito[Y == variables_interes[1] |
                           Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                            Y == variables_interes[2]])
    p.value_F <- leveneTest(mod)$`Pr(>F)`[1]
    p.value_norm <- shapiro.test(residuals(mod))$p.value
    
    
    if (p.value_F < 0.05 && p.value_norm > 0.05) {
      matriz_pvalores[i, j] <-
        t.test(metabolito[Y == variables_interes[1] |
                            Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                             Y == variables_interes[2]], var.equal =
                 F)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[Y == variables_interes[1]]), mean(metabolito[Y ==
                                                                              variables_interes[2]])))
      
    } else if (p.value_F < 0.05 && p.value_norm < 0.05) {
      matriz_pvalores[i, j] <-
        wilcox.test(metabolito[Y == variables_interes[1] |
                                 Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                                  Y == variables_interes[2]], exact =
                      F)$p.value
      matriz_location[i, j] <-
        diff(c(median(metabolito[Y == variables_interes[1]]), median(metabolito[Y ==
                                                                                  variables_interes[2]])))
      
    } else if (p.value_F > 0.05 && p.value_norm > 0.05) {
      matriz_pvalores[i, j] <-
        t.test(metabolito[Y == variables_interes[1] |
                            Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                             Y == variables_interes[2]], var.equal =
                 T)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[Y == variables_interes[1]]), mean(metabolito[Y ==
                                                                              variables_interes[2]])))
      
    } else if (p.value_F > 0.05 && p.value_norm < 0.05) {
      matriz_pvalores[i, j] <-
        wilcox.test(metabolito[Y == variables_interes[1] |
                                 Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                                  Y == variables_interes[2]], exact =
                      F)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[Y == variables_interes[1]]), mean(metabolito[Y ==
                                                                              variables_interes[2]])))
      
    }
    
    
  }
  
}

colnames(matriz_pvalores) <-
  c("Female-PCOS", "Female-Male", "PCOS-Male", "interaccion")
rownames(matriz_pvalores) <- colnames(metaboloma)
matriz_pvalores.adj <-
  as.data.frame(t(apply(matriz_pvalores, 1, function(x)
    c(p.adjust(x[1:3], "BH"), x[4]))))
colnames(matriz_location) <-
  c("mu Female-PCOS", "mu Female-Male", "mu PCOS-Male")
matriz_grupo_no_obesos <-
  cbind(matriz_location, matriz_pvalores.adj)


### datos grupo  OBESOS
combinacion <- t(combn(levels(grupo), 2))
matriz_pvalores <- matrix(NA, nrow = ncol(metaboloma), ncol = 4)
matriz_location <- matrix(NA, nrow = ncol(metaboloma), ncol = 3)

for (i in 1:ncol(metaboloma)) {
  metabolito <- metaboloma[obesidad == "Obese", i]
  Y <- grupo[obesidad == "Obese"]
  Y <- droplevels(Y)
  mod_1 <- lm(metabolito ~ Y)
  if (shapiro.test(residuals(mod_1))$p.value < 0.05) {
    matriz_pvalores[i, 4] <- kruskal.test(metabolito, Y)$p.value
  } else{
    matriz_pvalores[i, 4] <- Anova(mod_1)$`Pr(>F)`[1]
  }
  for (j in 2:(nrow(combinacion) + 1)) {
    j <- j - 1
    variables_interes <- combinacion[j, ]
    mod <- lm(metabolito[Y == variables_interes[1] |
                           Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                            Y == variables_interes[2]])
    p.value_F <- leveneTest(mod)$`Pr(>F)`[1]
    p.value_norm <- shapiro.test(residuals(mod))$p.value
    
    
    if (p.value_F < 0.05 && p.value_norm > 0.05) {
      matriz_pvalores[i, j] <-
        t.test(metabolito[Y == variables_interes[1] |
                            Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                             Y == variables_interes[2]], var.equal =
                 F)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[Y == variables_interes[1]]), mean(metabolito[Y ==
                                                                              variables_interes[2]])))
      
    } else if (p.value_F < 0.05 && p.value_norm < 0.05) {
      matriz_pvalores[i, j] <-
        wilcox.test(metabolito[Y == variables_interes[1] |
                                 Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                                  Y == variables_interes[2]], exact =
                      F)$p.value
      matriz_location[i, j] <-
        diff(c(median(metabolito[Y == variables_interes[1]]), median(metabolito[Y ==
                                                                                  variables_interes[2]])))
      
    } else if (p.value_F > 0.05 && p.value_norm > 0.05) {
      matriz_pvalores[i, j] <-
        t.test(metabolito[Y == variables_interes[1] |
                            Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                             Y == variables_interes[2]], var.equal =
                 T)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[Y == variables_interes[1]]), mean(metabolito[Y ==
                                                                              variables_interes[2]])))
      
    } else if (p.value_F > 0.05 && p.value_norm < 0.05) {
      matriz_pvalores[i, j] <-
        wilcox.test(metabolito[Y == variables_interes[1] |
                                 Y == variables_interes[2]] ~ Y[Y == variables_interes[1] |
                                                                  Y == variables_interes[2]], exact =
                      F)$p.value
      matriz_location[i, j] <-
        diff(c(mean(metabolito[Y == variables_interes[1]]), mean(metabolito[Y ==
                                                                              variables_interes[2]])))
      
    }
    
    
  }
  
}

colnames(matriz_pvalores) <-
  c("Female-PCOS", "Female-Male", "PCOS-Male", "interaccion")
rownames(matriz_pvalores) <- colnames(metaboloma)
matriz_pvalores.adj <-
  as.data.frame(t(apply(matriz_pvalores, 1, function(x)
    c(p.adjust(x[1:3], "BH"), x[4]))))
colnames(matriz_location) <-
  c("mu Female-PCOS", "mu Female-Male", "mu PCOS-Male")
matriz_grupo_obesos <- cbind(matriz_location, matriz_pvalores.adj)



aromatic <- colnames(metaboloma)[c(2, 3, 4, 1, 6, 35, 36, 34)]
other <- colnames(metaboloma)[c(21, 9, 16, 26, 22, 18, 13, 29, 30, 27)]
aa <- colnames(metaboloma)[c(33, 15, 20, 19, 23, 25)]
carbo <-
  colnames(metaboloma)[!colnames(metaboloma) %in% c(aromatic, other, aa)]


### resultados metaboloma

resultados_metaboloma <-
  list(
    obesidad = res_obesidad,
    grupo = matriz_grupo,
    grupo.obesidad = res.interaccion_grupo_obesidad,
    grupo.no.obesos = matriz_grupo_no_obesos,
    grupo.obesos = matriz_grupo_obesos
  )


METABOLOMA.OBESIDAD <- resultados_metaboloma$obesidad
METABOLOMA.OBESIDAD$explicacion <-
  ifelse(
    METABOLOMA.OBESIDAD$diferencia > 0 &
      METABOLOMA.OBESIDAD$p.adj.value < 0.05 ,
    "MAS CONCENTRACION NO OBESOS",
    NA
  )
METABOLOMA.OBESIDAD$explicacion <-
  ifelse(
    METABOLOMA.OBESIDAD$diferencia < 0 &
      METABOLOMA.OBESIDAD$p.adj.value < 0.05 ,
    "MENOS CONCENTRACION NO OBESOS",
    METABOLOMA.OBESIDAD$explicacion
  )
METABOLOMA.OBESIDAD$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
METABOLOMA.OBESIDAD$grupo_meta <-
  ifelse(colnames(metaboloma) %in% other,
         "OTHER",
         METABOLOMA.OBESIDAD$grupo_meta)
METABOLOMA.OBESIDAD$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aa,
         "DERIVADOS AA",
         METABOLOMA.OBESIDAD$grupo_meta)
METABOLOMA.OBESIDAD$grupo_meta <-
  ifelse(
    colnames(metaboloma) %in% carbo,
    "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
    METABOLOMA.OBESIDAD$grupo_meta
  )




METABOLOMA.GRUPO <- resultados_metaboloma$grupo
METABOLOMA.GRUPO$explicacion_interaccion <- NA
METABOLOMA.GRUPO$EXPLICACION_PCOS_FEMALE <- NA
METABOLOMA.GRUPO$EXPLICACION_PCOS_MALE <- NA
METABOLOMA.GRUPO$EXPLICACION_SEX <- NA
METABOLOMA.GRUPO$EXPLICACION_UNICA <- NA

##########
for(j in 1:ncol(metaboloma)){
  
  if(METABOLOMA.GRUPO$interaccion[j]<0.05){
    METABOLOMA.GRUPO$explicacion_interaccion[j] <- "DIFERENCIAS INTRAGRUPO"
    if(METABOLOMA.GRUPO$`Female-PCOS`[j]<0.05){
      
      if(METABOLOMA.GRUPO$`mu Female-PCOS`[j]<0){
        METABOLOMA.GRUPO$EXPLICACION_PCOS_FEMALE[j] <- "CONCENTRACION MAYOR EN PCOS"
      }else{
        METABOLOMA.GRUPO$EXPLICACION_PCOS_FEMALE[j] <- "CONCENTRACION MENOR EN PCOS"
        
      }
    }else{
      METABOLOMA.GRUPO$EXPLICACION_PCOS_FEMALE[j] <-""
    }
    if(METABOLOMA.GRUPO$`Female-Male`[j]<0.05){
      
      if(METABOLOMA.GRUPO$`mu Female-Male`[j]<0){
        METABOLOMA.GRUPO$EXPLICACION_SEX[j] <- "CONCENTRACION MAYOR EN MALE"
      }else{
        METABOLOMA.GRUPO$EXPLICACION_SEX[j] <- "CONCENTRACION MENOR EN MALE"
        
      }
    }else{
      METABOLOMA.GRUPO$EXPLICACION_SEX[j] <-""
    }
    if(METABOLOMA.GRUPO$`PCOS-Male`[j]<0.05){
      
      if(METABOLOMA.GRUPO$`mu PCOS-Male`[j]<0){
        METABOLOMA.GRUPO$EXPLICACION_PCOS_MALE[j] <- "CONCENTRACION MAYOR EN MALE vs PCOS"
      }else{
        METABOLOMA.GRUPO$EXPLICACION_PCOS_MALE[j] <- "CONCENTRACION MENOR EN MALE vs PCOS"
        
      }
    }else{
      METABOLOMA.GRUPO$EXPLICACION_PCOS_MALE[j] <-""
    }
    if (METABOLOMA.GRUPO$`Female-PCOS`[j] < 0.05 &
        METABOLOMA.GRUPO$`Female-Male`[j] > 0.05 &
        METABOLOMA.GRUPO$`PCOS-Male`[j] > 0.05) {
      METABOLOMA.GRUPO$EXPLICACION_UNICA[j] <- "SOLAMENTE SE OBSERVA ENTRE PCOS Y MUJERES"
    }else if(
      METABOLOMA.GRUPO$`Female-PCOS`[j] > 0.05 &
      METABOLOMA.GRUPO$`Female-Male`[j] < 0.05 &
      METABOLOMA.GRUPO$`PCOS-Male`[j] > 0.05
    ){
      METABOLOMA.GRUPO$EXPLICACION_UNICA[j] <- "SOLAMENTE SE OBSERVA ENTRE SEXOS"
      
    }else if(
      METABOLOMA.GRUPO$`Female-PCOS`[j] > 0.05 &
      METABOLOMA.GRUPO$`Female-Male`[j] > 0.05 &
      METABOLOMA.GRUPO$`PCOS-Male`[j] < 0.05
    ){
      METABOLOMA.GRUPO$EXPLICACION_UNICA[j] <- "SOLAMENTE SE OBSERVA ENTRE PCOS y HOMBRES"
      
    }else{
      METABOLOMA.GRUPO$EXPLICACION_UNICA[j] <- ""
      
    }
    
    
  }else{
    METABOLOMA.GRUPO$explicacion_interaccion[j] <- ""
    
  }
  
}

##########


METABOLOMA.GRUPO$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
METABOLOMA.GRUPO$grupo_meta <-
  ifelse(colnames(metaboloma) %in% other,
         "OTHER",
         METABOLOMA.OBESIDAD$grupo_meta)
METABOLOMA.GRUPO$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aa,
         "DERIVADOS AA",
         METABOLOMA.OBESIDAD$grupo_meta)
METABOLOMA.GRUPO$grupo_meta <-
  ifelse(
    colnames(metaboloma) %in% aromatic,
    "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
    METABOLOMA.OBESIDAD$grupo_meta
  )

METABOLOMA.grupo.obesidad <- resultados_metaboloma$grupo.obesidad

METABOLOMA.grupo.obesidad$explicacion_OBESIDAD <- NA
METABOLOMA.grupo.obesidad$explicacion_GRUPO <- NA
METABOLOMA.grupo.obesidad$explicacion_INTERACCION <- NA

METABOLOMA.grupo.obesidad$explicacion_PCOS_MEDIA <- NA
METABOLOMA.grupo.obesidad$explicacion_FEMALE_MEDIA <- NA
METABOLOMA.grupo.obesidad$explicacion_MALE_MEDIA <- NA


METABOLOMA.grupo.obesidad$explicacion_UNICA <- NA




for (j in 1:ncol(metaboloma)) {
  res_i <- METABOLOMA.grupo.obesidad[j, ]
  if (res_i$glm_obesidad < 0.05) {
    res_i$explicacion_OBESIDAD <- "DIFERENCIA ENTRE OBESOS"
  } else{
    res_i$explicacion_OBESIDAD <- ""
    
  }
  
  if (res_i$glm_grupo < 0.05) {
    res_i$explicacion_GRUPO <- "DIFERENCIA ENTRE GRUPOS"
  } else{
    res_i$explicacion_GRUPO <- ""
    
  }
  
  if (res_i$glm_interaccion < 0.05) {
    res_i$explicacion_INTERACCION <- "DIFERENCIA EN LA INTERACCION"
  } else{
    res_i$explicacion_INTERACCION <- ""
    
  }
  
  if (res_i$pcos < 0.05) {
    if (res_i$pcos_media < 0) {
      res_i$explicacion_PCOS_MEDIA <- "OBESIDAD MAYOR EN PCOS "
    } else{
      res_i$explicacion_PCOS_MEDIA <- "OBESIDAD MENOR EN PCOS"
    }
  } else{
    res_i$explicacion_PCOS_MEDIA <- ""
    
  }
  if (res_i$female < 0.05) {
    if (res_i$female_media < 0) {
      res_i$explicacion_FEMALE_MEDIA <- "OBESIDAD MAYOR EN FEMALE "
    } else{
      res_i$explicacion_FEMALE_MEDIA <- "OBESIDAD MENOR EN FEMALE"
    }
  } else{
    res_i$explicacion_FEMALE_MEDIA <- ""
    
  }
  if (res_i$male < 0.05) {
    if (res_i$male_media < 0) {
      res_i$explicacion_MALE_MEDIA <- "OBESIDAD MAYOR EN MALE "
    } else{
      res_i$explicacion_MALE_MEDIA <- "OBESIDAD MENOR EN MALE"
    }
  } else{
    res_i$explicacion_MALE_MEDIA <- ""
    
  }
  
  if (res_i$female < 0.05 && res_i$male > 0.05 && res_i$pcos > 0.05) {
    res_i$explicacion_UNICA <- "LA OBESIDAD SOLAMENtE SE VE EN MUJERES"
    
  } else  if (res_i$female > 0.05 &&
              res_i$male < 0.05 && res_i$pcos > 0.05) {
    res_i$explicacion_UNICA <- "LA OBESIDAD SOLAMENtE SE VE EN HOMBRES"
    
    
  } else  if (res_i$female > 0.05 &&
              res_i$male > 0.05 && res_i$pcos < 0.05) {
    res_i$explicacion_UNICA <- "LA OBESIDAD SOLAMENtE SE VE EN PCOS"
    
  } else{
    res_i$explicacion_UNICA <- ""
    
  }
  
  
  METABOLOMA.grupo.obesidad$explicacion_OBESIDAD[j] <-
    res_i$explicacion_OBESIDAD
  METABOLOMA.grupo.obesidad$explicacion_GRUPO[j] <-
    res_i$explicacion_GRUPO
  METABOLOMA.grupo.obesidad$explicacion_INTERACCION[j] <-
    res_i$explicacion_INTERACCION
  
  
  METABOLOMA.grupo.obesidad$explicacion_FEMALE_MEDIA[j] <-
    res_i$explicacion_FEMALE_MEDIA
  METABOLOMA.grupo.obesidad$explicacion_MALE_MEDIA[j] <-
    res_i$explicacion_MALE_MEDIA
  METABOLOMA.grupo.obesidad$explicacion_PCOS_MEDIA[j] <-
    res_i$explicacion_PCOS_MEDIA
  METABOLOMA.grupo.obesidad$explicacion_UNICA[j] <-
    res_i$explicacion_UNICA
  
}
METABOLOMA.grupo.obesidad$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
METABOLOMA.grupo.obesidad$grupo_meta <-
  ifelse(colnames(metaboloma) %in% other,
         "OTHER",
         METABOLOMA.grupo.obesidad$grupo_meta)
METABOLOMA.grupo.obesidad$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aa,
         "DERIVADOS AA",
         METABOLOMA.grupo.obesidad$grupo_meta)
METABOLOMA.grupo.obesidad$grupo_meta <-
  ifelse(
    colnames(metaboloma) %in% carbo,
    "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
    METABOLOMA.grupo.obesidad$grupo_meta
  )


#### obesidad por grupo

#### obesos
METABOLOMA.grupo_OBESOS <- resultados_metaboloma$grupo.obesos

METABOLOMA.grupo_OBESOS$INTERPRETACION_SEX <- NA
METABOLOMA.grupo_OBESOS$INTERPRETACION_FEMALE_PCOS <- NA
METABOLOMA.grupo_OBESOS$INTERPRETACION_MALE_PCOS <- NA
METABOLOMA.grupo_OBESOS$INTERPRETACION_DIFERENCIAS <- NA
METABOLOMA.grupo_OBESOS$INTERPRETACION_UNICA <- NA


for (j in 1:ncol(metaboloma)) {
  res_i <- METABOLOMA.grupo_OBESOS[j, ]
  if (res_i$interaccion < 0.05) {
    res_i$INTERPRETACION_DIFERENCIAS <- "HAY DIFERENCIAS INTRAGRUPO"
    if (res_i$`Female-PCOS` < 0.05) {
      if (res_i$`mu Female-PCOS` < 0) {
        res_i$INTERPRETACION_FEMALE_PCOS <-
          "LA CONCENTRACION ES MAYOR EN PCOS VS MUJERES"
        
      } else{
        res_i$INTERPRETACION_FEMALE_PCOS <-
          "LA CONCENTRACION ES MENO EN PCOS VS MUJERES"
        
      }
    } else{
      res_i$INTERPRETACION_FEMALE_PCOS <- ""
      
    }
    
    if (res_i$`Female-Male` < 0.05) {
      if (res_i$`mu Female-Male` < 0) {
        res_i$INTERPRETACION_SEX <-
          "LA CONCENTRACION ES MAYOR EN HOMBRES VS MUJERES"
        
      } else{
        res_i$INTERPRETACION_SEX <-
          "LA CONCENTRACION ES MENO EN HOMBRES VS MUJERES"
        
      }
    } else{
      res_i$INTERPRETACION_SEX <- ""
      
    }
    if (res_i$`PCOS-Male` < 0.05) {
      if (res_i$`mu PCOS-Male` < 0) {
        res_i$INTERPRETACION_MALE_PCOS <-
          "LA CONCENTRACION ES MAYOR EN HOMBRES VS PCOS"
        
      } else{
        res_i$INTERPRETACION_MALE_PCOS <-
          "LA CONCENTRACION ES MENO EN HOMBRES VS PCOS"
        
      }
    } else{
      res_i$INTERPRETACION_MALE_PCOS <- ""
      
    }
    
    if (res_i$`Female-PCOS` < 0.05 &&
        res_i$`Female-Male` > 0.05 && res_i$`mu PCOS-Male` > 0.05) {
      res_i$INTERPRETACION_UNICA <- "ENTRE MUJERES Y PCOS"
      
    } else if (res_i$`Female-PCOS` > 0.05 &&
               res_i$`Female-Male` < 0.05 && res_i$`mu PCOS-Male` > 0.05) {
      res_i$INTERPRETACION_UNICA <- "ENTRE SEXO"
      
      
    } else if (res_i$`Female-PCOS` > 0.05 &&
               res_i$`Female-Male` > 0.05 && res_i$`mu PCOS-Male` < 0.05) {
      res_i$INTERPRETACION_UNICA <- "ENTRE PCOS Y HOMBRES"
      
      
    } else{
      res_i$INTERPRETACION_UNICA <- ""
      
    }
    
  } else{
    res_i$INTERPRETACION_DIFERENCIAS <- ""
    
  }
  
  
  METABOLOMA.grupo_OBESOS$INTERPRETACION_DIFERENCIAS[j] <-
    res_i$INTERPRETACION_DIFERENCIAS
  METABOLOMA.grupo_OBESOS$INTERPRETACION_MALE_PCOS[j] <-
    res_i$INTERPRETACION_MALE_PCOS
  METABOLOMA.grupo_OBESOS$INTERPRETACION_SEX[j] <-
    res_i$INTERPRETACION_SEX
  METABOLOMA.grupo_OBESOS$INTERPRETACION_FEMALE_PCOS[j] <-
    res_i$INTERPRETACION_FEMALE_PCOS
  METABOLOMA.grupo_OBESOS$INTERPRETACION_UNICA[j] <-
    res_i$INTERPRETACION_UNICA
  
}

METABOLOMA.grupo_OBESOS$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
METABOLOMA.grupo_OBESOS$grupo_meta <-
  ifelse(colnames(metaboloma) %in% other,
         "OTHER",
         METABOLOMA.grupo_OBESOS$grupo_meta)
METABOLOMA.grupo_OBESOS$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aa,
         "DERIVADOS AA",
         METABOLOMA.grupo_OBESOS$grupo_meta)
METABOLOMA.grupo_OBESOS$grupo_meta <-
  ifelse(
    colnames(metaboloma) %in% carbo,
    "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
    METABOLOMA.grupo_OBESOS$grupo_meta
  )



### no obesis
METABOLOMA.grupo_NO_OBESOS <- resultados_metaboloma$grupo.no.obesos

METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_SEX <- NA
METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_FEMALE_PCOS <- NA
METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_MALE_PCOS <- NA
METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_DIFERENCIAS <- NA
METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_UNICA <- NA


for (j in 1:ncol(metaboloma)) {
  res_i <- METABOLOMA.grupo_NO_OBESOS[j, ]
  if (res_i$interaccion < 0.05) {
    res_i$INTERPRETACION_DIFERENCIAS <- "HAY DIFERENCIAS INTRAGRUPO"
    if (res_i$`Female-PCOS` < 0.05) {
      if (res_i$`mu Female-PCOS` < 0) {
        res_i$INTERPRETACION_FEMALE_PCOS <-
          "LA CONCENTRACION ES MAYOR EN PCOS VS MUJERES"
        
      } else{
        res_i$INTERPRETACION_FEMALE_PCOS <-
          "LA CONCENTRACION ES MENO EN PCOS VS MUJERES"
        
      }
    } else{
      res_i$INTERPRETACION_FEMALE_PCOS <- ""
      
    }
    
    if (res_i$`Female-Male` < 0.05) {
      if (res_i$`mu Female-Male` < 0) {
        res_i$INTERPRETACION_SEX <-
          "LA CONCENTRACION ES MAYOR EN HOMBRES VS MUJERES"
        
      } else{
        res_i$INTERPRETACION_SEX <-
          "LA CONCENTRACION ES MENO EN HOMBRES VS MUJERES"
        
      }
    } else{
      res_i$INTERPRETACION_SEX <- ""
      
    }
    if (res_i$`PCOS-Male` < 0.05) {
      if (res_i$`mu PCOS-Male` < 0) {
        res_i$INTERPRETACION_MALE_PCOS <-
          "LA CONCENTRACION ES MAYOR EN HOMBRES VS PCOS"
        
      } else{
        res_i$INTERPRETACION_MALE_PCOS <-
          "LA CONCENTRACION ES MENO EN HOMBRES VS PCOS"
        
      }
    } else{
      res_i$INTERPRETACION_MALE_PCOS <- ""
      
    }
    
    if (res_i$`Female-PCOS` < 0.05 &&
        res_i$`Female-Male` > 0.05 && res_i$`mu PCOS-Male` > 0.05) {
      res_i$INTERPRETACION_UNICA <- "ENTRE MUJERES Y PCOS"
      
    } else if (res_i$`Female-PCOS` > 0.05 &&
               res_i$`Female-Male` < 0.05 && res_i$`mu PCOS-Male` > 0.05) {
      res_i$INTERPRETACION_UNICA <- "ENTRE SEXO"
      
      
    } else if (res_i$`Female-PCOS` > 0.05 &&
               res_i$`Female-Male` > 0.05 && res_i$`mu PCOS-Male` < 0.05) {
      res_i$INTERPRETACION_UNICA <- "ENTRE PCOS Y HOMBRES"
      
      
    } else{
      res_i$INTERPRETACION_UNICA <- ""
      
    }
    
  } else{
    res_i$INTERPRETACION_DIFERENCIAS <- ""
    
  }
  
  
  METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_DIFERENCIAS[j] <-
    res_i$INTERPRETACION_DIFERENCIAS
  METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_MALE_PCOS[j] <-
    res_i$INTERPRETACION_MALE_PCOS
  METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_SEX[j] <-
    res_i$INTERPRETACION_SEX
  METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_FEMALE_PCOS[j] <-
    res_i$INTERPRETACION_FEMALE_PCOS
  METABOLOMA.grupo_NO_OBESOS$INTERPRETACION_UNICA[j] <-
    res_i$INTERPRETACION_UNICA
  
}

METABOLOMA.grupo_NO_OBESOS$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
METABOLOMA.grupo_NO_OBESOS$grupo_meta <-
  ifelse(colnames(metaboloma) %in% other,
         "OTHER",
         METABOLOMA.grupo_NO_OBESOS$grupo_meta)
METABOLOMA.grupo_NO_OBESOS$grupo_meta <-
  ifelse(colnames(metaboloma) %in% aa,
         "DERIVADOS AA",
         METABOLOMA.grupo_NO_OBESOS$grupo_meta)
METABOLOMA.grupo_NO_OBESOS$grupo_meta <-
  ifelse(
    colnames(metaboloma) %in% carbo,
    "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
    METABOLOMA.grupo_NO_OBESOS$grupo_meta
  )


##### RESULTADOS

METABOLOMA_RESULTADOS <- list(METABOLOMA.OBESIDAD=METABOLOMA.OBESIDAD,
                              METABOLOMA.GRUPO=METABOLOMA.GRUPO,
                              METABOLOMA.grupo.obesidad=METABOLOMA.grupo.obesidad,
                              METABOLOMA.grupo_NO_OBESOS=METABOLOMA.grupo_NO_OBESOS,
                              METABOLOMA.grupo_OBESOS=METABOLOMA.grupo_OBESOS)


##### PREPARACION GENERO #######
genus.abs<-genera.abs[,-c(grep("^ud-",colnames(genera.abs)),
                      grep("GROUP",colnames(genera.abs)),
                      grep("SEX",colnames(genera.abs)),
                      grep("OBESE",colnames(genera.abs)))]


phylum.abs<-phylum.abs[,-c(grep("^ud.",colnames(phylum.abs)),
                          grep("GROUP",colnames(phylum.abs)),
                          grep("SEX",colnames(phylum.abs)),
                          grep("OBESE",colnames(phylum.abs)))]
phylum.abs <-phylum.abs[, colSums(phylum.abs != 0) > 0]

## variables en comun

obesidad <- variables_in_bacteria$OBESE
grupo <- variables_in_bacteria$GROUP

library(ALDEx2)

phylum.aldex <- aldex.clr(t(phylum.abs))


mc.instances <- numMCInstances(phylum.aldex)
mc.all <- getMonteCarloInstances(phylum.aldex)
s.obesidad <- vector("numeric",length=mc.instances)
s.grupo <-  vector("numeric",length=mc.instances)
s.interaccion<-  vector("numeric",length=mc.instances)
p.val.obes <- vector("numeric",length=mc.instances)
loc.obes <- vector("numeric",length=mc.instances)
p.val.Sex <- vector("numeric",length=mc.instances)
p.val.Fms <- vector("numeric",length=mc.instances)
p.val.Hombs <- vector("numeric",length=mc.instances)
loc.Sex<-vector("numeric",length=mc.instances)
loc.Fms <-vector("numeric",length=mc.instances)
loc.Hombs <-vector("numeric",length=mc.instances)
p.val.Sex_ob <- vector("numeric",length=mc.instances)
p.val.Fms_ob <- vector("numeric",length=mc.instances)
p.val.Hombs_ob<- vector("numeric",length=mc.instances)
loc.Sex_ob<-vector("numeric",length=mc.instances)
loc.Fms_ob <-vector("numeric",length=mc.instances)
loc.Hombs_ob <-vector("numeric",length=mc.instances)

p.val.Sex_noob <- vector("numeric",length=mc.instances)
p.val.Fms_noob <- vector("numeric",length=mc.instances)
p.val.Hombs_noob <- vector("numeric",length=mc.instances)
loc.Sex_noob<-vector("numeric",length=mc.instances)
loc.Fms_noob <-vector("numeric",length=mc.instances)
loc.Hombs_noob <-vector("numeric",length=mc.instances)

s.obesidad2 <- vector("numeric",length=mc.instances)
s.grupo2 <-  vector("numeric",length=mc.instances)
s.interaccion2<-  vector("numeric",length=mc.instances)
p.val.obes2 <- vector("numeric",length=mc.instances)
loc.obes2 <- vector("numeric",length=mc.instances)
p.val.Sex2 <- vector("numeric",length=mc.instances)
p.val.Fms2 <- vector("numeric",length=mc.instances)
p.val.Hombs2 <- vector("numeric",length=mc.instances)
loc.Sex2<-vector("numeric",length=mc.instances)
loc.Fms2 <-vector("numeric",length=mc.instances)
loc.Hombs2 <-vector("numeric",length=mc.instances)
p.val.Sex_ob2 <- vector("numeric",length=mc.instances)
p.val.Fms_ob2 <- vector("numeric",length=mc.instances)
p.val.Hombs_ob2<- vector("numeric",length=mc.instances)
loc.Sex_ob2<-vector("numeric",length=mc.instances)
loc.Fms_ob2 <-vector("numeric",length=mc.instances)
loc.Hombs_ob2 <-vector("numeric",length=mc.instances)

p.val.Sex2_noob <- vector("numeric",length=mc.instances)
p.val.Fms2_noob <- vector("numeric",length=mc.instances)
p.val.Hombs2_noob <- vector("numeric",length=mc.instances)
loc.Sex2_noob<-vector("numeric",length=mc.instances)
loc.Fms2_noob <-vector("numeric",length=mc.instances)
loc.Hombs2_noob <-vector("numeric",length=mc.instances)
s.prueba <- vector("character",length=mc.instances)
s.prueba2 <- vector("character",length=mc.instances)

  for(mc.i in 1:mc.instances) {
    
    mc.i<-1
    t.input <- t(sapply(mc.all, function(y) {
      y[, mc.i]
    }))
    
      for(j in 1:ncol(t.input)){
        
        ## vemos grupo*obesidad
        bac_i <- t.input[,1]
        mod <- lm(bac_i ~ grupo*obesidad,contrasts = list(grupo=contr.sum,obesidad=contr.sum))
        p.norm <- shapiro.test(residuals(mod))$p.value
        
        ### grupo*obesos
        if(p.norm <0.05){
          
          s <- summary(aovp(bac_i ~ grupo*obesidad,contrasts = list(grupo=contr.sum,obesidad=contr.sum)))
          s.obesidad[j] <- s[[1]]$`Pr(Prob)`[2]
          s.grupo[j] <- s[[1]]$`Pr(Prob)`[1]
          s.interaccion[j] <- s[[1]]$`Pr(Prob)`[3]
          s.prueba[j] <- "noparam"
          
        }else{
          
          s <- Anova(mod,type=3)
          s.grupo[j] <- s$`Pr(>F)`[2]
          s.obesidad[j]  <-s$`Pr(>F)`[3]
          s.interaccion[j] <- s$`Pr(>F)`[4]
          s.prueba[j] <- "param"
          
        }
        ##vemos obesidad
        
        if(s.prueba=="param" && s.obesidad<0.05){
          
          mod <- lm(bac_i~obesidad)
          modhet <- leveneTest(mod)$`Pr(>F)`
          if(modhet<0.05){
            p.val.obes[j] <- t.test(bac_i ~obesidad,var.equal=F)$p.value
            loc.obes[j] <-diff(c(mean(bac_i[obesidad=="No Obese"]),mean(bac_i[obesiad=="Obese"])))
            
          }else{
            p.val.obes[j] <- t.test(bac_i ~obesidad,var.equal=F)$p.value
            loc.obes[j] <-diff(c(mean(bac_i[obesidad=="No Obese"]),mean(bac_i[obesiad=="Obese"])))
            
          }
          
          
        }else if(s.prueba=="noparam" && s.obesidad<0.05){
          p.val.obes[j] <- wilcox.test(bac_i ~obesidad,exact=F)$p.value
          loc.obes[j] <-diff(c(median(bac_i[obesidad=="No Obese"]),median(bac_i[obesiad=="Obese"])))
          
        }
        ##  grupo
        if(s.prueba=="param" && s.grupo < 0.05){
          
          mod <- lm(bac_i~grupo)
          modhet <- leveneTest(mod)$`Pr(>F)`
          if(modhet<0.05){
            p.val.Sex[j] <- t.test(bac_i[grupo!="PCOS"] ~grupo[grupo!="PCOS"],var.equal=F)$p.value
            p.val.Fms[j] <- t.test(bac_i[grupo!="Male"] ~grupo[grupo!="Male"],var.equal=F)$p.value
            p.val.Hombs[j] <- t.test(bac_i[grupo!="Female"] ~grupo[grupo!="Female"],var.equal=F)$p.value
            loc.Sex[j]<-diff(c(mean(bac_i[grupo!="PCOS"]),mean(bac_i[grupo!="PCOS"])))
            loc.Fms[j] <-diff(c(mean(bac_i[grupo!="Male"]),mean(bac_i[grupo!="Male"])))
            loc.Hombs[j] <-diff(c(mean(bac_i[grupo!="Female"]),mean(bac_i[grupo!="Female"])))
            

          }else{
            p.val.Sex[j] <- t.test(bac_i[grupo!="PCOS"] ~grupo[grupo!="PCOS"],var.equal=T)$p.value
            p.val.Fms[j] <- t.test(bac_i[grupo!="Male"] ~grupo[grupo!="Male"],var.equal=T)$p.value
            p.val.Hombs[j] <- t.test(bac_i[grupo!="Female"] ~grupo[grupo!="Female"],var.equal=T)$p.value
            loc.Sex[j]<-diff(c(mean(bac_i[grupo!="PCOS"]),mean(bac_i[grupo!="PCOS"])))
            loc.Fms[j] <-diff(c(mean(bac_i[grupo!="Male"]),mean(bac_i[grupo!="Male"])))
            loc.Hombs[j] <-diff(c(mean(bac_i[grupo!="Female"]),mean(bac_i[grupo!="Female"])))
            
            
          }
          
          
        }else if(s.prueba=="noparam" && s.grupo < 0.05){
          p.val.Sex[j] <- wilcox.test(bac_i[grupo!="PCOS"] ~grupo[grupo!="PCOS"],exact=F)$p.value
          p.val.Fms[j] <- wilcox.test(bac_i[grupo!="Male"] ~grupo[grupo!="Male"],exact=F)$p.value
          p.val.Hombs[j] <- wilcox.test(bac_i[grupo!="Female"] ~grupo[grupo!="Female"],exact=F)$p.value
          loc.Sex[j]<-diff(c(median(bac_i[grupo!="PCOS"]),median(bac_i[grupo!="PCOS"])))
          loc.Fms[j] <-diff(c(median(bac_i[grupo!="Male"]),median(bac_i[grupo!="Male"])))
          loc.Hombs[j] <-diff(c(median(bac_i[grupo!="Female"]),median(bac_i[grupo!="Female"])))
           }
        ### interaccion
        if(s.interaccion<0.05 ){
          
          NoObesos <- grupo[obesidad=="No Obese"]
          bac_i_noob <- bac_i[obesidad=="No Obese"]
          
          if(s.prueba=="param" && s.interaccion < 0.05){
            
            mod <- lm(bac_i~NoObesos)
            modhet <- leveneTest(mod)$`Pr(>F)`
            if(modhet<0.05){
              p.val.Sex_noob[j] <- t.test(bac_i_noob[NoObesos!="PCOS"] ~NoObesos[NoObesos!="PCOS"],var.equal=F)$p.value
              p.val.Fms_noob[j] <- t.test(bac_i_noob[NoObesos!="Male"] ~NoObesos[NoObesos!="Male"],var.equal=F)$p.value
              p.val.Hombs_noob[j] <- t.test(bac_i_noob[NoObesos!="Female"] ~NoObesos[NoObesos!="Female"],var.equal=F)$p.value
              loc.Sex_noob[j]<-diff(c(mean(bac_i_noob[NoObesos!="PCOS"]),mean(bac_i_noob[NoObesos!="PCOS"])))
              loc.Fms_noob[j] <-diff(c(mean(bac_i_noob[NoObesos!="Male"]),mean(bac_i_noob[NoObesos!="Male"])))
              loc.Hombs_noob[j] <-diff(c(mean(bac_i_noob[NoObesos!="Female"]),mean(bac_i_noob[NoObesos!="Female"])))
              
              
            }else{
              p.val.Sex_noob[j] <- t.test(bac_i_noob[NoObesos!="PCOS"] ~NoObesos[NoObesos!="PCOS"],var.equal=T)$p.value
              p.val.Fms_noob[j] <- t.test(bac_i_noob[NoObesos!="Male"] ~NoObesos[NoObesos!="Male"],var.equal=T)$p.value
              p.val.Hombs_noob[j] <- t.test(bac_i_noob[NoObesos!="Female"] ~NoObesos[NoObesos!="Female"],var.equal=T)$p.value
              loc.Sex_noob[j]<-diff(c(mean(bac_i_noob[NoObesos!="PCOS"]),mean(bac_i_noob[NoObesos!="PCOS"])))
              loc.Fms_noob[j] <-diff(c(mean(bac_i_noob[NoObesos!="Male"]),mean(bac_i_noob[NoObesos!="Male"])))
              loc.Hombs_noob[j] <-diff(c(mean(bac_i_noob[NoObesos!="Female"]),mean(bac_ibac_i_noob)))
              
              
            }
            
            
          }else if(s.prueba=="noparam" && s.interaccion <0.05 ){
            p.val.Sex <- wilcox.test(bac_i_noob[NoObesos!="PCOS"] ~NoObesos[NoObesos!="PCOS"],exact=F)$p.value
            p.val.Fms <- wilcox.test(bac_i_noob[NoObesos!="Male"] ~NoObesos[NoObesos!="Male"],exact=F)$p.value
            p.val.Hombs <- wilcox.test(bac_i[NoObesos!="Female"] ~NoObesos[NoObesos!="Female"],exact=F)$p.value
            loc.Sex<-diff(c(median(bac_i_noob[NoObesos!="PCOS"]),median(bac_i_noob[NoObesos!="PCOS"])))
            loc.Fms <-diff(c(median(bac_i_noob[NoObesos!="Male"]),median(bac_i_noob[NoObesos!="Male"])))
            loc.Hombs <-diff(c(median(bac_i_noob[NoObesos!="Female"]),median(bac_i_noob[NoObesos!="Female"])))
            
          }
          
          
          Obesos <- grupo[obesidad=="Obese"]
          bac_i_ob <- bac_i[obesidad=="Obese"]
          
          
          if(s.prueba=="param" && s.interaccion < 0.05){
            
            mod <- lm(bac_i~Obesos)
            modhet <- leveneTest(mod)$`Pr(>F)`
            if(modhet<0.05){
              p.val.Sex_ob[j] <- t.test(bac_i_ob[Obesos!="PCOS"] ~Obesos[Obesos!="PCOS"],var.equal=F)$p.value
              p.val.Fms_ob[j] <- t.test(bac_i_ob[Obesos!="Male"] ~Obesos[Obesos!="Male"],var.equal=F)$p.value
              p.val.Hombs_ob[j]<- t.test(bac_i_ob[Obesos!="Female"] ~Obesos[Obesos!="Female"],var.equal=F)$p.value
              loc.Sex_ob[j]<-diff(c(mean(bac_i_ob[Obesos!="PCOS"]),mean(bac_i_ob[Obesos!="PCOS"])))
              loc.Fms_ob[j] <-diff(c(mean(bac_i_ob[Obesos!="Male"]),mean(bac_i_ob[Obesos!="Male"])))
              loc.Hombs_ob[j] <-diff(c(mean(bac_i_ob[Obesos!="Female"]),mean(bac_i_ob[Obesos!="Female"])))
              
              
            }else{
              p.val.Sex_ob[j] <- t.test(bac_i_ob[Obesos!="PCOS"] ~Obesos[Obesos!="PCOS"],var.equal=T)$p.value
              p.val.Fms_ob[j] <- t.test(bac_i_ob[Obesos!="Male"] ~Obesos[Obesos!="Male"],var.equal=T)$p.value
              p.val.Hombs_ob[j] <- t.test(bac_i_ob[Obesos!="Female"] ~Obesos[Obesos!="Female"],var.equal=T)$p.value
              loc.Sex_ob[j]<-diff(c(mean(bac_i_ob[Obesos!="PCOS"]),mean(bac_i_ob[Obesos!="PCOS"])))
              loc.Fms_ob[j] <-diff(c(mean(bac_i_ob[Obesos!="Male"]),mean(bac_i_ob[Obesos!="Male"])))
              loc.Hombs_ob[j] <-diff(c(mean(bac_i_ob[Obesos!="Female"]),mean(bac_ibac_i_ob)))
              
              
            }
            
            
          }else if(s.prueba=="noparam" && s.interaccion <0.05 ){
            p.val.Sex_ob[j] <- wilcox.test(bac_i_ob[Obesos!="PCOS"] ~Obesos[Obesos!="PCOS"],exact=F)$p.value
            p.val.Fms_ob[j] <- wilcox.test(bac_i_ob[Obesos!="Male"] ~Obesos[Obesos!="Male"],exact=F)$p.value
            p.val.Hombs[j] <- wilcox.test(bac_i[Obesos!="Female"] ~Obesos[Obesos!="Female"],exact=F)$p.value
            loc.Sex_ob[j]<-diff(c(median(bac_i_ob[Obesos!="PCOS"]),median(bac_i_ob[Obesos!="PCOS"])))
            loc.Fms_ob[j] <-diff(c(median(bac_i_ob[Obesos!="Male"]),median(bac_i_ob[Obesos!="Male"])))
            loc.Hombs_ob[j] <-diff(c(median(bac_i_ob[Obesos!="Female"]),median(bac_i_ob[Obesos!="Female"])))
            
          }
          if(s.interaccion<0.05 ){
            
            Obesos <- grupo[obesidad=="Obese"]
            bac_i_ob <- bac_i[obesidad=="Obese"]
            
            
            if(s.prueba=="param" && s.interaccion < 0.05){
              
              mod <- lm(bac_i~Obesos)
              modhet <- leveneTest(mod)$`Pr(>F)`
              if(modhet<0.05){
                p.val.Sex <- t.test(bac_i_ob[Obesos!="PCOS"] ~Obesos[Obesos!="PCOS"],var.equal=F)$p.value
                p.val.Fms <- t.test(bac_i_ob[Obesos!="Male"] ~Obesos[Obesos!="Male"],var.equal=F)$p.value
                p.val.Hombs <- t.test(bac_i_ob[Obesos!="Female"] ~Obesos[Obesos!="Female"],var.equal=F)$p.value
                loc.Sex<-diff(c(mean(bac_i_ob[Obesos!="PCOS"]),mean(bac_i_ob[Obesos!="PCOS"])))
                loc.Fms <-diff(c(mean(bac_i_ob[Obesos!="Male"]),mean(bac_i_ob[Obesos!="Male"])))
                loc.Hombs <-diff(c(mean(bac_i_ob[Obesos!="Female"]),mean(bac_i_ob[Obesos!="Female"])))
                
                
              }else{
                p.val.Sex <- t.test(bac_i_ob[Obesos!="PCOS"] ~Obesos[Obesos!="PCOS"],var.equal=T)$p.value
                p.val.Fms <- t.test(bac_i_ob[Obesos!="Male"] ~Obesos[Obesos!="Male"],var.equal=T)$p.value
                p.val.Hombs <- t.test(bac_i_ob[Obesos!="Female"] ~Obesos[Obesos!="Female"],var.equal=T)$p.value
                loc.Sex<-diff(c(mean(bac_i_ob[Obesos!="PCOS"]),mean(bac_i_ob[Obesos!="PCOS"])))
                loc.Fms <-diff(c(mean(bac_i_ob[Obesos!="Male"]),mean(bac_i_ob[Obesos!="Male"])))
                loc.Hombs <-diff(c(mean(bac_i_ob[Obesos!="Female"]),mean(bac_ibac_i_ob)))
                
                
              }
              
              
            }else if(s.prueba=="noparam" && s.interaccion <0.05 ){
              p.val.Sex <- wilcox.test(bac_i_ob[Obesos!="PCOS"] ~Obesos[Obesos!="PCOS"],exact=F)$p.value
              p.val.Fms <- wilcox.test(bac_i_ob[Obesos!="Male"] ~Obesos[Obesos!="Male"],exact=F)$p.value
              p.val.Hombs <- wilcox.test(bac_i[Obesos!="Female"] ~Obesos[Obesos!="Female"],exact=F)$p.value
              loc.Sex<-diff(c(median(bac_i_ob[Obesos!="PCOS"]),median(bac_i_ob[Obesos!="PCOS"])))
              loc.Fms <-diff(c(median(bac_i_ob[Obesos!="Male"]),median(bac_i_ob[Obesos!="Male"])))
              loc.Hombs <-diff(c(median(bac_i_ob[Obesos!="Female"]),median(bac_i_ob[Obesos!="Female"])))
              
            }
            
            
            
          }
          
        }
        ## obesos*grupo
        mod <- lm(bac_i ~ obesidad*grupo,contrasts = list(grupo=contr.sum,obesidad=contr.sum))
        p.norm <- shapiro.test(residuals(mod))$p.value
        

        
      }
    
  }





