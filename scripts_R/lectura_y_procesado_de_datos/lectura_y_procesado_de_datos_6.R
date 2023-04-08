









######
# SCRIPT PARA PREPARAR VARIABLES #

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
    "car",
    "zCompositions"
  )

### Preparamos librerias BIoconductor
list.of.bioc.packages <-
  c("phyloseq",
    "pcaMethods",
    "ropls",
    "mixOmics",
    "limma",
    "DESeq2",
    "ALDEx2")

### chequeamos que todo esté bien
res <- unlist(lapply(list.of.packages, require, character.only = T))
### chqueamos que todo este bien en bioconducot
res.bioc <-
  unlist(lapply(list.of.bioc.packages, require, character.only = T))

### Procesamos los datos de las bacterias

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
                                      datos.comunes[datos.comunes$interaccion == grupo[i], "Paciente"], ]), 2, function(x)
                                        sum(x) == 0))
  }
  
  w2 <- vector(mode = "list", length = 6)
  names(w2) <- levels(interaccion)
  for (i in 1:6) {
    w2[[names(w2)[i]]] <-
      which(apply(as.matrix(datos_i[rownames(datos_i) %in%
                                      datos.comunes[datos.comunes$interaccion == grupo[i], "Paciente"], ]), 2, function(x)
                                        sum(x) > 0))
  }
  
  w3 <- lapply(w2, names)
  common_variables <- Reduce(intersect, x = w3)
  
  return(common_variables)
}
pretratado_clinico <-
  function(general_data,
           variables_in_bacteria,
           comunes) {
    if (comunes) {
      sujetos <- variables_in_bacteria$Paciente
      Clin.tmp <-
        variables_in_bacteria[, -c(1, 2, 3, 4, 5)] # get rid off orden, paciente, sex, group and obese
      auc_clin <-
        colnames(Clin.tmp)[grep("AUC", colnames(Clin.tmp))]
      mean_basal_clin <-
        colnames(Clin.tmp)[grep("_B$", colnames(Clin.tmp))]
      sog <-
        colnames(Clin.tmp)[grep("^SO[GLP]", colnames(Clin.tmp))]
      variables_clin <- colnames(Clin.tmp)
      testosterona <-
        colnames(Clin.tmp)[grep("TEST", colnames(Clin.tmp))]
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
      Clin <- Clin.tmp[, -w]
      rownames(Clin) <- sujetos
      colnames(Clin)
      
    } else{
      sujetos <- general_data$Paciente
      Clin.tmp <-
        general_data[, -c(1, 2, 3, 4, 5)] # get rid off orden, paciente, sex, group and obese
      auc_clin <-
        colnames(Clin.tmp)[grep("AUC", colnames(Clin.tmp))]
      mean_basal_clin <-
        colnames(Clin.tmp)[grep("_B$", colnames(Clin.tmp))]
      sog <-
        colnames(Clin.tmp)[grep("^SO[GLP]", colnames(Clin.tmp))]
      variables_clin <- colnames(Clin.tmp)
      testosterona <-
        colnames(Clin.tmp)[grep("TEST", colnames(Clin.tmp))]
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
      Clin <- Clin.tmp[, -w]
      rownames(Clin) <- sujetos
      colnames(Clin)
      
    }
    return(Clin)
    
  }
pretratado_IP <-
  function(IP,
           variables_in_bacteria,
           common,
           general_data) {
    if (common) {
      sujetos <- variables_in_bacteria$Paciente
      IP.tmp <-
        IP[IP$Paciente %in% variables_in_bacteria$Paciente, ]
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
      
    } else{
      sujetos <- general_data$Paciente
      IP.tmp <- IP
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
      
    }
    return(IP.basal)
  }


pretratado_metaboloma <-
  function(metabolome,
           variables_in_bacteria,
           general_data,
           common) {
    ## si queremos las observaciones comunes y procesamos
    if (common) {
      sujetos <- variables_in_bacteria$Paciente
      grupo <- variables_in_bacteria$GROUP
      obesidad <- variables_in_bacteria$OBESE
      metabolome.basal.tmp <-
        metabolome[metabolome$Order %in% variables_in_bacteria$Orden,]
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
      
      metaboloma <- metabolome.basal_mean
      
    } else{
      sujetos <- general_data$Paciente
      grupo <- general_data$GROUP
      obesidad <- general_data$OBESE
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
      metabolome.basal_mean <- metabolome.basal.tmp[, -w]
      rownames(metabolome.basal_mean) <- sujetos
      
      metaboloma <- metabolome.basal_mean
      
      
    }
    ### hacemos doble normalizacion con voom y el metodo escogido
    if (common) {
      if (nrow(metaboloma) == nrow(variables_in_bacteria)) {
        metaboloma <- t(metaboloma)
      }
      met37 <- (metaboloma)
     
    } else{
      if (nrow(metaboloma) == nrow(general_data)) {
        metaboloma <- t(metaboloma)
      }
      metaboloma <- (metaboloma)
      met37 <- (metaboloma)
    
      
      
    }
    
    
    datos.metaboloma <- met37
    
    
    return(datos.metaboloma)
  }

post_procesado_aldex.t <- function(genera.abs,
                                   phylum.abs,
                                   grupo.,
                                   obesidad.,
                                   sexo.,
                                   condiciones) {
  # hacemos el procesado con aldex y todas las observaciones
  mc.samples <- 1000
  genus.abs <- genera.abs[, -c(
    grep("^ud-", colnames(genera.abs)),
    grep("GROUP", colnames(genera.abs)),
    grep("SEX", colnames(genera.abs)),
    grep("OBESE", colnames(genera.abs))
  )]
  
  phylum.abs <- phylum.abs[, -c(
    grep("^ud.", colnames(phylum.abs)),
    grep("GROUP", colnames(phylum.abs)),
    grep("SEX", colnames(phylum.abs)),
    grep("OBESE", colnames(phylum.abs))
  )]
  
  phylum.abs <- phylum.abs[, colSums(phylum.abs != 0) > 0]
  if (condiciones == "obesidad") {
    conds <- obesidad.
    
    phylum.aldex <-
      aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
    genus.aldex <-
      aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
    
    
  } else if (condiciones == "grupo") {
    conds <- grupo.
    
    phylum.aldex <-
      aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
    genus.aldex <-
      aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
    
  } else if (condiciones == "interaccion") {
    conds <- interaction(grupo., obesidad.)
    
    phylum.aldex <-
      aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
    genus.aldex <-
      aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
    
  } else if (condiciones == "sex") {
    conds <-sexo.
    
    phylum.aldex <-
      aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
    genus.aldex <-
      aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
    
  } else{
    phylum.aldex <-
      aldex.clr(t(phylum.abs), mc.samples = mc.samples)
    genus.aldex <-
      aldex.clr(t(genus.abs), mc.samples = mc.samples)
    
    
  }

  
  
  return(
    list(
      genero = genus.aldex,
      filo = phylum.aldex,
      genus.abs = genus.abs,
      phylum.abs = phylum.abs
    )
  )
  
}


post_procesado_aldex.t.filtrado <-
  function(genera.abs,
           phylum.abs,
           grupo.,
           obesidad.,
           sexo.,
           condiciones) {
    # aldex filtrado normal
    mc.samples <- 1000
    genus.abs <- genera.abs[, -c(
      grep("^ud-", colnames(genera.abs)),
      grep("GROUP", colnames(genera.abs)),
      grep("SEX", colnames(genera.abs)),
      grep("OBESE", colnames(genera.abs))
    )]
    
    
    phylum.abs <- phylum.abs[, -c(
      grep("^ud.", colnames(phylum.abs)),
      grep("GROUP", colnames(phylum.abs)),
      grep("SEX", colnames(phylum.abs)),
      grep("OBESE", colnames(phylum.abs))
    )]
    
    phylum.abs <- phylum.abs[, colSums(phylum.abs != 0) > 0]
    phylum.abs <-
      low.count.removal.counts(phylum.abs, 5)$data.filter
    
    genus.abs <-
      low.count.removal.counts(genus.abs, 120)$data.filter
    
    if (condiciones == "obesidad") {
      conds <- obesidad.
      
      phylum.aldex <-
        aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
      genus.aldex <-
        aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
      
      
    } else if (condiciones == "grupo") {
      conds <- grupo.
      
      phylum.aldex <-
        aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
      genus.aldex <-
        aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
      
    } else if (condiciones == "interaccion") {
      conds <- interaction(grupo., obesidad.)
      
      phylum.aldex <-
        aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
      genus.aldex <-
        aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
      
    } else if (condiciones == "sex") {
      conds <- sexo.
      
      phylum.aldex <-
        aldex.clr(t(phylum.abs), conds = conds, mc.samples = mc.samples)
      genus.aldex <-
        aldex.clr(t(genus.abs), conds = conds, mc.samples = mc.samples)
      
    } else{
      phylum.aldex <-
        aldex.clr(t(phylum.abs), mc.samples = mc.samples)
      genus.aldex <-
        aldex.clr(t(genus.abs), mc.samples = mc.samples)
      
    }

    return(
      list(
        genero = genus.aldex,
        filo = phylum.aldex,
        genus.abs = genus.abs,
        phylum.abs = phylum.abs
      )
    )
    
  }


post_procesado_bayesian <- function(genera.abs, phylum.abs) {
  # procesado bayesiano filtrado porque tiene que ser asi
  # clr e ilr
  genus.abs <- genera.abs[, -c(
    grep("^ud-", colnames(genera.abs)),
    grep("GROUP", colnames(genera.abs)),
    grep("SEX", colnames(genera.abs)),
    grep("OBESE", colnames(genera.abs))
  )]
  
  
  phylum.abs <- phylum.abs[, -c(
    grep("^ud.", colnames(phylum.abs)),
    grep("GROUP", colnames(phylum.abs)),
    grep("SEX", colnames(phylum.abs)),
    grep("OBESE", colnames(phylum.abs))
  )]
  
  phylum.abs <- phylum.abs[, colSums(phylum.abs != 0) > 0]
  phylum.abs <- low.count.removal.counts(phylum.abs, 5)$data.filter
  phylum.bayes <- cmultRepl((phylum.abs), output = "p-counts")
  phylum.bayes.clr <-
    matrix((mixOmics::logratio.transfo(phylum.bayes, "CLR")),
           nrow = nrow(phylum.abs),
           ncol = ncol(phylum.abs)
    )
  rownames(phylum.bayes.clr) <- rownames(phylum.abs)
  colnames(phylum.bayes.clr) <- colnames(phylum.abs)
  
  
  genus.abs <- low.count.removal.counts(genus.abs, 120)$data.filter
  
  genus.bayes <- cmultRepl((genus.abs), output = "p-counts")
  genus.bayes.clr <-
    matrix((mixOmics::logratio.transfo(genus.bayes, "CLR")),
           nrow = nrow(genus.abs),
           ncol = ncol(genus.abs)
    )
  rownames(genus.bayes.clr) <- rownames(genus.abs)
  colnames(genus.bayes.clr) <- colnames(genus.abs)
  
  
  return(list(genero.clr = genus.bayes.clr,
              filo.clr = phylum.bayes.clr))
}



read_integromics <- function(ruta, condiciones_) {
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

  general_data <-
    read.xlsx(file.path(ruta, "Integromics_1.xlsx"), sheetIndex = 1)
  general_data$SEX <-
    factor(general_data$SEX,
           levels = c(0, 2),
           labels = c("Female", "Male"))
  general_data$OBESE <-
    factor(
      general_data$OBESE,
      levels = c(0, 1),
      labels = c("No Obese", "Obese")
    )
  general_data$GROUP <-
    factor(
      general_data$GROUP,
      levels = c(0, 1, 2),
      labels = c("Female",  "PCOS", "Male")
    )
  IP_file <-
    file.path(file.path(ruta, "Integromics_IPmarkers.xlsx")) ## file path of the data
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
  
  metabolome_file <-
    file.path(file.path(ruta, "Integromics_Metabolome.xlsx")) ## file path of the data
  metabolome.tmp <-
    read.xlsx(metabolome_file, sheetIndex = 1) # read the data
  metabolites_name <-
    read.xlsx(metabolome_file, sheetIndex = 2)$Metabolitos # read the names of the metabolites
  any(is.na(metabolome.tmp)) # there is a row of missing values.
  metabolome <-
    metabolome.tmp[-nrow(metabolome.tmp), ] ## remove the missing value row
  metabolome$GROUP <- general_data$GROUP ## add GROUP variable
  metabolome$OBESE <- general_data$OBESE ## add OBESE variables.
  
  bacteria_phylum.abs <-
    as.data.frame(read.xlsx(
      file.path(ruta, "integromics_microbiota.xlsx"),
      sheetIndex = 1
    ))
  
  bacteria_genera.abs <-
    as.data.frame(t(read.xlsx(
      file.path(ruta, "integromics_microbiota.xlsx"),
      sheetIndex = 3
    )))
  phylum.list <-
    proces_bacteria(bacteria_phylum.abs, tipo = "phylum")
  genera.list <-
    proces_bacteria(bacteria_genera.abs, tipo = "genera")
  
  phylum.abs <- phylum.list$abs
  phylum.rel <- phylum.list$rel
  genera.abs <- genera.list$abs
  genera.rel <- genera.list$rel
  
  
  variables_in_bacteria <-
    general_data[general_data$Orden %in% bacteria_phylum.abs$Order,]
  sexo.d <- variables_in_bacteria$SEX
  obesidad.d <- variables_in_bacteria$OBESE
  grupo.d <- variables_in_bacteria$GROUP
  
  
  
  datos.clinicos.comunes <-
    pretratado_clinico(general_data, variables_in_bacteria, comunes = T)
  datos.clinicos.totales <-
    pretratado_clinico(general_data, variables_in_bacteria, comunes = F)
  datos.IP.comunes <-
    pretratado_IP(IP, variables_in_bacteria, common = T, general_data)
  datos.IP.totales <-
    pretratado_IP(IP, variables_in_bacteria , common = F, general_data)
  ## microbioma
  datos.microbiota <-
    post_procesado_aldex.t(
      genera.abs,
      phylum.abs,
      grupo. = variables_in_bacteria$GROUP,
      obesidad. = variables_in_bacteria$OBESE,
      sexo. = variables_in_bacteria$SEX,
      condiciones = condiciones_
    )
  datos.microbiota.filt <-
    post_procesado_aldex.t.filtrado(
      genera.abs,
      phylum.abs,
      grupo. = variables_in_bacteria$GROUP,
      obesidad. = variables_in_bacteria$OBESE,
      sexo. = variables_in_bacteria$SEX,
      condiciones = condiciones_
    )
  datos.microbiota.bayesianos <-
    post_procesado_bayesian(genera.abs, phylum.abs)
  
  genero.clr <- datos.microbiota.bayesianos$genero.clr
  filo.clr <- datos.microbiota.bayesianos$filo.clr
  
  
  genero.aldex <- datos.microbiota$genero
  filo.aldex <- datos.microbiota$filo
  genus.abs <- datos.microbiota$genus.abs
  phylum.abs <- datos.microbiota$phylum.abs
  
  genus.abs.filt <- datos.microbiota.filt$genus.abs
  phylum.abs.filt <- datos.microbiota.filt$phylum.abs
  genero.aldex.filt <- datos.microbiota.filt$genero
  filo.aldex.filt <- datos.microbiota.filt$filo
  
  microbiota <- list(
    bayes.clr = list(filo.clr = filo.clr,
                     genero.clr = genero.clr),
    aldex = list(filo.aldex = filo.aldex,
                 genero.aldex = genero.aldex),
    aldex.filt = list(filo.aldex = filo.aldex.filt,
                      genero.aldex = genero.aldex.filt),
    microbiota.abs = list(genero = genus.abs,
                          filo = phylum.abs)
    
  )
  
  metaboloma.comun <- vector(mode = "list", length = 5)
  metaboloma.total <- vector(mode = "list", length = 5)
  

    metaboloma.comun <- pretratado_metaboloma(
      metabolome = metabolome,
      variables_in_bacteria = variables_in_bacteria,
      general_data = general_data,
      common = T)
    metaboloma.total<- pretratado_metaboloma(
      metabolome = metabolome,
      variables_in_bacteria = variables_in_bacteria,
      general_data = general_data,
      common = F)
  

  metaboloma <- list(metaboloma.comun = metaboloma.comun,
                     metaboloma.total = metaboloma.total)
  
  
  totales <- list(
    clinicos = datos.clinicos.totales,
    ip = datos.IP.totales,
    metaboloma = metaboloma.total,
    grupo = general_data$GROUP,
    obesidad = general_data$OBESE,
    general_data = general_data
  )
  comunes <- list(
    clinicos = datos.clinicos.comunes,
    ip = datos.IP.comunes,
    metaboloma = metaboloma.comun,
    microbiota = microbiota,
    grupo = variables_in_bacteria$GROUP,
    obesidad = variables_in_bacteria$OBESE,
    variables_in_bacteria = variables_in_bacteria
  )
  
  
  
  retorno <- list(totales = totales, comunes = comunes)
  return(retorno)
}


conds <- c("obesidad", "grupo", "interaccion", "sex", "none")

for (cond in conds) {
  print(paste(cond))
  datos.integromics <-
    read_integromics("../../datos", condiciones_  = cond)
  
  if (!dir.exists("../../datos/resultados_analisis_datos_procesados_1000_novoom")) {
    dir.create("../../datos/resultados_analisis_datos_procesados_1000_novoom")
  }
  
  saveRDS(datos.integromics,
          file.path(
            paste0(
              "../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom",
              "_",
              cond,
              "_",
              ".rds"
            )
          ))
  
  
  
}
