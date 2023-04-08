library(car)
library(dplyr)
library(lmPerm)
fn.paired <- function(x, g, termino1, termino2, name_g) {
  x1 <- x[g == termino1]
  x2 <- x[g == termino2]
  mod <- lm(x ~ g)
  if (length(unique(residuals(mod))) == 1) {
    prueba.normalidad <- NA
    prueba.homocedasticidad <- NA
  } else{
    prueba.normalidad <- shapiro.test(residuals(mod))$p.value
    prueba.homocedasticidad <- leveneTest(mod)$`Pr(>F)`[1]
  }
  
  if (length(unique(x1)) == 1 ||
      length(unique(x2)) == 1 || is.na(prueba.normalidad)) {
    if (length(unique(x1)) == 1) {
      p.prueba <- NA
      p.location <- NA
      tipo.prueba <- paste0("Ausencia de ", termino1, " en ", name_g)
      
    } else if (length(unique(x2)) == 1) {
      p.prueba <- NA
      p.location <- NA
      tipo.prueba <- paste0("Ausencia de ", termino2, " en ", name_g)
      
    } else if (length(unique(x1)) == 1 && length(unique(x2)) == 1) {
      p.prueba <- NA
      p.location <- NA
      tipo.prueba <-
        paste0("Ausencia de ", termino2, " y ", termino1, " en ", name_g)
      
    }
  }
  else if (prueba.normalidad < 0.05) {
    p.prueba <- wilcox.test(x ~ g, exact = F)$p.value
    p.location <- diff(c(median(x1), median(x2)))
    tipo.prueba <- 0
    
  } else if (prueba.normalidad > 0.05 &&
             prueba.homocedasticidad > 0.05) {
    p.prueba <- t.test(x ~ g, var.equal = T)$p.value
    p.location <- diff(c(mean(x1), mean(x2)))
    tipo.prueba <- 1
    
  } else if (prueba.normalidad > 0.05 &&
             prueba.homocedasticidad < 0.05) {
    p.prueba <- t.test(x ~ g, var.equal = F)$p.value
    p.location <- diff(c(mean(x1), mean(x2)))
    tipo.prueba <- 1
  }
  
  return(bind_cols(
    p.prueba = p.prueba,
    p.location = p.location,
    tipo.prueba = tipo.prueba
  ))
}


inferencia_univariante <- function(X,grupo,obesidad,progress=F){
  
  
  p <- ncol(X)
  n <- nrow(X)
  

  ### VARIABLES
  obesidad.mod <- vector("numeric", length = p)
  s.obesidad.mod1 <- vector("numeric", length = p)
  s.grupo.mod1 <- vector("numeric", length = p)
  s.interaccion.mod1 <- vector("numeric", length = p)
  s.prueba.mod1 <-  vector("numeric", length = p)
  
  
  
  obesidad.pvalue <-  vector("numeric", length = p)
  diferencia_loc_obesidad <- vector("numeric", length = p)
  obesidad.prueba <- vector("numeric", length = p)
  
  
  
  n.sexo <- length(grupo[grupo != "PCOS"])
  sexo.pvalue <-  vector("numeric", length = p)
  diferencia_loc_sexo <- vector("numeric", length = p)
  sexo.prueba <- vector("numeric", length = p)
  
  n.malespcos <- length(grupo[grupo != "Female"])
  malespcos.pvalue <-  vector("numeric", length = p)
  diferencia_loc_malespcos <- vector("numeric", length = p)
  malespcos.prueba <- vector("numeric", length = p)
  
  n.femalespcos <- length(grupo[grupo != "Male"])
  femalespcos.pvalue <-  vector("numeric", length = p)
  diferencia_loc_femalespcos <-
    vector("numeric", length = p)
  femalespcos.prueba <- vector("numeric", length = p)
  
  
  
  
  sexo.pvalue.no.obese <-  vector("numeric", length = p)
  diferencia_loc_sexo.no.obese <-
    vector("numeric", length = p)
  sexo.prueba.no.obese <-
    vector("numeric", length = p)
  
  
  malespcos.pvalue.no.obese <-
    vector("numeric", length = p)
  diferencia_loc_malespcos.no.obese <-
    vector("numeric", length = p)
  malespcos.prueba.no.obese <-
    vector("numeric", length = p)
  
  
  femalespcos.pvalue.no.obese <-
    vector("numeric", length = p)
  diferencia_loc_femalespcos.no.obese <-
    vector("numeric", length = p)
  femalespcos.prueba.no.obese <-
    vector("numeric", length = p)
  
  
  
  sexo.pvalue.obese <-  vector("numeric", length = p)
  diferencia_loc_sexo.obese <-
    vector("numeric", length = p)
  sexo.prueba.obese <- vector("numeric", length = p)
  
  
  malespcos.pvalue.obese <-
    vector("numeric", length = p)
  diferencia_loc_malespcos.obese <-
    vector("numeric", length = p)
  malespcos.prueba.obese <-
    vector("numeric", length = p)
  
  
  femalespcos.pvalue.obese <-
    vector("numeric", length = p)
  diferencia_loc_femalespcos.obese <-
    vector("numeric", length = p)
  femalespcos.prueba.obese <-
    vector("numeric", length = p)
  
  
  
  
  sexo.pvalue.no.obese2 <-
    vector("numeric", length = p)
  diferencia_loc_sexo.no.obese2 <-
    vector("numeric", length = p)
  sexo.prueba.no.obese2 <-
    vector("numeric", length = p)
  
  malespcos.pvalue.no.obese2 <-
    vector("numeric", length = p)
  diferencia_loc_malespcos.no.obese2 <-
    vector("numeric", length = p)
  malespcos.prueba.no.obese2 <-
    vector("numeric", length = p)
  
  
  femalespcos.pvalue.no.obese2 <-
    vector("numeric", length = p)
  diferencia_loc_femalespcos.no.obese2 <-
    vector("numeric", length = p)
  femalespcos.prueba.no.obese2 <-
    vector("numeric", length = p)
  
  
  sexo.pvalue.obese2 <-  vector("numeric", length = p)
  diferencia_loc_sexo.obese2 <-
    vector("numeric", length = p)
  sexo.prueba.obese2 <- vector("numeric", length = p)
  
  
  malespcos.pvalue.obese2 <-
    vector("numeric", length = p)
  diferencia_loc_malespcos.obese2 <-
    vector("numeric", length = p)
  malespcos.prueba.obese2 <-
    vector("numeric", length = p)
  
  femalespcos.pvalue.obese2 <-
    vector("numeric", length = p)
  diferencia_loc_femalespcos.obese2 <-
    vector("numeric", length = p)
  femalespcos.prueba.obese2 <-
    vector("numeric", length = p)
  
  
  
  pcos.p.value.obesidad <-
    vector("numeric", length = p)
  pcos.diferencia.obesidad <-
    vector("numeric", length = p)
  pcos.prueba.obesidad <-
    vector("numeric", length = p)
  
  fems.p.value.obesidad <- vector("numeric", length = p)
  fems.diferencia.obesidad <-
    vector("numeric", length = p)
  fems.prueba.obesidad <- vector("numeric", length = p)
  
  males.p.value.obesidad <-
    vector("numeric", length = p)
  males.diferencia.obesidad <-
    vector("numeric", length = p)
  males.prueba.obesidad <-
    vector("numeric", length = p)
  

for (j in 1:p) {
  if(progress){
    print(paste("variable", colnames(X)[j], "numero ", j))
    
  }
  var_i <- X[, j]
  modelo1 <-
    lm(var_i ~ grupo * obesidad,
       contrasts = list(grupo = contr.sum, obesidad = contr.sum))
  ### vemos las suposiciones de los modelos
  
  prueba.normalidad <- shapiro.test(residuals(modelo1))$p.value
  prueba.homocedasticidad <- leveneTest(modelo1)$`Pr(>F)`[1]
  
  ### si es cumple los criterios de normalidad
  
  if (prueba.normalidad < 0.05) {
    s <-
      summary(aovp(
        var_i ~ grupo * obesidad,
        contrasts = list(grupo = contr.sum, obesidad = contr.sum)
        ,
        settings = F
      ))
    s.obesidad.mod1[j] <- s[[1]]$`Pr(Prob)`[2]
    s.grupo.mod1[j] <- s[[1]]$`Pr(Prob)`[1]
    s.interaccion.mod1[j] <- s[[1]]$`Pr(Prob)`[3]
    s.prueba.mod1[j] <- 0
    
    
  } else{
    s <- Anova(modelo1, type = 3)
    s.obesidad.mod1[j] <- s$`Pr(>F)`[2]
    s.grupo.mod1[j] <- s$`Pr(>F)`[3]
    s.interaccion.mod1[j] <- s$`Pr(>F)`[4]
    s.prueba.mod1[j] <- 1
    
  }
  
  
  ### veamos la obesidad del primer modelo
  
  
  vector_resultados <-
    fn.paired(
      x = var_i,
      g = obesidad,
      termino1 = "No Obese",
      termino2 = "Obese",
      name_g = "obesidad"
    )
  obesidad.pvalue[j] = vector_resultados[1]
  diferencia_loc_obesidad[j] = vector_resultados[2]
  obesidad.prueba[j] = vector_resultados[3]
  
  
  
  ### grupo modelo 1
  
  ### sexo
  var_i_sexo <- var_i[grupo != "PCOS"]
  grupo_sexo <- grupo[grupo != "PCOS"]
  vector_resultados <-
    fn.paired(
      x = var_i_sexo,
      g = grupo_sexo,
      termino1 = "Female",
      termino2 = "Male",
      name_g = "sexo"
    )
  sexo.pvalue[j] = vector_resultados[1]
  diferencia_loc_sexo[j] = vector_resultados[2]
  sexo.prueba[j] = vector_resultados[3]
  ### hombres vs mujeres PCOS
  
  var_i_malespcos <- var_i[grupo != "Female"]
  grupo_malespcos <- grupo[grupo != "Female"]
  vector_resultados <-
    fn.paired(
      x = var_i_malespcos,
      g = grupo_malespcos,
      termino1 = "Male",
      termino2 = "PCOS",
      name_g = "pcos-hombres"
    )
  malespcos.pvalue[j] = vector_resultados[1]
  diferencia_loc_malespcos[j] = vector_resultados[2]
  malespcos.prueba[j] = vector_resultados[3]
  
  ### mujeres vs mujeres PCOS
  var_i_femalespcos <- var_i[grupo != "Male"]
  grupo_femalespcos <- grupo[grupo != "Male"]
  vector_resultados <-
    fn.paired(
      x = var_i_femalespcos,
      g = grupo_femalespcos,
      termino1 = "Female",
      termino2 = "PCOS",
      name_g = "pcos-mujeres"
    )
  femalespcos.pvalue[j] = vector_resultados[1]
  diferencia_loc_femalespcos[j] = vector_resultados[2]
  femalespcos.prueba[j] = vector_resultados[3]
  
  
  
  ### obesidad|grupo
  ### Vemos el sexo en las personas no obesas
  
  No.Obese_var_i_sexo <-
    var_i[obesidad == "No Obese" & grupo != "PCOS"]
  No.Obese_grupo_sexo <-
    grupo[obesidad == "No Obese" & grupo != "PCOS"]
  vector_resultados <-
    fn.paired(
      x = No.Obese_var_i_sexo,
      g = No.Obese_grupo_sexo,
      termino1 = "Female",
      termino2 = "Male",
      name_g = "sexo No Obesos"
    )
  sexo.pvalue.no.obese[j] = vector_resultados[1]
  diferencia_loc_sexo.no.obese[j] = vector_resultados[2]
  sexo.prueba.no.obese[j] = vector_resultados[3]
  ###### Vemos hombres vs mujeres PCOS no obesas
  
  No.Obese_var_i_malespcos <-
    var_i[obesidad == "No Obese" & grupo != "Female"]
  No.Obese_grupo_malespcos <-
    grupo[obesidad == "No Obese" & grupo != "Female"]
  vector_resultados <-
    fn.paired(
      x = No.Obese_var_i_malespcos,
      g = No.Obese_grupo_malespcos,
      termino1 = "Male",
      termino2 = "PCOS",
      name_g = "pcos-hombres no obesos"
    )
  malespcos.pvalue.no.obese[j] = vector_resultados[1]
  diferencia_loc_malespcos.no.obese[j] = vector_resultados[2]
  malespcos.prueba.no.obese[j] = vector_resultados[3]
  
  ### No obesidad en mujeres vs mujeres PCOS
  No.Obese_var_i_femalespcos <-
    var_i[obesidad == "No Obese" & grupo != "Male"]
  No.Obese_grupo_femalespcos <-
    grupo[obesidad == "No Obese" & grupo != "Male"]
  vector_resultados <-
    fn.paired(
      x = No.Obese_var_i_femalespcos,
      g = No.Obese_grupo_femalespcos,
      termino1 = "Female",
      termino2 = "PCOS",
      name_g = "pcos-mujeres no obesos"
    )
  femalespcos.pvalue.no.obese[j] = vector_resultados[1]
  diferencia_loc_femalespcos.no.obese[j] = vector_resultados[2]
  femalespcos.prueba.no.obese[j] = vector_resultados[3]
  
  
  ### obesidad|grupo
  ### Vemos el sexo en las personas obesas
  
  Obese_var_i_sexo <- var_i[obesidad == "Obese" & grupo != "PCOS"]
  Obese_grupo_sexo <-
    grupo[obesidad == "Obese" & grupo != "PCOS"]
  vector_resultados <-
    fn.paired(
      x = Obese_var_i_sexo,
      g = Obese_grupo_sexo,
      termino1 = "Female",
      termino2 = "Male",
      name_g = "sexo obsos"
    )
  sexo.pvalue.obese[j] = vector_resultados[1]
  diferencia_loc_sexo.obese[j] = vector_resultados[2]
  sexo.prueba.obese[j] = vector_resultados[3]
  ###### Vemos hombres vs mujeres PCOS obesas
  
  Obese_var_i_malespcos <-
    var_i[obesidad == "Obese" & grupo != "Female"]
  Obese_grupo_malespcos <-
    grupo[obesidad == "Obese" & grupo != "Female"]
  vector_resultados <-
    fn.paired(
      x = Obese_var_i_malespcos,
      g = Obese_grupo_malespcos,
      termino1 = "Male",
      termino2 = "PCOS",
      name_g = "pcos male obesos"
    )
  malespcos.pvalue.obese[j] = vector_resultados[1]
  diferencia_loc_malespcos.obese[j] = vector_resultados[2]
  malespcos.prueba.obese[j] = vector_resultados[3]
  
  ### obesidad en mujeres vs mujeres PCOS
  Obese_var_i_femalespcos <-
    var_i[obesidad == "Obese" & grupo != "Male"]
  Obese_grupo_femalespcos <-
    grupo[obesidad == "Obese" & grupo != "Male"]
  vector_resultados <-
    fn.paired(
      x = Obese_var_i_femalespcos,
      g = Obese_grupo_femalespcos,
      termino1 = "Female",
      termino2 = "PCOS",
      name_g = "pcos mujeres obesos"
    )
  femalespcos.pvalue.obese[j] = vector_resultados[1]
  diferencia_loc_femalespcos.obese[j] = vector_resultados[2]
  femalespcos.prueba.obese[j] = vector_resultados[3]
  
  
  ### AHORA VEREMOS EL GRUPO|OBESIDAD
  
  pcos_var_i <- var_i[grupo == "PCOS"]
  pcos.obesidad <- obesidad[grupo == "PCOS"]
  
  fems_var_i <- var_i[grupo == "Female"]
  fems.obesidad <- obesidad[grupo == "Female"]
  
  males_var_i <- var_i[grupo == "Male"]
  males.obesidad <- obesidad[grupo == "Male"]
  
  vector_resultados <-
    fn.paired(
      x = pcos_var_i,
      g = pcos.obesidad,
      termino1 = "No Obese",
      termino2 = "Obese",
      name_g = "pcos obesidad"
    )
  pcos.p.value.obesidad[j] <- vector_resultados[1]
  pcos.diferencia.obesidad[j]  <- vector_resultados[2]
  pcos.prueba.obesidad[j]  <- vector_resultados[3]
  
  vector_resultados <-
    fn.paired(
      x = fems_var_i,
      g = fems.obesidad,
      termino1 = "No Obese",
      termino2 = "Obese",
      name_g = "females obesidad"
    )
  fems.p.value.obesidad[j]  <- vector_resultados[1]
  fems.diferencia.obesidad[j]  <- vector_resultados[2]
  fems.prueba.obesidad[j]  <- vector_resultados[3]
  
  
  vector_resultados <-
    fn.paired(
      x = males_var_i,
      g = males.obesidad,
      termino1 = "No Obese",
      termino2 = "Obese",
      name_g = "hombres obesidad"
    )
  males.p.value.obesidad[j]  <- vector_resultados[1]
  males.diferencia.obesidad[j]  <- vector_resultados[2]
  males.prueba.obesidad[j]  <- vector_resultados[3]
  
}
  
  retorno_pvalores <- data.frame(
    cbind(
      obesidad = s.obesidad.mod1,
      p.value_OBESIDAD = obesidad.pvalue,
      grupo = s.grupo.mod1,
      interaccion = s.interaccion.mod1,
      p_value_Females_PCOS = femalespcos.pvalue,
      p_value_SEXO = sexo.pvalue,
      p_value_Males_PCOS = malespcos.pvalue,
      interaccion = s.interaccion.mod1,
      p_value_Females_PCOS_No_Obese = femalespcos.pvalue.no.obese,
      p_value_SEXO_No_Obese = sexo.pvalue.no.obese,
      p_value_Males_PCOS_No_Obese = malespcos.pvalue.no.obese,
      p_value_Females_PCOS_Obese = femalespcos.pvalue.obese,
      p_value_SEXO_Obese = sexo.pvalue.obese,
      p_value_Males_PCOS_Obese = malespcos.pvalue.obese,
      pcos_obesidad = pcos.p.value.obesidad,
      males_obesidad = males.p.value.obesidad,
      females_obesidad = fems.p.value.obesidad
      
    )
  )
  rownames(retorno_pvalores) <- colnames(X)
  retorno_tipo_prueba <- data.frame(
    cbind(
      tipo_prueba = s.prueba.mod1,
      tipo.prueba_OBESIDAD = obesidad.prueba,
      sexo_prueba = sexo.prueba,
      femalespcos_prueba = femalespcos.prueba,
      malespcos_prueba = malespcos.prueba,
      sexo_prueba_No_Obese = sexo.prueba.no.obese,
      femalespcos_prueba_No_Obese = femalespcos.prueba.no.obese,
      malespcos_prueba_No_Obese = malespcos.prueba.no.obese,
      sexo_prueba_Obese = sexo.prueba.obese,
      femalespcos_prueba_Obese = femalespcos.prueba.obese,
      malespcos_prueba_Obese = malespcos.prueba.obese,
      pcos_prueba_obesidad = pcos.prueba.obesidad,
      females_prueba_obesidad = fems.prueba.obesidad,
      males_prueba_obesidad = males.prueba.obesidad
    )
  )
  
  
  retorno_tipo_prueba[retorno_tipo_prueba == 0] <- "no parametrica"
  retorno_tipo_prueba[retorno_tipo_prueba == 1] <- "parametrica"
  rownames(retorno_tipo_prueba) <- colnames(X)
  
  
  
  
  diferencias <- data.frame(
    cbind(
      diferencia_OBESIDAD = diferencia_loc_obesidad,
      mu_Female_PCOS = diferencia_loc_femalespcos,
      mu_SEXO = diferencia_loc_sexo,
      mu_Males_PCOS = diferencia_loc_malespcos,
      mu_Female_PCOS_No_Obese = diferencia_loc_femalespcos.no.obese,
      mu_SEXO_No_Obese = diferencia_loc_sexo.no.obese,
      mu_Males_PCOS_No_Obese = diferencia_loc_malespcos.no.obese,
      mu_Female_PCOS_Obese = diferencia_loc_femalespcos.obese,
      mu_SEXO_Obese = diferencia_loc_sexo.obese,
      mu_Males_PCOS_Obese = diferencia_loc_malespcos.obese,
      pcos_diferencia_obesidad = pcos.diferencia.obesidad,
      females_diferencia_obesidad = fems.diferencia.obesidad,
      males_diferencia_obesidad = males.diferencia.obesidad
    )
  )
  
  rownames(diferencias) <- colnames(X)
  glm_modelo_1 <-
    bind_cols(diferencias, retorno_pvalores, retorno_tipo_prueba)
  
  rownames(glm_modelo_1) <- colnames(X)
  
  
  #### INTERPRETACION FINAL
  
  glm_modelo <- glm_modelo_1
  
  
  retorno <- list(p.values=retorno_pvalores,
                  diferencias=diferencias,
                  tipo_prueba=retorno_tipo_prueba,
                  glm_modelo=glm_modelo,
                  X=X)
  
}



funcion_interpretacion <- function(glm_modelo,X){
  
  p <- ncol(X)
  n <- nrow(X)
  
  interpretacion_final <- vector("character", length = p)
  interpretacion_obesidad <- vector("character", length = p)
  interpretacion_obesidad_v <- vector("character", length = p)
  interpretacion_grupo <- vector("character", length = p)
  interpretacion_fems_pcos <- vector("character", length = p)
  interpretacion_males_pcos <- vector("character", length = p)
  interpretacion_sex <- vector("character", length = p)
  interpretacion_interaccion <- vector("character", length = p)
  interpretacion_fems_pcos_no_obese <- vector("character", length = p)
  interpretacion_males_pcos_no_obese <-
    vector("character", length = p)
  interpretacion_sex_no_obese <- vector("character", length = p)
  interpretacion_fems_pcos_obese <- vector("character", length = p)
  interpretacion_males_pcos_obese <- vector("character", length = p)
  interpretacion_sex_obese <- vector("character", length = p)
  interpretacion_pcos_obesidad <- vector("character", length = p)
  interpretacion_fems_obesidad <- vector("character", length = p)
  interpretacion_males_obesidad <- vector("character", length = p)
  interpretacion_unica_grupo <- vector("character", length = p)
  interpretacion_unica_grupo_obesidad <-
    vector("character", length = p)
  interpretacion_unica_obesidad_grupo <-
    vector("character", length = p)
  
  
  for (j in 1:p) {
    name <- colnames(X)[j]
    inter_name <- paste("La variable", name, "es", j)
    print(inter_name)
    inter_i <- glm_modelo[j, ]
    if (inter_i$obesidad < 0.05 & !is.na(inter_i$obesidad)) {
      ## NO OBESOS - OBESOS
      if (as.numeric(inter_i$diferencia_OBESIDAD) < 0 &
          !is.na(inter_i$diferencia_OBESIDAD)) {
        interpretacion_obesidad[j] <-
          paste(inter_name,
                "MAYOR en OBESOS",
                inter_i$tipo.prueba_OBESIDAD)
      } else{
        interpretacion_obesidad[j] <-
          paste(inter_name,
                "MENOR en OBESOS",
                inter_i$tipo.prueba_OBESIDAD)
        
      }
      
      
    }
    if (inter_i$grupo < 0.05 & !is.na(inter_i$grupo)) {
      interpretacion_grupo[j] <-
        paste(name, "HAY DIFERENCIAS ENTRE LA VARIABLE GRUPO")
      if (inter_i$p_value_Females_PCOS < 0.05 &
          !is.na(inter_i$p_value_Females_PCOS)) {
        if (inter_i$mu_Female_PCOS < 0) {
          interpretacion_fems_pcos[j] <-
            paste(inter_name,
                  "MAYOR EN PCOS vs MUJERES",
                  inter_i$femalespcos_prueba)
        } else{
          interpretacion_fems_pcos[j] <-
            paste(inter_name,
                  "MAYOR EN MUJERES vs PCOS",
                  inter_i$femalespcos_prueba)
          
        }
      }
      if (inter_i$p_value_SEXO < 0.05 &
          !is.na(inter_i$p_value_SEXO)) {
        if (inter_i$mu_SEXO < 0) {
          interpretacion_sex[j] <-
            paste(inter_name,
                  "MAYOR EN HOMBRES vs MUJERES",
                  inter_i$sexo_prueba)
        } else{
          interpretacion_sex[j] <-
            paste(inter_name,
                  "MAYOR EN MUEJRES vs HOMBRES",
                  inter_i$sexo_prueba)
          
        }
      }
      if (inter_i$p_value_Males_PCOS < 0.05 &
          !is.na(inter_i$p_value_Males_PCOS)) {
        if (inter_i$mu_Males_PCOS < 0) {
          interpretacion_males_pcos[j] <-
            paste(inter_name,
                  "MAYOR EN PCOS vs HOMBRES",
                  inter_i$malespcos_prueba)
        } else{
          interpretacion_males_pcos[j] <-
            paste(inter_name,
                  "MAYOR EN HOMBRES vs PCOS",
                  inter_i$malespcos_prueba)
          
        }
      }
      if ((inter_i$p_value_Males_PCOS < 0.05 &
           !is.na(inter_i$p_value_Males_PCOS)) &
          (inter_i$p_value_Females_PCOS > 0.05 |
           is.na(inter_i$p_value_Females_PCOS)) &
          (inter_i$p_value_SEXO > 0.05 |
           is.na(inter_i$p_value_SEXO))) {
        interpretacion_unica_grupo[j] <-
          paste(inter_name, "SOLO ENTRE HOMBRES Y PCOS")
      } else if ((inter_i$p_value_Males_PCOS > 0.05 |
                  is.na(inter_i$p_value_Males_PCOS)) &
                 (inter_i$p_value_Females_PCOS < 0.05 &
                  !is.na(inter_i$p_value_Females_PCOS)) &
                 (inter_i$p_value_SEXO > 0.05 |
                  is.na(inter_i$p_value_SEXO))) {
        interpretacion_unica_grupo[j] <-
          paste(inter_name, "SOLO ENTRE MUJERES")
        
      } else if ((inter_i$p_value_Males_PCOS > 0.05 |
                  is.na(inter_i$p_value_Males_PCOS)) &
                 (inter_i$p_value_Females_PCOS > 0.05 |
                  is.na(inter_i$p_value_Females_PCOS)) &
                 (inter_i$p_value_SEXO < 0.05 &
                  !is.na(inter_i$p_value_SEXO))) {
        interpretacion_unica_grupo[j] <- paste(inter_name, "SOLO SEXO")
        
      }
      
    }
    if (inter_i$interaccion < 0.05 & !is.na(inter_i$interaccion)) {
      interpretacion_interaccion[j] <-
        paste(inter_name, "EXISTE INTERACCION")
      if (inter_i$p_value_Females_PCOS_No_Obese < 0.05 &
          !is.na(inter_i$p_value_Females_PCOS_No_Obese)) {
        if (inter_i$mu_Female_PCOS_No_Obese < 0.05) {
          interpretacion_fems_pcos_no_obese[j] <-
            paste(inter_name, "MAYOR EN PCOS NO OBESAS vs MUJERES SANAS")
        } else{
          interpretacion_fems_pcos_no_obese[j] <-
            paste(inter_name, "MENOR EN PCOS NO OBESAS vs MUJERES SANAS")
          
        }
      }
      if (inter_i$p_value_Males_PCOS_No_Obese < 0.05 &
          !is.na(inter_i$p_value_Males_PCOS_No_Obese)) {
        if (inter_i$mu_Males_PCOS_No_Obese < 0.05) {
          interpretacion_males_pcos_no_obese[j] <-
            paste(
              inter_name,
              "MAYOR EN PCOS NO OBESAS vs HOMBRES SANOS",
              inter_i$malespcos_prueba_No_Obese
            )
        } else{
          interpretacion_males_pcos_no_obese[j] <-
            paste(
              inter_name,
              "MENOR EN PCOS NO OBESAS vs HOMBRES SANOS",
              inter_i$malespcos_prueba_No_Obese
            )
          
        }
      }
      if (inter_i$p_value_SEXO_No_Obese < 0.05 &
          !is.na(inter_i$p_value_SEXO_No_Obese)) {
        if (inter_i$mu_SEXO_No_Obese < 0.05) {
          interpretacion_sex_no_obese[j] <-
            paste(
              inter_name,
              "MAYOR EN HOMBRES NO OBESOS vs MUJERES SANAS",
              inter_i$sexo_prueba_No_Obese
            )
        } else{
          interpretacion_sex_no_obese[j] <-
            paste(
              inter_name,
              "MENOR EN HOMBRES NO OBESOS vs MUJERES SANAS",
              inter_i$sexo_prueba_No_Obese
            )
          
        }
      }
      if (inter_i$p_value_Females_PCOS_Obese < 0.05  &
          !is.na(inter_i$p_value_Females_PCOS_Obese)) {
        if (inter_i$mu_Female_PCOS_Obese < 0.05) {
          interpretacion_fems_pcos_obese[j] <-
            paste(
              inter_name,
              "MAYOR EN PCOS OBESAS vs MUJERES OBESAS",
              inter_i$femalespcos_prueba_Obese
            )
        } else{
          interpretacion_fems_pcos_obese[j] <-
            paste(
              inter_name,
              "MENOR EN PCOS  OBESAS vs MUJERES OBESAS",
              inter_i$femalespcos_prueba_Obese
            )
          
        }
      }
      if (inter_i$p_value_Males_PCOS_Obese < 0.05  &
          !is.na(inter_i$p_value_Males_PCOS_Obese)) {
        if (inter_i$mu_Males_PCOS_Obese < 0.05) {
          interpretacion_males_pcos_obese[j] <-
            paste(
              inter_name,
              "MAYOR EN PCOS  OBESAS vs HOMBRES OBESOS",
              inter_i$malespcos_prueba_Obese
            )
        } else{
          interpretacion_males_pcos_obese[j] <-
            paste(inter_name,
                  "MENOR EN PCOS  OBESAS vs HOMBRES",
                  inter_i$malespcos_prueba_Obese)
          
        }
      }
      if (inter_i$p_value_SEXO_Obese < 0.05  &
          !is.na(inter_i$p_value_SEXO_Obese)) {
        if (inter_i$mu_SEXO_Obese < 0.05) {
          interpretacion_sex_obese[j] <-
            paste(
              inter_name,
              "MAYOR EN HOMBRES OBESOS vs MUJERES OBESAS",
              inter_i$sexo_prueba_Obese
            )
        } else{
          interpretacion_sex_obese[j] <-
            paste(
              inter_name,
              "MENOR EN HOMBRES  OBESOS vs MUJERES OBESAS",
              inter_i$sexo_prueba_Obese
            )
          
        }
      }
      
      if (inter_i$pcos_obesidad < 0.05 &
          !is.na(inter_i$pcos_obesidad)) {
        if (inter_i$pcos_diferencia_obesidad < 0) {
          interpretacion_pcos_obesidad[j] <-
            paste(inter_name, "MAYOR EN PCOS OBESAS vs PCOS NO OBESAS")
        } else{
          interpretacion_pcos_obesidad[j] <-
            paste(inter_name, "MENOR EN PCOS OBESAS vs PCOS NO OBESAS")
          
        }
        
      }
      if (inter_i$females_obesidad < 0.05 &
          !is.na(inter_i$females_obesidad)) {
        if (inter_i$females_diferencia_obesidad < 0) {
          interpretacion_fems_obesidad[j] <-
            paste(inter_name, "MAYOR EN MUJERES OBESAS vs MUJERES NO OBESAS")
        } else{
          interpretacion_fems_obesidad[j] <-
            paste(inter_name, "MENOR EN MUJERES  OBESAS vs MUERES NO OBESAS")
          
        }
        
      }
      if (inter_i$males_obesidad < 0.05 &
          !is.na(inter_i$males_obesidad)) {
        if (inter_i$males_diferencia_obesidad < 0) {
          interpretacion_males_obesidad[j] <-
            paste(inter_name, "MAYOR EN HOMBRES OBESAS vs HOMBRES NO OBESAS")
        } else{
          interpretacion_males_obesidad[j] <-
            paste(inter_name, "MENOR EN HOMBRES  OBESAS vs HOMBRES NO OBESAS")
          
        }
        
      }
      
      if ((
        inter_i$p_value_Females_PCOS_No_Obese < 0.05 &
        !is.na(inter_i$p_value_Females_PCOS_No_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_No_Obese)
      ) |
      (inter_i$p_value_SEXO_No_Obese > 0.05 |
       is.na(inter_i$p_value_SEXO_No_Obese)) &
      (
        inter_i$p_value_Females_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_Obese)
      ) &
      (inter_i$p_value_SEXO_Obese >  0.05 |
       is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_unica_grupo_obesidad[j] <-
          paste(inter_name, "UNICO ENTRE MUJERES SANAS Y PCOS")
        
      } else if ((
        inter_i$p_value_Females_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_No_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_No_Obese < 0.05 &
        !is.na(inter_i$p_value_Males_PCOS_No_Obese)
      ) &
      (inter_i$p_value_SEXO_No_Obese > 0.05 |
       is.na(inter_i$p_value_SEXO_No_Obese)) &
      (
        inter_i$p_value_Females_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_Obese)
      ) &
      (inter_i$p_value_SEXO_Obese >  0.05 |
       is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_unica_grupo_obesidad[j] <-
          paste(inter_name, "UNICO ENTRE HOMBRES SANOS Y PCOS")
        
        
      } else if ((
        inter_i$p_value_Females_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_No_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_No_Obese > 0.05 &
        is.na(inter_i$p_value_Males_PCOS_No_Obese)
      ) &
      (inter_i$p_value_SEXO_No_Obese < 0.05  &
       !is.na(inter_i$p_value_SEXO_No_Obese)) &
      (
        inter_i$p_value_Females_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_Obese)
      ) &
      (inter_i$p_value_SEXO_Obese >  0.05 |
       is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_unica_grupo_obesidad[j] <-
          paste(inter_name, "ENTRE SEXO SANOS")
        
        
      } else if ((
        inter_i$p_value_Females_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_No_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_No_Obese > 0.05  |
        is.na(inter_i$p_value_Males_PCOS_No_Obese)
      ) &
      (inter_i$p_value_SEXO_No_Obese > 0.05 |
       is.na(inter_i$p_value_SEXO_No_Obese)) &
      (
        inter_i$p_value_Females_PCOS_Obese < 0.05 &
        is.na(inter_i$p_value_Females_PCOS_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_Obese)
      ) &
      (inter_i$p_value_SEXO_Obese >  0.05 |
       is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_unica_grupo_obesidad[j] <-
          paste(inter_name, "ENTRE MUJERES OBESAS Y PCOS")
        
        
      } else if ((
        inter_i$p_value_Females_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_No_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_No_Obese)
      ) &
      (inter_i$p_value_SEXO_No_Obese > 0.05 |
       is.na(inter_i$p_value_SEXO_No_Obese)) &
      (
        inter_i$p_value_Females_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_Obese < 0.05 &
        !is.na(inter_i$p_value_Males_PCOS_Obese)
      ) &
      (inter_i$p_value_SEXO_Obese >  0.05 |
       is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_unica_grupo_obesidad[j] <-
          paste(inter_name, "ENTRE HOMBRES Y PCOS OBESOS")
        
        
      } else if ((
        inter_i$p_value_Females_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_No_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_No_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_No_Obese)
      ) &
      (inter_i$p_value_SEXO_No_Obese > 0.05 |
       is.na(inter_i$p_value_SEXO_No_Obese)) &
      (
        inter_i$p_value_Females_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Females_PCOS_Obese)
      ) &
      (
        inter_i$p_value_Males_PCOS_Obese > 0.05 |
        is.na(inter_i$p_value_Males_PCOS_Obese)
      ) &
      (inter_i$p_value_SEXO_Obese <  0.05 &
       !is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_unica_grupo_obesidad[j] <-
          paste(inter_name, "ENTRE SEXO OBESO")
        
        
      }
      if (inter_i$pcos_obesidad < 0.05 &
          !is.na(inter_i$pcos_obesidad) & (
            inter_i$females_obesidad > 0.05 | is.na(inter_i$females_obesidad) |
            inter_i$males_obesidad > 0.05 |
            is.na(inter_i$males_obesidad)
          )) {
        interpretacion_unica_obesidad_grupo[j] <-
          paste(inter_name, "SOLO SE VE EN PCOS")
        
        
      } else if ((inter_i$pcos_obesidad > 0.05 |
                  is.na(inter_i$pcos_obesidad)) &
                 (inter_i$females_obesidad < 0.05 &
                  !is.na(inter_i$females_obesidad)) &
                 (inter_i$males_obesidad > 0.05 |
                  is.na(inter_i$males_obesidad))) {
        interpretacion_unica_obesidad_grupo[j] <-
          paste(inter_name, "SOLO SE VE EN MUJERES")
        
        
      } else if ((inter_i$pcos_obesidad > 0.05 |
                  is.na(inter_i$pcos_obesidad)) &
                 (inter_i$females_obesidad > 0.05 |
                  is.na(inter_i$females_obesidad)) &
                 (inter_i$males_obesidad < 0.05 &
                  !is.na(inter_i$males_obesidad))) {
        interpretacion_unica_obesidad_grupo[j] <-
          paste(inter_name, "SOLO SE VE EN HOMBRES")
        
        
      }
      if ((inter_i$pcos_obesidad < 0.05 &
           !is.na(inter_i$pcos_obesidad)) &
          (inter_i$females_obesidad > 0.05 |
           is.na(inter_i$females_obesidad)) &
          (inter_i$males_obesidad > 0.05 |
           is.na(inter_i$males_obesidad)) &
          (
            inter_i$p_value_Females_PCOS_No_Obese < 0.05 &
            !is.na(inter_i$p_value_Females_PCOS_No_Obese)
          ) &
          (
            inter_i$p_value_Females_PCOS_Obese < 0.05 &
            !is.na(inter_i$p_value_Females_PCOS_Obese)
          ) &
          (
            inter_i$p_value_Males_PCOS_No_Obese < 0.05 &
            !is.na(inter_i$p_value_Males_PCOS_No_Obese)
          ) &
          (
            inter_i$p_value_Males_PCOS_Obese < 0.05 &
            !is.na(inter_i$p_value_Males_PCOS_Obese)
          ) &
          (inter_i$p_value_SEXO_No_Obese > 0.05 |
           is.na(inter_i$p_value_SEXO_No_Obese)) &
          (inter_i$p_value_SEXO_Obese > 0.05 |
           is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_final[j] <-
          "La interaccion solamente se ve en la obesidad PCOS"
        
        
      } else if ((inter_i$pcos_obesidad > 0.05 |
                  is.na(inter_i$pcos_obesidad)) &
                 (inter_i$females_obesidad < 0.05 &
                  !is.na(inter_i$females_obesidad)) &
                 (inter_i$males_obesidad > 0.05 &
                  is.na(inter_i$males_obesidad)) &
                 (
                   inter_i$p_value_Females_PCOS_No_Obese < 0.05 &
                   !is.na(inter_i$p_value_Females_PCOS_No_Obese)
                 ) &
                 (
                   inter_i$p_value_Females_PCOS_Obese < 0.05 &
                   !is.na(inter_i$p_value_Females_PCOS_Obese)
                 ) &
                 (
                   inter_i$p_value_Males_PCOS_No_Obese > 0.05 |
                   is.na(inter_i$p_value_Males_PCOS_No_Obese)
                 ) &
                 (
                   inter_i$p_value_Males_PCOS_Obese > 0.05 |
                   is.na(inter_i$p_value_Males_PCOS_Obese)
                 ) &
                 (inter_i$p_value_SEXO_No_Obese < 0.05 &
                  !is.na(inter_i$p_value_SEXO_No_Obese)) &
                 (inter_i$p_value_SEXO_Obese < 0.05 &
                  !is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_final[j] <-
          "La interaccion solamente se ve en la obesidad MUJERES"
        
        
      } else if ((inter_i$pcos_obesidad > 0.05 |
                  is.na(inter_i$pcos_obesidad)) &
                 (inter_i$females_obesidad > 0.05 |
                  is.na(inter_i$females_obesidad)) &
                 (inter_i$males_obesidad < 0.05 &
                  !is.na(inter_i$males_obesidad)) &
                 (
                   inter_i$p_value_Females_PCOS_No_Obese > 0.05 |
                   is.na(inter_i$p_value_Females_PCOS_No_Obese)
                 ) &
                 (
                   inter_i$p_value_Females_PCOS_Obese > 0.05 |
                   is.na(inter_i$p_value_Females_PCOS_Obese)
                 ) &
                 (
                   inter_i$p_value_Males_PCOS_No_Obese < 0.05 &
                   !is.na(inter_i$p_value_Males_PCOS_No_Obese)
                 ) &
                 (
                   inter_i$p_value_Males_PCOS_Obese < 0.05 &
                   !is.na(inter_i$p_value_Males_PCOS_Obese)
                 ) &
                 (inter_i$p_value_SEXO_No_Obese < 0.05 &
                  !is.na(inter_i$p_value_SEXO_No_Obese)) &
                 (inter_i$p_value_SEXO_Obese < 0.05 &
                  !is.na(inter_i$p_value_SEXO_Obese))) {
        interpretacion_final[j] <-
          "La interaccion solamente se ve en la obesidad HOMBRES"
        
        
      }
      
      
    }
    
    
  }
  
  
  matriz_pre_resultados <- data.frame(
    cbind(
      p_value_obesidad = glm_modelo$obesidad,
      p_value_obesidad2 = glm_modelo$p.value_OBESIDAD,
      p_value_females = glm_modelo$p_value_Females_PCOS,
      p_value_pcosmales = glm_modelo$p_value_Males_PCOS,
      p_value_sexo = glm_modelo$p_value_SEXO,
      p_value_grupo = glm_modelo$grupo,
      p_value_females_no_obese = glm_modelo$p_value_Females_PCOS_No_Obese,
      p_value_females_obese = glm_modelo$p_value_Females_PCOS_Obese,
      p_value_pcosmales_no_obese = glm_modelo$p_value_Males_PCOS_No_Obese,
      p_value_pcosmales_obese = glm_modelo$p_value_Males_PCOS_Obese,
      p_value_sexo_no_obese = glm_modelo$p_value_SEXO_No_Obese,
      p_value_sexo_obese = glm_modelo$p_value_SEXO_Obese,
      p_value_pcos_obesidad = glm_modelo$pcos_obesidad,
      p_value_females_obesidad = glm_modelo$females_obesidad,
      p_value_males_obesidad = glm_modelo$males_obesidad,
      interpretacion_obesidad = interpretacion_obesidad,
      interpretacion_fems_pcos = interpretacion_fems_pcos,
      interpretacion_males_pcos = interpretacion_males_pcos,
      interpretacion_sexo = interpretacion_sex,
      interpretacion_grupo = interpretacion_unica_grupo,
      intepretacion_females_no_obese = interpretacion_fems_pcos_no_obese,
      intepretacion_females_obese = interpretacion_fems_pcos_obese,
      interpretracion_pcosmales_no_obese = interpretacion_males_pcos_no_obese,
      interpretracion_pcosmales_obese = interpretacion_males_pcos_obese,
      interpretacion_sex_no_obese = interpretacion_sex_no_obese,
      interpretacion_sex_obese=interpretacion_sex_obese,
      interpretacion_females_obesidad = interpretacion_fems_obesidad,
      interpretacion_males_obesidad = interpretacion_males_obesidad,
      intepretracion_pcos_obesidad = interpretacion_pcos_obesidad,
      interpretacion_unica_obesidad_grupo = interpretacion_unica_obesidad_grupo,
      interpretacion_unica_grupo_obesdiad = interpretacion_unica_grupo_obesidad,
      interpretracion_final = interpretacion_final
      
    )
  )
  
  matriz_resultados <- matriz_pre_resultados
  
  matriz_resultados$interpretacion_obesidad <- ifelse(is.na(matriz_resultados$p_value_obesidad2),
                                                      retorno_tipo_prueba$tipo.prueba_OBESIDAD,
                                                      matriz_resultados$interpretacion_obesidad)
  matriz_resultados$interpretacion_fems_pcos <- ifelse(is.na(matriz_resultados$p_value_females),
                                                       retorno_tipo_prueba$femalespcos_prueba,
                                                       matriz_resultados$interpretacion_fems_pcos)
  
  matriz_resultados$interpretacion_males_pcos <- ifelse(is.na(matriz_resultados$p_value_pcosmales),
                                                        retorno_tipo_prueba$malespcos_prueba,
                                                        matriz_resultados$interpretacion_males_pcos)
  
  matriz_resultados$interpretacion_sexo <- ifelse(is.na(matriz_resultados$p_value_sexo),
                                                  retorno_tipo_prueba$sexo_prueba,
                                                  matriz_resultados$interpretacion_sexo)
  
  matriz_resultados$intepretacion_females_no_obese <- ifelse(is.na(matriz_resultados$p_value_females_no_obese),
                                                             retorno_tipo_prueba$femalespcos_prueba_No_Obese,
                                                             matriz_resultados$intepretacion_females_no_obese)
  
  
  matriz_resultados$interpretracion_pcosmales_no_obese <- ifelse(is.na(matriz_resultados$p_value_pcosmales_no_obese),
                                                                 retorno_tipo_prueba$malespcos_prueba_No_Obese,
                                                                 matriz_resultados$interpretracion_pcosmales_no_obese)
  
  
  
  matriz_resultados$interpretacion_sex_no_obese <- ifelse(is.na(matriz_resultados$p_value_sexo_no_obese),
                                                          retorno_tipo_prueba$sexo_prueba_No_Obese,
                                                          matriz_resultados$interpretacion_sex_no_obese)
  
  matriz_resultados$intepretacion_females_obese <- ifelse(is.na(matriz_resultados$p_value_females_obese),
                                                          retorno_tipo_prueba$femalespcos_prueba_Obese,
                                                          matriz_resultados$intepretacion_females_obese)
  
  
  matriz_resultados$interpretracion_pcosmales_obese <- ifelse(is.na(matriz_resultados$p_value_pcosmales_obese),
                                                              retorno_tipo_prueba$malespcos_prueba_Obese,
                                                              matriz_resultados$interpretracion_pcosmales_obese)
  
  
  
  matriz_resultados$interpretacion_sex_obese <- ifelse(is.na(matriz_resultados$p_value_sexo_obese),
                                                       retorno_tipo_prueba$sexo_prueba_Obese,
                                                       matriz_resultados$interpretacion_sex_obese)
  
  
  matriz_resultados$intepretracion_pcos_obesidad<- ifelse(is.na(matriz_resultados$p_value_pcos_obesidad),
                                                          retorno_tipo_prueba$pcos_prueba_obesidad,
                                                          matriz_resultados$intepretracion_pcos_obesidad)
  
  matriz_resultados$interpretacion_females_obesidad <- ifelse(is.na(matriz_resultados$p_value_females_obesidad),
                                                              retorno_tipo_prueba$females_prueba_obesidad,
                                                              matriz_resultados$interpretacion_females_obesidad)
  
  matriz_resultados$interpretacion_males_obesidad <- ifelse(is.na(matriz_resultados$p_value_males_obesidad),
                                                            retorno_tipo_prueba$males_prueba_obesidad,
                                                            matriz_resultados$interpretacion_males_obesidad)
  
  
  
  
  
  
  rownames(matriz_resultados) <- colnames(X)
  
  return(matriz_resultados)
  
}


