

library(ggplot2)
library(dplyr)
library(reshape2)
library(lmPerm)
library(data.table)
library(DescTools)


source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
panova <-
  function(x)
    summary(aovp(
      x ~ grupo * obesidad,
      perm = "Exact",
      contrasts = list(grupo = contr.sum(n = 3), obesidad = contr.sum(n = 2)),
      center = T,
      seqs = F
    ))[[1]]$`Pr(Prob)`[-4]

triada <- function(X_list, g_, variable) {
  X.list <- split(X_list[[variable]], g_[[variable]])
  X <- X.list[[1]]
  Y <- X.list[[2]]
  p <- ncol(X)
  return(unlist(lapply(1:p, function(x)
    HodgesLehmann(X[, x], Y[, x]))))
  
  
}

triada_test <- function(X_list, g_, variable) {
  X.list <- split(X_list[[variable]], g_[[variable]])
  X <- X.list[[1]]
  Y <- X.list[[2]]
  p <- ncol(X)
  return(unlist(lapply(1:p, function(x)
    wilcox.test(X[, x], Y[, x], exact = T, correct = T)$p.value)))
  
  
}

get_wilcox_stat <- function(X, grupo, obesidad) {
  X.grupo <- split(X, grupo)
  X.obesidad <- split(X, obesidad)
  
  
  X_list <-
    list(
      Female_vs_PCOS_O = X[grupo != "Male" & obesidad == "Obese", ],
      Female_vs_Male_O = X[grupo != "PCOS" &
                             obesidad == "Obese", ],
      PCOS_vs_Male_O = X[grupo != "Male"  &
                           obesidad == "Obese", ],
      Female_vs_PCOS_N = X[grupo != "Male" &
                             obesidad == "No Obese", ],
      Female_vs_Male_N = X[grupo != "PCOS" &
                             obesidad == "No Obese", ],
      PCOS_vs_Male_N = X[grupo != "Male"  &
                           obesidad == "No Obese", ],
      NO_Female = X[grupo == "Female", ],
      NO_PCOS = X[grupo == "PCOS", ],
      NO_Male = X[grupo == "Male", ],
      Female_vs_PCOS = X[grupo != "Male" , ],
      Female_vs_Male = X[grupo != "PCOS" , ],
      PCOS_vs_Male = X[grupo != "Male" , ]
    )
  
  g_ <-
    list(
      Female_vs_PCOS_O = as.factor(as.character(grupo[grupo != "Male" &
                                                        obesidad == "Obese"])),
      Female_vs_Male_O =  as.factor(as.character(grupo[grupo !=
                                                         "PCOS" & obesidad == "Obese"])),
      PCOS_vs_Male_O =  as.factor(as.character(grupo[grupo !=
                                                       "Male"  & obesidad == "Obese"])),
      Female_vs_PCOS_N =  as.factor(as.character(grupo[grupo !=
                                                         "Male" & obesidad == "No Obese"])),
      Female_vs_Male_N =  as.factor(as.character(grupo[grupo !=
                                                         "PCOS" & obesidad == "No Obese"])),
      PCOS_vs_Male_N =  as.factor(as.character(grupo[grupo !=
                                                       "Male"  & obesidad == "No Obese"])),
      NO_Female =  as.factor(as.character(obesidad[grupo == "Female"])),
      NO_PCOS = as.factor(as.character(obesidad[grupo == "PCOS"])),
      NO_Male =  as.factor(as.character(obesidad[grupo == "Male"])),
      Female_vs_PCOS = as.factor(as.character(grupo[grupo != "Male"])),
      Female_vs_Male = as.factor(as.character(grupo[grupo != "PCOS"])),
      PCOS_vs_Male = as.factor(as.character(grupo[grupo != "Female"]))
    )
  
  
  
  
  retorno <-
    data.frame(
      Female_vs_PCOS = triada(X_list, g_, "Female_vs_PCOS"),
      Female_vs_Male = triada(X_list = X_list, g_, "Female_vs_Male"),
      PCOS_vs_Male = triada(X_list = X_list, g_, "PCOS_vs_Male"),
      NO_Female_vs_PCOS = triada(X_list = X_list, g_, "Female_vs_PCOS_N"),
      NO_Female_vs_Male = triada(X_list = X_list, g_, "Female_vs_Male_N"),
      NO_PCOS_Vs_Male = triada(X_list = X_list, g_, "PCOS_vs_Male_N"),
      O_Female_vs_PCOS = triada(X_list = X_list, g_, "Female_vs_PCOS_O"),
      O_Female_vs_Male = triada(X_list = X_list, g_, "Female_vs_Male_O"),
      O_PCOS_Vs_Male = triada(X_list = X_list, g_, "PCOS_vs_Male_O"),
      PCOS = triada(X_list = X_list, g_, variable = "NO_PCOS"),
      Female = triada(X_list = X_list, g_, "NO_Female"),
      Male = triada(X_list = X_list, g_, "NO_Male")
    )
  rownames(retorno) <- colnames(X_list$Female_vs_PCOS_O)
  
  
  return(retorno)
  
  
}


get_wilcox_p <- function(X, grupo, obesidad) {
  X.grupo <- split(X, grupo)
  X.obesidad <- split(X, obesidad)
  
  
  X_list <-
    list(
      Female_vs_PCOS_O = X[grupo != "Male" & obesidad == "Obese", ],
      Female_vs_Male_O = X[grupo != "PCOS" &
                             obesidad == "Obese", ],
      PCOS_vs_Male_O = X[grupo != "Male"  &
                           obesidad == "Obese", ],
      Female_vs_PCOS_N = X[grupo != "Male" &
                             obesidad == "No Obese", ],
      Female_vs_Male_N = X[grupo != "PCOS" &
                             obesidad == "No Obese", ],
      PCOS_vs_Male_N = X[grupo != "Male"  &
                           obesidad == "No Obese", ],
      NO_Female = X[grupo == "Female", ],
      NO_PCOS = X[grupo == "PCOS", ],
      NO_Male = X[grupo == "Male", ],
      Female_vs_PCOS = X[grupo != "Male" , ],
      Female_vs_Male = X[grupo != "PCOS" , ],
      PCOS_vs_Male = X[grupo != "Male" , ]
    )
  
  g_ <-
    list(
      Female_vs_PCOS_O = as.factor(as.character(grupo[grupo != "Male" &
                                                        obesidad == "Obese"])),
      Female_vs_Male_O =  as.factor(as.character(grupo[grupo != "PCOS" &
                                                         obesidad == "Obese"])),
      PCOS_vs_Male_O =  as.factor(as.character(grupo[grupo != "Male"  &
                                                       obesidad == "Obese"])),
      Female_vs_PCOS_N =  as.factor(as.character(grupo[grupo != "Male" &
                                                         obesidad == "No Obese"])),
      Female_vs_Male_N =  as.factor(as.character(grupo[grupo != "PCOS" &
                                                         obesidad == "No Obese"])),
      PCOS_vs_Male_N =  as.factor(as.character(grupo[grupo != "Male"  &
                                                       obesidad == "No Obese"])),
      NO_Female =  as.factor(as.character(obesidad[grupo == "Female"])),
      NO_PCOS = as.factor(as.character(obesidad[grupo == "PCOS"])),
      NO_Male =  as.factor(as.character(obesidad[grupo == "Male"])),
      Female_vs_PCOS = as.factor(as.character(grupo[grupo != "Male"])),
      Female_vs_Male = as.factor(as.character(grupo[grupo != "PCOS"])),
      PCOS_vs_Male = as.factor(as.character(grupo[grupo != "Female"]))
    )
  
  
  
  
  
  
  retorno <-
    data.frame(
      Female_vs_PCOS = triada_test(X_list, g_, "Female_vs_PCOS"),
      Female_vs_Male = triada_test(X_list = X_list, g_, "Female_vs_Male"),
      PCOS_vs_Male = triada_test(X_list = X_list, g_, "PCOS_vs_Male"),
      NO_Female_vs_PCOS = triada_test(X_list = X_list, g_, "Female_vs_PCOS_N"),
      NO_Female_vs_Male = triada_test(X_list = X_list, g_, "Female_vs_Male_N"),
      NO_PCOS_Vs_Male = triada_test(X_list = X_list, g_, "PCOS_vs_Male_N"),
      O_Female_vs_PCOS = triada_test(X_list = X_list, g_, "Female_vs_PCOS_O"),
      O_Female_vs_Male = triada_test(X_list = X_list, g_, "Female_vs_Male_O"),
      O_PCOS_Vs_Male = triada_test(X_list = X_list, g_, "PCOS_vs_Male_O"),
      PCOS = triada_test(X_list = X_list, g_, variable = "NO_PCOS"),
      Female = triada_test(X_list = X_list, g_, "NO_Female"),
      Male = triada_test(X_list = X_list, g_, "NO_Male")
    )
  rownames(retorno) <- colnames(X_list$Female_vs_PCOS_O)
  
  
  return(retorno)
  
}


get_anova_p <- function(X) {
  panova_Res <- t(apply(X, 2, panova))
  colnames(panova_Res) <- c("Grupo", "Obesidad", "Interaccion")
  panova_Res <- as.data.frame(panova_Res)
  return(panova_Res)
}
datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")

fun_correct.BH <- function(X_) {
  X_[, 1:3] <- t(apply(X_[, 1:3], 1, function(z)
    p.adjust(z, "BH")))
  X_[, 4:ncol(X_)] <-
    t(apply(X_[, 4:ncol(X_)], 1, function(z)
      p.adjust(z, "BH")))
  
  X_.c <- apply(X_, 2, function(z)
    p.adjust(z, "BH"))
  return(X_.c)
}



inter_niveles <- function(niveles,ob) {
  Female_vs_PCOS <- niveles[1]
  Female_vs_Male <- niveles[2]
  PCOS_vs_Male <- niveles[3]
  # Definir los valores
  
  # Evaluar si la obesidad es menor a 0.05
  if(!is.na(ob)) {
    o <- paste(ob)
    # Identificar cuáles valores son menores a 0.05
    if(!is.na(Female_vs_PCOS)) {
      fp <- (Female_vs_PCOS)
    }else{
      fp <- "-"
      
    }
    if(!is.na(Female_vs_Male)) {
      fm <- (Female_vs_Male)
      
    }else{
      fm <- "-"
    }
    if(!is.na(PCOS_vs_Male)) {
      pm <- (PCOS_vs_Male)
    }else{
      pm <- "-"
    }
  } else {
    o <- "No hay diferencias en la obesidad, pero si en "
    # Identificar cuáles valores son mayores o iguales a 0.05
    if(!is.na(Female_vs_PCOS)) {
      fp <- (Female_vs_PCOS)
    }else{
      fp <- "-"
    }
    if(!is.na(Female_vs_Male)) {
      fm <- (Female_vs_Male)
      
    }else{
      fm <- "-"
    }
    if(!is.na(PCOS_vs_Male)) {
      pm <- (PCOS_vs_Male)
    }else{
      pm <- "-"
    }
  }
  
  retorno <- paste(o,fp,fm,pm)
  return(retorno)
  
}

inter_Ob <- function(niveles,ob) {
  Female_vs_PCOS <- niveles[1]
  Female_vs_Male <- niveles[2]
  PCOS_vs_Male <- niveles[3]
  # Definir los valores
  
  # Evaluar si la obesidad es menor a 0.05
  if(!is.na(ob)) {
    o <- paste(ob)
    # Identificar cuáles valores son menores a 0.05
    if(!is.na(Female_vs_PCOS)) {
      fp <- (Female_vs_PCOS)
    }else{
      fp <- "-"
      
    }
    if(!is.na(Female_vs_Male)) {
      fm <- (Female_vs_Male)
      
    }else{
      fm <- "-"
    }
    if(!is.na(PCOS_vs_Male)) {
      pm <- (PCOS_vs_Male)
    }else{
      pm <- "-"
    }
  } else {
    o <- "No hay diferencias en la obesidad, pero si en "
    # Identificar cuáles valores son mayores o iguales a 0.05
    if(!is.na(Female_vs_PCOS)) {
      fp <- (Female_vs_PCOS)
    }else{
      fp <- "-"
    }
    if(!is.na(Female_vs_Male)) {
      fm <- (Female_vs_Male)
      
    }else{
      fm <- "-"
    }
    if(!is.na(PCOS_vs_Male)) {
      pm <- (PCOS_vs_Male)
    }else{
      pm <- "-"
    }
  }
  
  retorno <- paste(fp,fm,pm)
  return(retorno)
  
}

aux_gsub <- function(texto){
  nuevo_texto <- gsub("-", ". ", texto)
  
  # Eliminar guiones consecutivos
  nuevo_texto <- gsub("\\s*-+\\s*", ". ", nuevo_texto)
  return(nuevo_texto)
}








