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
  X.list <-
    lapply(split(X_list[[variable]],
                 g_[[variable]]),
           matrix, ncol = ncol(X_list[[variable]]))
  X <- X.list[[1]]
  Y <- X.list[[2]]
  p <- ncol(X)
  return(unlist(lapply(1:p, function(x)
    HodgesLehmann(X[, x], Y[, x]))))
  
  
}

triada_test <- function(X_list, g_, variable) {
  X.list <-
    lapply(split(X_list[[variable]],
                 g_[[variable]]),
           matrix, ncol = ncol(X_list[[variable]]))
  X <- X.list[[1]]
  Y <- X.list[[2]]
  p <- ncol(X)
  return(unlist(lapply(1:p, function(x)
    wilcox.test(X[, x], Y[, x], exact = T, correct = T)$p.value)))
  
  
}

get_wilcox_stat <- function(X, grupo, obesidad) {
  X_list <-
    list(
      Female_vs_PCOS_O = X[grupo != "Male" & obesidad == "Obese",],
      Female_vs_Male_O = X[grupo != "PCOS" &
                             obesidad == "Obese",],
      PCOS_vs_Male_O = X[grupo != "Female"  &
                           obesidad == "Obese",],
      Female_vs_PCOS_N = X[grupo != "Male" &
                             obesidad == "No Obese",],
      Female_vs_Male_N = X[grupo != "PCOS" &
                             obesidad == "No Obese",],
      PCOS_vs_Male_N = X[grupo != "Female"  &
                           obesidad == "No Obese",],
      NO_Female = X[grupo == "Female",],
      NO_PCOS = X[grupo == "PCOS",],
      NO_Male = X[grupo == "Male",],
      Female_vs_PCOS = X[grupo != "Male" ,],
      Female_vs_Male = X[grupo != "PCOS" ,],
      PCOS_vs_Male = X[grupo != "Female" ,]
    )
  
  g_ <-
    list(
      Female_vs_PCOS_O = as.factor(as.character(grupo[grupo != "Male" &
                                                        obesidad == "Obese"])),
      Female_vs_Male_O =  as.factor(as.character(grupo[grupo != "PCOS" &
                                                         obesidad == "Obese"])),
      PCOS_vs_Male_O =  as.factor(as.character(grupo[grupo != "Female"  &
                                                       obesidad == "Obese"])),
      Female_vs_PCOS_N =  as.factor(as.character(grupo[grupo != "Male" &
                                                         obesidad == "No Obese"])),
      Female_vs_Male_N =  as.factor(as.character(grupo[grupo != "PCOS" &
                                                         obesidad == "No Obese"])),
      PCOS_vs_Male_N =  as.factor(as.character(grupo[grupo != "Female"  &
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
      Female_vs_PCOS_O = X[grupo != "Male" & obesidad == "Obese",],
      Female_vs_Male_O = X[grupo != "PCOS" &
                             obesidad == "Obese",],
      PCOS_vs_Male_O = X[grupo != "Female"  &
                           obesidad == "Obese",],
      Female_vs_PCOS_N = X[grupo != "Male" &
                             obesidad == "No Obese",],
      Female_vs_Male_N = X[grupo != "PCOS" &
                             obesidad == "No Obese",],
      PCOS_vs_Male_N = X[grupo != "Female"  &
                           obesidad == "No Obese",],
      NO_Female = X[grupo == "Female",],
      NO_PCOS = X[grupo == "PCOS",],
      NO_Male = X[grupo == "Male",],
      Female_vs_PCOS = X[grupo != "Male" ,],
      Female_vs_Male = X[grupo != "PCOS" ,],
      PCOS_vs_Male = X[grupo != "Female" ,]
    )
  
  g_ <-
    list(
      Female_vs_PCOS_O = as.factor(as.character(grupo[grupo != "Male" &
                                                        obesidad == "Obese"])),
      Female_vs_Male_O =  as.factor(as.character(grupo[grupo != "PCOS" &
                                                         obesidad == "Obese"])),
      PCOS_vs_Male_O =  as.factor(as.character(grupo[grupo != "Female"  &
                                                       obesidad == "Obese"])),
      Female_vs_PCOS_N =  as.factor(as.character(grupo[grupo != "Male" &
                                                         obesidad == "No Obese"])),
      Female_vs_Male_N =  as.factor(as.character(grupo[grupo != "PCOS" &
                                                         obesidad == "No Obese"])),
      PCOS_vs_Male_N =  as.factor(as.character(grupo[grupo != "Female"  &
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

fun_correct.BH <- function(X_) {
  X_[, 1:3] <- t(apply(X_[, 1:3], 1, function(z)
    p.adjust(z, "BH")))
  X_[, 4:ncol(X_)] <-
    t(apply(X_[, 4:ncol(X_)], 1, function(z)
      p.adjust(z, "BH")))
  
  # X_.c <- apply(X_, 2, function(z)
  # p.adjust(z, "BH"))
  return(X_)
}




concatenando_grupo <- function(x) {

  Grupo <- x[1]
  Obesidad <- x[2]
  Interaccion <- x[3]
  Female_vs_PCOS <- x[4]
  Female_vs_Male <-x[5]
  PCOS_vs_Male <-x[6]
  NO_Female_vs_PCOS <- x[7]
  NO_Female_vs_Male <- x[8]
  NO_PCOS_Vs_Male <- x[9]
  O_Female_vs_PCOS <-x[10]
  O_Female_vs_Male <- x[11]
  O_PCOS_Vs_Male <-x[12]
  PCOS <- x[13]
  Female <-x[14]
  Male <-x[15]
  # Definir los valores
  if(!is.na(Obesidad)){
    ob <- Obesidad
  }else{
    ob <- ""
  }
  # Evaluar si el grupo es menor a 0.05
  if (!is.na(Grupo)) {
    o <- Grupo
    # Identificar cuáles valores son menores a 0.05
    if (!is.na(Female_vs_PCOS)) {
      fp <- (Female_vs_PCOS)
    } else{
      fp <- "-"
      
    }
    if (!is.na(Female_vs_Male)) {
      fm <- (Female_vs_Male)
      
    } else{
      fm <- "-"
    }
    if (!is.na(PCOS_vs_Male)) {
      pm <- (PCOS_vs_Male)
    } else{
      pm <- "-"
    }
  } else {
    o <- ""
    fp <- ""
    fm <- ""
    pm <- ""
  }
  
  if(!is.na(Interaccion)){
    i <- Interaccion
    if(!is.na(NO_Female_vs_PCOS)){
      fp.no <- NO_Female_vs_PCOS
    }else{
      fp.no <-""
    }
    if (!is.na(NO_Female_vs_Male)) {
      
      fm.no <- NO_Female_vs_Male
    }else{
      fm.no <- ""
    }
    if (!is.na(NO_PCOS_Vs_Male)) {
      pm.no <- NO_PCOS_Vs_Male
    }else{
      pm.no <- ""
    }
    if(!is.na(O_Female_vs_PCOS)){
      fp.o <- O_Female_vs_PCOS
    }else{
      fp.o <- ""
    }
    if (!is.na(O_Female_vs_Male)) {
      
      fm.o <-O_Female_vs_Male
      
    }else{
      fm.o <-""
    }
    if(!is.na(O_PCOS_Vs_Male)){
      pm.o <-O_PCOS_Vs_Male
    }else{
      pm.o <- ""
    }
    if(!is.na(PCOS)){
      p <- PCOS
    }else{
      p<-""
    }
    if(!is.na(Female)){
      f <- Female
    }else{
      f<- ""
    }
    if(!is.na(Male)){
      m<-Male
    }else{
      m<-""
    }
    
  }else{
    i <- ""
    fp.no <- ""
    fm.no <- ""
    pm.no <- ""
    fp.o <- ""
    fm.o <- ""
    pm.o <- ""
    p <- ""
    f <- ""
    m <- ""
  }
  
  
  retorno <- paste(o,"GRUPO(",fm,fp,pm,")",i,
                   "INTERACCION[",
                   "NO OBESIDAD(",fp.no,fm.no,pm.no,")",
                   "OBESIDAD(",fp.o,fm.o,pm.o,")",
                   "PCOS(",p,")",
                   "Female(",f,")",
                   "Male(",m,")","]")
  if(ob!="" & o !="" & i !=""){
    retorno <- paste(ob,o,"GRUPO(",fm,fp,pm,")",i,
                     "INTERACCION[",
                     "NO OBESIDAD(",fp.no,fm.no,pm.no,")",
                     "OBESIDAD(",fp.o,fm.o,pm.o,")",
                     "PCOS(",p,")",
                     "Female(",f,")",
                     "Male(",m,")","]")
    
  }else if(ob!="" & o !="" & i==""){
    ## obesidad y grupo no interaccion
    retorno <- paste(ob,o,"GRUPO(",fm,fp,pm,")")
    
  }else if(ob!="" & o =="" & i!=""){
    ## obesidad e interaccion no grupo
    retorno <- paste(ob,i,
                     "INTERACCION[",
                     "NO OBESIDAD(",fp.no,fm.no,pm.no,")",
                     "OBESIDAD(",fp.o,fm.o,pm.o,")",
                     "PCOS(",p,")",
                     "Female(",f,")",
                     "Male(",m,")","]")
  }else if(ob=="" & o !="" & i!=""){
    ## grupo e interaccion no obesidad
    
    retorno <- paste(o,"GRUPO(",fm,fp,pm,")",i,
                     "INTERACCION[",
                     "NO OBESIDAD(",fp.no,fm.no,pm.no,")",
                     "OBESIDAD(",fp.o,fm.o,pm.o,")",
                     "PCOS(",p,")",
                     "Female(",f,")",
                     "Male(",m,")","]")
    
  }else if(ob!="" & o !="" & i==""){
    ## solo interaccion
    
    retorno <- paste(i,
                     "INTERACCION[",
                     "NO OBESIDAD(",fp.no,fm.no,pm.no,")",
                     "OBESIDAD(",fp.o,fm.o,pm.o,")",
                     "PCOS(",p,")",
                     "Female(",f,")",
                     "Male(",m,")","]")
    
  }else if(ob=="" & o !="" & i==""){
    ## solo grupo
    retorno <- paste(o,"GRUPO(",fm,fp,pm,")")
    
  }else if(ob!="" & o =="" & i==""){
    ## solo obesidad
    retorno <- ob
  }else if(ob=="" & o =="" & i==""){
    ## solo ineteraccion
    retorno <- ""
  }else if(ob!="" & o !="" & i!=""){
    ## todos
    retorno <- paste(ob,o,"GRUPO(",fm,fp,pm,")",i,
                     "INTERACCION[",
                     "NO OBESIDAD(",fp.no,fm.no,pm.no,")",
                     "OBESIDAD(",fp.o,fm.o,pm.o,")",
                     "PCOS(",p,")",
                     "Female(",f,")",
                     "Male(",m,")","]")
  }
  return(retorno)
  
}





aux_gsub <- function(texto) {
  nuevo_texto <- gsub("-", ". ", texto)
  
  # Eliminar guiones consecutivos
  nuevo_texto <- gsub("\\s*-+\\s*", ". ", nuevo_texto)
  return(nuevo_texto)
}

analyze_data <-  function(X,grupo,obesidad,correccion=2){
  ## correccion 1 2 3 diferentes tipos de correcion la mas adecuada 3
  res_w_p <- (get_wilcox_p(X, grupo, obesidad))
  res_aovp <- get_anova_p(X)
  res_w_t <-
    cbind(NO_vs_O = (apply(X, 2, function(x)
      HodgesLehmann(x[which(obesidad == "Obese")], x[which(obesidad == "No Obese")]))), get_wilcox_stat(X, grupo, obesidad))
  
  res_aovp.p_adj <- apply(res_aovp, 2, function(x)
    p.adjust(x, "BH"))
  
  res_w_p.adj_T <- fun_correct.BH(res_w_p)
  if (correccion == 1) {
    res_p <- as.data.frame(cbind(res_aovp, res_w_p))
  }else if(correccion==2){
    res_p <- as.data.frame(cbind(res_aovp, res_w_p.adj_T))
    
  }else{
    res_p <- as.data.frame(cbind(res_aovp.p_adj, apply(res_w_p.adj_T,2,function(x) p.adjust(x,"BH"))))
    
  }
  
 
  
  Female_vs_PCOS <- ifelse(
    res_p$Female_vs_PCOS < 0.05 &
      res_p$Grupo < 0.05  &
      res_w_t$Female_vs_PCOS != 0,
    ifelse(
      res_w_t$Female_vs_PCOS > 0,
      "Disminuye PCOS vs Female",
      "Aumenta PCOS vs Female"
    ),
    NA
  )
  
  
  Female_vs_Male <- ifelse(
    res_p$Female_vs_Male < 0.05 &
      res_p$Grupo < 0.05  &
      res_w_t$Female_vs_Male != 0,
    ifelse(
      res_w_t$Female_vs_Male > 0,
      "Disminuye Male vs Female",
      "Aumenta Male vs Female"
    ),
    NA
  )
  
  
    PCOS_vs_Male <- ifelse(
    res_p$PCOS_vs_Male < 0.05 &
      res_p$Grupo < 0.05  &
      res_w_t$PCOS_vs_Male != 0,
    ifelse(
      res_w_t$PCOS_vs_Male > 0,
      "Aumenta PCOS vs Male",
      "Disminuye PCOS vs Male"
    ),
    NA
  )
  
  
  
  NO_Female_vs_PCOS <-
    ifelse(
      res_p$NO_Female_vs_PCOS < 0.05 &
        res_p$Interaccion < 0.05  &
        res_w_t$NO_Female_vs_PCOS != 0,
      ifelse(
        res_w_t$NO_Female_vs_PCOS > 0,
        "Delgados Disminuye PCOS vs Female",
        "Delgados Aumenta PCOS vs Female"
      ),
      NA
    )
  
  
  NO_Female_vs_Male <-
    ifelse(
      res_p$NO_Female_vs_Male < 0.05 &
        res_p$Interaccion < 0.05  &
        res_w_t$NO_Female_vs_Male != 0,
      ifelse(
        res_w_t$NO_Female_vs_Male > 0,
        "Delgados Disminuye Male vs Female",
        "Delgados Aumenta Male vs Female"
      ),
      NA
    )
  
  
  NO_PCOS_Vs_Male <- ifelse(
    res_p$NO_PCOS_Vs_Male < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$NO_PCOS_Vs_Male != 0,
    ifelse(
      res_w_t$NO_PCOS_Vs_Male > 0,
      "Delgados Aumenta PCOS vs Male",
      "Delgados  Disminuye PCOS vs Male"
    ),
    NA
  )
  
  
  
  
  O_Female_vs_PCOS <- ifelse(
    res_p$O_Female_vs_PCOS < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$O_Female_vs_PCOS != 0,
    ifelse(
      res_w_t$O_Female_vs_PCOS > 0,
      "OBESOS Disminuye PCOS vs Female",
      "OBESOS Aumenta PCOS vs Female"
    ),
    NA
  )
  
  
  O_Female_vs_Male <- ifelse(
    res_p$O_Female_vs_Male < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$O_Female_vs_Male != 0,
    ifelse(
      res_w_t$O_Female_vs_Male > 0,
      "OBESOS Disminuye Male vs Female",
      "OBESOS Aumenta Male vs Female"
    ),
    NA
  )
  
  
  O_PCOS_Vs_Male <- ifelse(
    res_p$O_PCOS_Vs_Male < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$O_PCOS_Vs_Male != 0,
    ifelse(
      res_w_t$O_PCOS_Vs_Male > 0,
      "OBESOS Aumenta PCOS vs Male",
      "OBESOS  Disminuye PCOS vs Male"
    ),
    NA
  )
  
  
  PCOS <- ifelse(
    res_p$PCOS < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$PCOS != 0,
    ifelse(
      res_w_t$PCOS > 0,
      "Disminuye en OBESOS PCOS",
      "Aumenta en OBESOS PCOS"
    ),
    NA
  )
  
  
  
  Female <- ifelse(
    res_p$Female < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$Female != 0,
    ifelse(
      res_w_t$Female > 0,
      "Disminuye en OBESOS Female",
      "Aumenta en OBESOS Female"
    ),
    NA
  )
  
  
  
  Male <- ifelse(
    res_p$Male < 0.05 &
      res_p$Interaccion < 0.05  &
      res_w_t$Male != 0,
    ifelse(
      res_w_t$Male > 0,
      "Disminuye en OBESOS Male",
      "Aumenta en OBESOS Male"
    ),
    NA
  )
  
  
  
  Obesidad <- ifelse(
    res_p$Obesidad < 0.05 &
      res_w_t$NO_vs_O != 0,
    ifelse(res_w_t$NO_vs_O > 0, "Aumenta en OBESOS", "Disminuye en OBESOS"),
    NA
  )
  
  
  Grupo <-
    ifelse(res_p$Grupo < 0.05, "Hay diferencias entre grupos", NA)
  
  Interaccion <-
    ifelse(res_p$Interaccion < 0.05,
           "Hay diferencias entre grupo y obesidad",
           NA)
  
  res_interp <- as.data.frame(bind_cols(Grupo=Grupo,
                          Obesidad=Obesidad,
                          Interaccion=Interaccion,
                          Female_vs_PCOS=Female_vs_PCOS,
                          Female_vs_Male=Female_vs_Male,
                          PCOS_vs_Male=PCOS_vs_Male,
                          NO_Female_vs_PCOS=NO_Female_vs_PCOS,
                          NO_Female_vs_Male=NO_Female_vs_Male,
                          NO_PCOS_Vs_Male=NO_PCOS_Vs_Male,
                          O_Female_vs_PCOS=O_Female_vs_PCOS,
                          O_Female_vs_Male=O_Female_vs_Male,
                          O_PCOS_Vs_Male=O_PCOS_Vs_Male,
                          PCOS=PCOS,
                          Female=Female,
                          Male=Male))
  rownames(res_interp) <- rownames(res_p)
  res_interp <- (apply(res_interp,2,as.character))
  retornos <- vector("character",length=nrow(res_interp))
  for(row in 1:nrow(res_interp)){
    retornos[row] <- concatenando_grupo(res_interp[row,])
    
  }
  
  retornos.df <- data.frame(retornos=retornos)
  rownames(retornos.df) <- rownames(res_p)

  retorno <- list(bloque = X,
                  p_valores = res_p,
                  interpretacion = retornos.df)
  
  
  return(retorno)
}
