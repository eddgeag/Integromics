
library(vegan)
library(car)
library(Hotelling)
library(psych)
vanValen.test <- function(X, categorica) {
  X.standard <- scale(X)
  X.standard <- data.frame(X.standard, categorica = categorica)
  X.dat <-
    split(X.standard[, which(unlist(lapply(X.standard, is.numeric)) == T)],
          X.standard[, which(unlist(lapply(X.standard, is.factor)) ==
                               T)])
  X.dat <- lapply(X.dat, as.matrix)
  medianas <- lapply(X.dat, function(x)
    apply(x, 2, median))
  abs.dev.factor <- vector("list", length = length(X.dat))
  for (j in 1:length(levels(categorica))) {
    y <- X.dat[[j]]
    mediana_x <- medianas[[j]]
    abs.dev.factor[[j]] <- abs(sweep(y, 2, mediana_x))
    
  }
  
  abs.dev.all <- Reduce("rbind", abs.dev.factor)
  abs.dev.all <- data.frame(categorica, abs.dev.all)
  norma.euclidea <- function(x)
    sqrt(sum(x ^ 2))
  distancias <-
    data.frame(categorica, dist = apply(abs.dev.all[, -1], 1, norma.euclidea))
  norm.test <-
    shapiro.test(residuals(lm(distancias[, 2] ~ distancias[, 1])))$p.value
  vari.test <-
    leveneTest(distancias[, 2], group = distancias[, 1])$`Pr(>F)`[1]
  
  if (length(levels(categorica)) == 2) {
    if (norm.test > 0.05 && vari.test > 0.05) {
      p.valor <-
        t.test(distancias[, 2] ~ distancias[, 1], var.equal = T)$p.value
      
    } else if (norm.test > 0.05 && vari.test < 0.05) {
      p.valor <-
        t.test(distancias[, 2] ~ distancias[, 1], var.equal = F)$p.value
      
    } else if (norm.test < 0.05) {
      p.valor <-
        wilcox.test(distancias[, 2] ~ distancias[, 1], exact = F)$p.value
      
    }
    
  } else if (length(levels(categorica)) > 2) {
    pregunta <- as.data.frame(table(categorica))$Freq
    mod <- lm(distancias[, 2] ~ distancias[, 1])
    if (all(pregunta)) {
      if (norm.test > 0.05 && vari.test > 0.05) {
        p.valor <- Anova(mod)$`Pr(>F)`[1]
        
      } else if (norm.test > 0.05 && vari.test < 0.05) {
        p.valor <- Anova(mod, white.adjust = T)$`Pr(>F)`[1]
        
      } else if (norm.test < 0.05) {
        p.valor <-
          kruskal.test(distancias[, 2] ~ distancias[, 1])$p.value
        
      }
    }
    
  }
  return(p.valor)
  
}



aux_mardia <- function(X){
  X <- as.matrix(X)
  if(length(try(mardia(X,plot=F),T))==1){
    pvals <- c(NA,NA)
  }else{
    pvals.mardia <- mardia(X,plot=F)
    pvals <- c(pvals.mardia$p.kurt,pvals.mardia$p.skew)
  }
  return(pvals)
}


normalidad_multivariante <- function(X,var_i){

  X.split <- split(X,var_i)
  X.split <-lapply(X.split, as.data.frame)
  pvals <- lapply(X.split, function(x) aux_mardia(x[,-ncol(x)]))
  
  isna_2 <- function(X){
    if(all(is.na(X))){
      return(0)
    }else if(all(is.na(X))==F && anyNA(X)==T){
      
      return(1)
      
    }else if(all(!is.na(X)==T )&& all(X>0.05)){
      
      return(2)
      
    }else if(all(!is.na(X)==T) && any(X<0.05) ){
      return(3)
    }else if(all(is.na(X))==T){
      return(4)
    }
      
  }
  
 retorno <- unlist(lapply(pvals, isna_2))
  
  return(retorno)
  
  
  
}

normalidad_multivariante2 <- function(X_,var_i){

  X.split <- split(X_,var_i)
  X.split <-lapply(X.split, as.data.frame)
  pvals <- lapply(X.split, function(x) aux_mardia(x[,-ncol(x)]))

 
  
  retorno <- unlist(pvals)
  return(retorno)
  
  
  
}


funcion_multivariante2 <-
  function(X_2,
           categorica,
           categoricas,
           nombre_analisis) {
    # X_2 <- analisis_2
    # categorica <- obesidad
    # categoricas <- list(grupo,obesidad)
    # nombre_analisis<-"preuba"
    n <- nrow(X_2)
    p <- ncol(X_2)
    X_2[is.na(X_2)] <- 0
    categorica1 <- categoricas[[1]]
    categorica2 <- categoricas[[2]]
    interaccion <- interaction(categorica1,categorica2)
    norm_multi_i <- normalidad_multivariante2(X_2,categorica)
    norm_multi_interaccion <- normalidad_multivariante2(X_2,interaccion)

    if(all(!is.na(norm_multi_i)) && all(norm_multi_i>0.05)){
      
      if(length(levels(categorica))==2){
        
        p.vanvalen <- vanValen.test(X_2, categorica)
        
        if (p.vanvalen > 0.05 ) {
          X.split <- split(as.data.frame(X_2), categorica)
          pval <-
            hotelling.test(as.matrix(X.split[[1]]),
                           as.matrix(X.split[[2]]),
                           var.equal = T)$pval
          retorno <- pval
        }else{
          X.split <- split(as.data.frame(X_2), categorica)
          
          pval <-
            hotelling.test(as.matrix(X.split[[1]]),
                           as.matrix(X.split[[2]]),
                           var.equal = F,shrinkage = T,perm = F)$pval
          retorno <- pval
          
        }
        
      }else if(length(levels(categorica))>2){
        p.vanvalen <- vanValen.test(X_2, categorica)
     
        if (p.vanvalen > 0.05) {
          pval <- summary(manova(as.matrix(X_2)~categorica))$stats[1, 6]
          
          retorno <- pval
          
        } else{
          pval <- adonis2(as.matrix(X_2) ~ categorica, method = "euclidean")$`Pr(>F)`[1]
          retorno <- pval
          

        
      }
      
    }
    

    }else{
      
      retorno <- c(NA)
    
      
    }
    
    if(all(!is.na(norm_multi_interaccion)) && all(norm_multi_interaccion>0.05)){
      
      
      p.vanvalen <- vanValen.test(X_2, interaccion)
      if (p.vanvalen > 0.05) {
        X_2 <- as.matrix(X_2)
        
        pval <- summary(manova(X_2~categorica1*categorica2))$stats[1, 6]
        
        retorno2 <- pval
        
      } else{
        pval <- adonis2(X_2~categorica1*categorica2, method = "euclidean")$`Pr(>F)`[3]
        retorno2 <- pval
        
        
        
      }
      
    }else{
      retorno2 <- NA
    }
    
    retornos <- list(retorno, retorno2)
    names(retornos) <-
      c(paste0(nombre_analisis, ":", paste(levels(categorica), collapse = "_")),
        paste0(nombre_analisis, "_interaccion"))
    
    return(retornos)
}

### esta funcion me mira la normalidad en los residuos
funcion_multivariante <-
  function(X,
           categorica,
           categoricas,
           nombre_analisis) {
    n <- nrow(X)
    p <- ncol(X)
    X[is.na(X)] <- 0
    X <- as.matrix(X)
    categorica1 <- categoricas[[1]]
    categorica2 <- categoricas[[2]]


    
    if (length(levels(categorica)) == 2) {
      ## primero checamos las suposiciones
      modelo <- manova(X ~ categorica)
      residuos <- residuals(modelo)
      if (qr(residuos)$rank < p |
          length(class(try(mardia(X, plot = F), T)
          )) == 1)  {
        p.mardia <- c(NA, NA)
      } else{
        tmp <- mardia(X, plot = F)
        kurtosis <- tmp$p.kurt
        skewness <- tmp$p.skew
        p.mardia <- c(kurtosis = kurtosis, skewness = skewness)
      }
      if (anyNA(p.mardia)) {
        retorno <- NA
      } else{
        if (all(!is.na(p.mardia))) {
          p.vanvalen <- vanValen.test(X, categorica)
          if (p.vanvalen > 0.05 && all(p.mardia>0.05)) {
            X.split <- split(as.data.frame(X), categorica)
            pval <-
              hotelling.test(as.matrix(X.split[[1]]),
                             as.matrix(X.split[[2]]),
                             var.equal = T)$pval
            retorno <- pval
          } else{
            X.split <- split(as.data.frame(X), categorica)
            pval <-
              hotelling.test(
                as.matrix(X.split[[1]]),
                as.matrix(X.split[[2]]),
                var.equal = F,
                shrinkage = T
              )$pval
            retorno <- pval
          }
          
        } else{
          retorno <- NA
        }
        
      }
      
      
    } else if (length(levels(categorica)) > 2) {
      modelo2_ <- manova(X ~ categorica)
      residuos <- residuals(modelo2_)
      if (qr(residuos)$rank < p |
          length(class(try(mardia(X, plot = F), T)
          )) == 1) {
        p.mardia <- c(NA, NA)
      } else{
        tmp <- mardia(X, plot = F)
        kurtosis <- tmp$p.kurt
        skewness <- tmp$p.skew
        p.mardia <- c(kurtosis = kurtosis, skewness = skewness)
      }
      if (anyNA(p.mardia)) {
        retorno <- NA
      } else{
        if (all(p.mardia > 0.05 ) && all(p.mardia>0.05)) {
          p.vanvalen <- vanValen.test(X, categorica)
          if (p.vanvalen > 0.05) {
            pval <- summary(modelo2_)$stats[1, 6]
            
            retorno <- pval
            
          } else{
            pval <- adonis2(X ~ categorica, method = "euclidean")$`Pr(>F)`[1]
            retorno <- pval
            
          }
          
        } else{
          retorno <- NA
        }
        
      }
      
    }
    
    modelo2 <- manova(X ~ categorica1 * categorica2)
    residuos <- residuals(modelo2)
    if (qr(residuos)$rank < p |
        length(class(try(mardia(X, plot = F), T)
        )) == 1) {
      p.mardia <- c(NA, NA)
    } else{
      tmp <- mardia(X, plot = F)
      kurtosis <- tmp$p.kurt
      skewness <- tmp$p.skew
      p.mardia <- c(kurtosis = kurtosis, skewness = skewness)
    }
    if (anyNA(p.mardia)) {
      retorno2 <- NA
    } else{
      if (all(p.mardia > 0.05 )){
        p.vanvalen <-
          vanValen.test(X, interaction(categorica1, categorica2))
        if (p.vanvalen > 0.05) {
          pval <- summary(modelo2)$stats[3, 6]
          
          retorno2 <- pval
        } else{
          pval <-
            adonis2(X ~ categorica1 * categorica2, method = "euclidean")$`Pr(>F)`[3]
          retorno2 <- pval
          
        }
        
      } else{
        retorno2 <- NA
      }
      
    }
    
    retornos <- list(retorno, retorno2)
    names(retornos) <-
      c(paste0(nombre_analisis, ":", paste(levels(categorica), collapse = "_")),
        paste0(nombre_analisis, "_interaccion"))
    return(retornos)
  }


funcion_multivariante_perm <-
  function(X,
           categorica,
           categoricas,
           nombre_analisis) {
    n <- nrow(X)
    p <- ncol(X)
    X[is.na(X)] <- 0
    X <- as.matrix(X)
    categorica1 <- categoricas[[1]]
    categorica2 <- categoricas[[2]]
    
    retorno <-
      adonis2(X ~ categorica, method = "euclidean", permutations = 10000)$`Pr(>F)`[1]
    retorno2 <-
      adonis2(X ~ categorica1 * categorica2,
              method = "euclidean",
              permutations = 10000)$`Pr(>F)`[1]
    retornos <- list(retorno, retorno2)
    names(retornos) <-
      c(paste0(nombre_analisis, ":", paste(levels(categorica), collapse = "_")),
        paste0(nombre_analisis, "_interaccion"))
    return(retornos)
  }