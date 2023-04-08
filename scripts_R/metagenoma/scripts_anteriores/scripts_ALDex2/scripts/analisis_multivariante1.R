
library(vegan)
fichero <- "../../datos/preprocesado_05_02_23/novoom.rds"
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/otras_funciones_utiles.R")
datos <- readRDS(fichero)

sexo <- datos$comunes$variables_in_bacteria$SEX
obesidad <- datos$comunes$variables_in_bacteria$OBESE
grupo <- datos$comunes$variables_in_bacteria$GROUP

convert_aldex <- function(X){
  
  mc.all <- getMonteCarloInstances(X)
  n <- length(mc.all)
  p <- nrow(mc.all[[1]])
  aux <- do.call(rbind,mc.all)
  my_list <- split(aux,rep(1:ncol(aux),each=nrow(aux)))
  my_list <- lapply(my_list,function(x) matrix(x,nrow=p,ncol=n))
  return(my_list)
  
}


filo <- datos$comunes$microbiota$filo
genero <- datos$comunes$microbiota$genero
filo.convert <- scale(t(mean_aldex(filo)))
genero.convert <- scale(t(mean_aldex(genero)))

multi_fun <- function(x){
  x <- scale((x))
  pmanova <- adonis2(x ~ grupo*obesidad,method="euclidean")
  controles <- as.factor(as.character(grupo[grupo!="PCOS"]))
  herma <- as.factor(as.character(grupo[grupo!="Female"]))
  mujeres <- as.factor(as.character(grupo[grupo!="Male"]))
  controles_t <- adonis2(x[grupo!="PCOS",]~controles,method="euclidean")$`Pr(>F)`[1]
  herma_t <- adonis2(x[grupo!="Female",]~herma,method="euclidean")$`Pr(>F)`[1]
  mujeres_t <- adonis2(x[grupo!="Male",]~mujeres,method="euclidean")$`Pr(>F)`[1]
  sexo <- adonis2(x~sexo,method="euclidean")$`Pr(>F)`[1]
  if(pmanova$`Pr(>F)`[1]<0.05){
    p.grupo <- pmanova$`Pr(>F)`[1]
    return(c(grupo=p.grupo,
             mujeres=mujeres_t,
             herma=herma_t,
             controles=controles_t,
             sexo=sexo,
             obesidad=pmanova$`Pr(>F)`[2],
             interaccion =pmanova$`Pr(>F)`[3]))
    
  }else{
    p.grupo <- pmanova$`Pr(>F)`[1]
    p.obesidad=pmanova$`Pr(>F)`[2]
    p.interaccion=pmanova$`Pr(>F)`[3]
    return(c(grupo=p.grupo,
           obesidad=p.obesidad,
          interaccion=p.interaccion))
  }
  
  
}


res.genero <- multi_fun(genero.convert)

adj.grupo <- p.adjust(res.genero[2:4])

res.filo <- multi_fun(filo.convert)





