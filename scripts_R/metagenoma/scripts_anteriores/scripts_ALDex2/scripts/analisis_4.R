

library(ALDEx2)
library(dplyr)
library(ggplot2)
library(zCompositions)

source("./12-2-22/scripts_analysis_reunion/scripts_funciones/otras_funciones_utiles.R")
funttest <- function(x,categotica) apply(x, 2, function(x) t.test(x ~ categotica)$p.value)
funwilcos <- function(x,categotica) apply(x, 2, function(x) wilcox.test(x ~ categotica)$p.value)
funaov <- function(x,categotica) apply(x, 2, function(x) summary(aov(x ~ categotica))[[1]]$`Pr(>F)`[1])
funkw <- function(x,categotica) apply(x, 2, function(x) kruskal.test(x ~ categotica)$p.value)

fun_summary <- function(datos,taxon){
  
  sexo <- datos$comunes$variables_in_bacteria$SEX
  obesidad <- datos$comunes$variables_in_bacteria$OBESE
  grupo <- datos$comunes$variables_in_bacteria$GROUP
  
  sexottest <- colMeans(do.call(rbind,lapply(taxon, function(y) funttest(t(y),sexo))))
  sexowilox <- colMeans(do.call(rbind,lapply(taxon, function(y) funwilcos(t(y),sexo))))
  obesidadttest <- colMeans(do.call(rbind,lapply(taxon, function(y) funttest(t(y),obesidad))))
  obesidadwilox <- colMeans(do.call(rbind,lapply(taxon, function(y) funwilcos(t(y),obesidad))))
  grupoparam <- colMeans(do.call(rbind,lapply(taxon, function(y) funaov(t(y),grupo))))
  gruponoparam <- colMeans(do.call(rbind,lapply(taxon, function(y) funkw(t(y),grupo))))
  
  
  res_sin_scale <- data.frame(sexott = sexottest,sexow = sexowilox,
                              obesidadtt = obesidadttest,obesidadwilox=obesidadwilox,
                              grupoaov = grupoparam,grupokw = gruponoparam)
  
  return(res_sin_scale)
}

convert_aldex <- function(X){
  
  mc.all <- getMonteCarloInstances(X)
  n <- length(mc.all)
  p <- nrow(mc.all[[1]])
  aux <- do.call(rbind,mc.all)
  my_list <- split(aux,rep(1:ncol(aux),each=nrow(aux)))
  my_list <- lapply(my_list,function(x) matrix(x,nrow=p,ncol=n))
  return(my_list)
  
}

datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")

genero <- (datos$comunes$microbiota$genero.abs)
filo <- (datos$comunes$microbiota$filo.abs)

grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE
sexo <- datos$comunes$variables_in_bacteria$SEX


#### por obesidad

taxon <- filo
iteraciones=1000
condiciones <- c("Obese","No Obese")
lista_obesidad <- vector("list",length=2)
for(i in 1:2){
  
  condicion <- condiciones[i]
  print(condicion)
  grupo.obesidad <- grupo[obesidad==condicion]
  grupo.condicion <- as.factor(as.character(grupo[obesidad==condicion & grupo!="Male"]))
  control.condicion <- as.factor(as.character(grupo[obesidad==condicion & grupo !="PCOS"]))
  males.condicion <- as.factor(as.character(grupo[obesidad==condicion & grupo!="Female"]))
  taxon.grupo <- taxon[obesidad==condicion & grupo!="Male",]
  taxon.males <- taxon[obesidad==condicion & grupo!="Female",]
  taxon.control <- taxon[obesidad==condicion & grupo!="PCOS",]
  
  condicion_grupo <- aldex(t(taxon.grupo),conditions = grupo.condicion,mc.samples = iteraciones)
  bacterias_grupo <- rownames(condicion_grupo[condicion_grupo$wi.ep<0.05 | condicion_grupo$we.ep<0.05,])
  if(length(bacterias_grupo)==0){
    resultado.grupo <- NA
  }else{
    taxon.grupo.res <- taxon.grupo[,bacterias_grupo]
    taxon.grupo.res$condicion <- grupo.condicion
    resultado.grupo <- t(taxon.grupo.res %>% 
                           group_by(condicion) %>% 
                           summarise(across(everything(), sum)))
    colnames(resultado.grupo) <- c("condicion",bacterias_grupo)
    resultado.grupo <- as.data.frame(resultado.grupo)
    resultado.grupo$pvalor <- condicion_grupo[rownames(condicion_grupo) %in% bacterias_grupo,c("wi.ep","we.ep","wi.eBH","we.eBH")]
    resultado.grupo$grupo.interes <- paste(condicion,"Female.vs.PCOS")
    
  }
    
  condicion_males <- aldex(t(taxon.males),conditions = males.condicion,mc.samples = iteraciones)
  bacterias_males <- rownames(condicion_males[condicion_males$wi.ep<0.05 | condicion_males$we.ep<0.05,])
  if(length(bacterias_males)==0){
    resultado.males <- NA
  }else{
    taxon.males.res <- taxon.males[,bacterias_males]
    taxon.males.res$condicion <- males.condicion
    resultado.males <- t(taxon.males.res %>% 
                           group_by(condicion) %>% 
                           summarise(across(everything(), sum)))
    colnames(resultado.males) <- c("condicion",bacterias_males)
    resultado.males <- as.data.frame(resultado.males)
    resultado.males$pvalor <- condicion_males[rownames(condicion_males) %in% bacterias_males,c("wi.ep","we.ep","wi.eBH","we.eBH")]
    resultado.males$grupo.interes <- paste(condicion,"PCOS.vs.Male")
    
  }
    
  condicion_control <- aldex(t(taxon.control),conditions = control.condicion,mc.samples = iteraciones)
  bacterias_control <- rownames(condicion_control[condicion_control$wi.ep<0.05 | condicion_control$we.ep<0.05,])
  if(length(bacterias_control)==0){
    resultado.control <- NA
  }else{
    
  
  taxon.control.res <- as.data.frame(taxon.control[,bacterias_control])
  taxon.control.res$condicion <- control.condicion
  resultado.control <- (taxon.control.res %>% 
                           group_by(condicion) %>% 
                           summarise(across(everything(), sum)))
  colnames(resultado.control) <- c("condicion",bacterias_control)
  resultado.control <-(resultado.control)
  resultado.control$pvalor <- condicion_control[rownames(condicion_control) %in% bacterias_control,c("wi.ep","we.ep","wi.eBH","we.eBH")]
  resultado.control$grupo.interes <- paste(condicion,"Male.vs.Female")
  }
  resultados <- list(females=resultado.grupo,
                     males=resultado.males,
                     control=resultado.control)
  lista_obesidad[[i]] <- resultados
  
  
  
}
names(lista_obesidad) <- condiciones
##### por grupo
condiciones <- c("PCOS","Male","Female")
lista_grupo <- vector("list",length=3)
for( i in 1:3){
  condicion2 <- condiciones[i]
  print(condicion2)
  grupo.condicion <- as.factor(as.character(obesidad[grupo==condicion2]))
  taxon.grupo <- taxon[grupo==condicion2 ,]
  
  condicion_grupo <- aldex(t(taxon.grupo),conditions = grupo.condicion,mc.samples = iteraciones)
  bacterias_grupo <- rownames(condicion_grupo[condicion_grupo$wi.ep<0.05 | condicion_grupo$we.ep<0.05,])
  if(length(bacterias_grupo)==0){
    resultado.grupo <- NA
  }else{
    

  taxon.grupo.res <- as.data.frame(taxon.grupo[,bacterias_grupo])
  taxon.grupo.res$condicion <- grupo.condicion
  resultado.grupo <- (taxon.grupo.res %>% 
                         group_by(condicion) %>% 
                         summarise(across(everything(), sum)))
  colnames(resultado.grupo) <- c("condicion",bacterias_grupo)
  resultado.grupo <- (resultado.grupo)
  resultado.grupo$pvalor <- condicion_grupo[rownames(condicion_grupo) %in% bacterias_grupo,c("wi.ep","we.ep","wi.eBH","we.eBH")]
  resultado.grupo$grupo.interes <- paste(condicion2,"Obese vs NO Obese")
  }
  lista_grupo[[i]] <- resultado.grupo
}
names(lista_grupo) <- condiciones

resultados.sexo <- aldex(t(taxon),conditions = sexo,mc.samples = iteraciones)

bacterias_sexo <- rownames(resultados.sexo[resultados.sexo$wi.ep<0.05 | resultados.sexo$we.ep<0.05,])
if(length(bacterias_sexo)==0){
  resultado.sexo <- NA
}else{
  
taxon.sexo.res <- taxon[,bacterias_sexo]
taxon.sexo.res$condicion <- sexo
resultado.sexo <- t(taxon.sexo.res %>% 
                      group_by(condicion) %>% 
                      summarise(across(everything(), sum)))
colnames(resultado.sexo) <- resultado.sexo[1,]
resultado.sexo<-resultado.sexo[-1,]
resultado.sexo <- as.data.frame(resultado.sexo)
resultado.sexo$pvalor <- resultados.sexo[rownames(resultados.sexo) %in% bacterias_sexo,c("wi.ep","we.ep" ,"wi.eBH","we.eBH")]
resultado.sexo$grupo.interes <- "Male vs Female"

}


resultados.obesidad <- aldex(t(taxon),conditions = obesidad,mc.samples = iteraciones)

bacterias_obesidad <- rownames(resultados.obesidad[resultados.obesidad$wi.ep<0.05 | resultados.obesidad$we.ep<0.05,])
if(length(bacterias_obesidad)==0){
  resultado.obesidad<- NA
}else{
taxon.obesidad.res <- taxon[,bacterias_obesidad]
taxon.obesidad.res$condicion <- obesidad
resultado.obesidad <- t(taxon.obesidad.res %>% 
                          group_by(condicion) %>% 
                          summarise(across(everything(), sum)))
colnames(resultado.obesidad) <- resultado.obesidad[1,]
resultado.obesidad<-resultado.obesidad[-1,]
resultado.obesidad <- as.data.frame(resultado.obesidad)
resultado.obesidad$pvalor <- resultados.obesidad[rownames(resultados.obesidad) %in% bacterias_obesidad,c("wi.ep","we.ep","wi.eBH","we.eBH")]
resultado.obesidad$grupo.interes <- "Obese vs NO Obese"
}
grupo_function <- function(grupo.interes){
  taxon.grupo.interes <- taxon[grupo!=grupo.interes,]
  condiciones <- as.factor(as.character(grupo[grupo!=grupo.interes]))
  resultados.grupo.interes <- aldex(t(taxon.grupo.interes),conditions = condiciones,mc.samples = iteraciones,test = "t")
  
  bacterias_grupo <- rownames(resultados.grupo.interes[resultados.grupo.interes$wi.ep<0.05 | resultados.grupo.interes$we.ep<0.05,])
  if(length(bacterias_grupo)==0){
    resultado.grupo.interes<- NA
  }else{
   taxon.grupo.res <- as.data.frame(taxon.grupo.interes[,bacterias_grupo])
  taxon.grupo.res$condicion <- condiciones
  resultado.grupo.interes <- as.data.frame(t(taxon.grupo.res %>% 
                                               group_by(condicion) %>% 
                                               summarise(across(everything(), sum))))
  colnames(resultado.grupo.interes) <- resultado.grupo.interes[1,]
  resultado.grupo.interes<-resultado.grupo.interes[-1,]
  resultado.grupo.interes <- as.data.frame(resultado.grupo.interes)
  resultado.grupo.interes$pvalor <- resultados.grupo.interes[rownames(resultados.grupo.interes) %in% bacterias_grupo,c("wi.ep","we.ep","wi.eBH","we.eBH")]
  resultado.grupo.interes$grupo.interes <- paste(unique(as.character(grupo[grupo!=grupo.interes])),collapse = " ")
  rownames(resultado.grupo.interes) <- bacterias_grupo
  }
  return(resultado.grupo.interes)
}

PCOS.vs.Female <- grupo_function("Male")
Male.vs.Female <- grupo_function("PCOS")
PCOS.vs.Male <- grupo_function("Female")


resultados.grupo <- aldex(t(taxon),conditions = grupo,mc.samples = iteraciones,test = "kw")

bacterias_grupo <- rownames(resultados.grupo[resultados.grupo$kw.ep<0.05 | resultados.grupo$glm.ep<0.05,])
if(length(bacterias_grupo)==0){
  resultado.grupo <- NA
}else{
  taxon.grupo.res <- taxon[,bacterias_grupo]
  taxon.grupo.res$condicion <- grupo
  resultado.grupo <- t(taxon.grupo.res %>% 
                         group_by(condicion) %>% 
                         summarise(across(everything(), sum)))
  colnames(resultado.grupo) <- resultado.grupo[1,]
  resultado.grupo<-resultado.grupo[-1,]
  resultado.grupo <- as.data.frame(resultado.grupo)
  resultado.grupo$pvalor <- resultados.grupo[rownames(resultados.grupo) %in% bacterias_grupo,c("kw.ep","glm.ep","kw.eBH","glm.eBH")]
  resultado.grupo$grupo.interes <-"grupo"
  
}


directorio <- "./12-2-22/scripts_analysis_reunion/scripts_metagenoma_novoom/resultados_filo"
dir.create(directorio)
write.csv(resultado.sexo,file=file.path(directorio,"sexo.csv"))
write.csv(resultado.obesidad,file=file.path(directorio,"obesidad.csv"))
write.csv(resultado.grupo,file=file.path(directorio,"diferencias_grupo.csv"))
write.csv(lista_obesidad$Obese$females,file=file.path(directorio,"obesidad_mujeres_pcos.csv"))
write.csv(lista_obesidad$Obese$males,file=file.path(directorio,"obesidad_hombresPCOS.csv"))
write.csv(lista_obesidad$Obese$control,file=file.path(directorio,"obesidad_sexo_control.csv"))

write.csv(lista_obesidad$`No Obese`$females,file = file.path(directorio,"NoObesidad_mujeres.csv"))
write.csv(lista_obesidad$`No Obese`$males,file = file.path(directorio,"Noobesidad_hombresPCOS.csv"))
write.csv(lista_obesidad$`No Obese`$control,file = file.path(directorio,"NoObesidad_sexo_control.csv"))
write.csv(lista_grupo$PCOS,file=file.path(directorio,"Obesidad_PCOS.csv"))
write.csv(lista_grupo$Female,file = file.path(directorio,"obesidad_mujeres.csv"))
write.csv(lista_grupo$Male,file=file.path(directorio,"obesidad_hombres.csv"))
write.csv(PCOS.vs.Female,file=file.path(directorio,"PCOS.vs.Female.csv"))
write.csv(PCOS.vs.Male,file=file.path(directorio,"PCOS.vs.Male.csv"))
write.csv(Male.vs.Female,file=file.path(directorio,"Male.vs.Female.csv"))

# # 
