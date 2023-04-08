
library(reshape2)
library(dplyr)
library(ggplot2)
library(limma)
# library(MetaboAnalystR)
datos <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_none_.rds")

sexo <- datos$comunes$variables_in_bacteria$SEX
obesidad <- datos$comunes$variables_in_bacteria$OBESE
grupo <- datos$comunes$variables_in_bacteria$GROUP

comunes <- scale(t(datos$comunes$metaboloma))

# comunes <- comunes-apply(comunes,2,median)
# # transformation()
# # transformation()
# comunes <- (comunes / rowSums(comunes))
# comunes <- -t(apply(comunes, 1, log2))



# 
# 
comunes.plot <- as.data.frame((comunes))
comunes.plot$sexo <- sexo
comunes.plot$grupo <- grupo
comunes.plot$obesidad <- obesidad

comunes.plot.melt <- melt(comunes.plot)


ggplot(comunes.plot.melt,aes(value,color=sexo))+geom_density()


pcx <- prcomp((comunes))
plotear <- bind_cols(pcx$x,obesidad=obesidad)

ggplot(plotear,aes(PC1,PC2,color=obesidad,shape=grupo))+geom_point(size=2)



# 
# 
# 
# 
# 
# 
# 
# 
