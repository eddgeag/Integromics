library(reshape2)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(ggrepel)
library(ALDEx2)
library(vegan)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source(
  "./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R"
)
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")
datos <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")

directorio.totales <-
  "./scripts_R/CLINICOS/resultados_totales_bivariante"

if (!dir.exists(directorio.totales)) {
  dir.create(directorio.totales)
}

directorio.comunes <-
  "./scripts_R/CLINICOS/resultados_comunes_bivariante"

if (!dir.exists(directorio.comunes)) {
  dir.create(directorio.comunes)
}


comunes <- F

if (comunes) {
  directorio <- directorio.comunes
  
  obesidad <- datos$comunes$variables_in_bacteria$OBESE
  grupo <- datos$comunes$variables_in_bacteria$GROUP
  sexo <- datos$comunes$variables_in_bacteria$SEX
  X <- scale(log2(as.matrix(datos$comunes$clinicos)))
  
} else{
  directorio <- directorio.totales
  
  obesidad <- datos$totales$general_data$OBESE
  grupo <- datos$totales$general_data$GROUP
  sexo <- datos$totales$general_data$SEX
  X <- scale(log2(as.matrix(datos$totales$clinicos)))
  
}



objcor <- cor.test(X,method = "spearman")

# 
# cordata <- objcor$r
# pdata <- objcor$p
# cordata[c(ncol(cordata):(ncol(cordata) - length(top.vars.interes) + 1), ncol(cordata)), ] <-
#   NA
# cordata[, 1:length(top.vars.interes)] <- NA
# cordata <-
#   matrix(cordata[!is.na(cordata)],
#          ncol = length(top.vars.interes),
#          byrow = T)
# rownames(cordata) <- top.vars.interes
# colnames(cordata) <- paste0("PC", 1:length(top.vars.interes))
# pdata <- objcor$p
# pdata[c(ncol(pdata):(ncol(pdata) - length(top.vars.interes) + 1), ncol(pdata)), ] <-
#   NA
# pdata[, 1:length(top.vars.interes)] <- NA
# pdata <-
#   matrix(pdata[!is.na(pdata)], ncol = length(top.vars.interes), byrow = T)
# rownames(pdata) <- top.vars.interes
# colnames(pdata) <- paste0("PC", 1:length(top.vars.interes))
# 
# if (any(pdata < 0.05)) {
#   p1 <-
#     ggcorrplot::ggcorrplot(
#       cordata,
#       p.mat = 1 - pdata,
#       hc.order = F,
#       sig.level = 0.95,
#       insig = "pch",
#       pch.cex = 3
#     )
#   
# } else{
#   p1 <-
#     ggcorrplot::ggcorrplot(
#       cordata,
#       p.mat = NULL,
#       hc.order = F,
#       sig.level = 0.95,
#       insig = "pch",
#       pch.cex = 3
#     )
#   
#   
# }
# p1
# ggsave(file.path(directorio, "correlacion_Vars_Pcs.jpg"), p1)
# 
# # Contributions of variables to PC1
# p1 <- fviz_contrib(pcx,
#                    choice = "var",
#                    axes = 1,
#                    top = 10)
# # Contributions of variables to PC2
# p2 <- fviz_contrib(pcx,
#                    choice = "var",
#                    axes = 2,
#                    top = 10)
# 
# p3 <- ggarrange(p1, p2)
# p3
# 
# ggsave(
#   filename = file.path(directorio, "graficos_contribuciones_PCA.jpg"),
#   plot = p3
# )
