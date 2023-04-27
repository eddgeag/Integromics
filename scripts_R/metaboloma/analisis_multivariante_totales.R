library(reshape2)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(ggrepel)
library(ALDEx2)
library(vegan)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")
datos <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")

directorio.totales <- "./scripts_R/metaboloma/resultados_totales_multivariante"

if(!dir.exists(directorio.totales)){
  dir.create(directorio.totales)
}

directorio.comunes <- "./scripts_R/metaboloma/resultados_comunes_multivariante"

if(!dir.exists(directorio.comunes)){
  dir.create(directorio.comunes)
}

transformacion <- function(X){
  return((t(log2(X/colSums(X)))))
}



directorio <- directorio.totales

obesidad <- datos$totales$general_data$OBESE
grupo <- datos$totales$general_data$GROUP
sexo <- datos$totales$general_data$SEX
metaboloma.total <- datos$totales$metaboloma


X <- transformacion(metaboloma.total)


columnas_X <- add_classificacion_metaboloma(data.frame(columnas=colnames(X)))
library(limma)
csa <- add_classificacion_metaboloma(as.data.frame(strsplit2(colnames(X),"_")))
csa$grupo_meta <- strsplit2(csa$grupo_meta,"_")[,1]
colnames(X) <- apply(csa,1,function(x) paste(x[1],x[3],sep = "_"))

X <- X[,order(csa$grupo_meta)]
resultado <-
  analyze_multivariate(X, grupo, obesidad, sexo = sexo)



X.df <- as.data.frame(as.matrix((dist(as.matrix(X)))))
X.df$obesidad <- obesidad
X.df$grupo <- grupo
X.df$sexo <- sexo
X.melt <- reshape2::melt(X.df)

all_grupos <-
  X.melt %>% group_by(obesidad, grupo) %>% dplyr::summarise(sum =mean(value))
all_grupos$sum <- 100 * all_grupos$sum / sum(all_grupos$sum)

all_obesidad <-
  X.melt %>% group_by(obesidad) %>% dplyr::summarise(sum = mean(value))
all_obesidad$sum <- 100 * all_obesidad$sum / sum(all_obesidad$sum)

all_grupo <-
  X.melt %>% group_by(grupo) %>% dplyr::summarise(sum = mean(value))
all_grupo$sum <- 100 * all_grupo$sum / sum(all_grupo$sum)

resultados_prop <- bind_rows(all_obesidad, all_grupo, all_grupos)



write.csv(resultados_prop,
          file.path(directorio, "resultados_proporciones.csv"))

write.csv(resultado,
          file.path(directorio, "resultados_multivariantes.csv"))


barras_all <-
  X.melt %>% group_by(obesidad, grupo) %>% dplyr::summarise(suma = mean(value))

p1 <-
  ggplot(data = barras_all, aes(x = grupo, y = suma, fill = obesidad)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = c('steelblue', 'red')) + xlab("") + ylab("") 
barras_obesidad <-
  X.melt %>% group_by(obesidad) %>% dplyr::summarise(suma = mean(value))

p2 <-
  ggplot(data = barras_obesidad, aes(
    x = obesidad,
    y = suma,
    fill = "steelblue"
  )) + geom_bar(stat = "identity",
                position = "dodge") +
  scale_fill_manual(values = c('steelblue')) + xlab("") + ylab("%") + theme(legend.position = "none")

barras_grupo <-
  X.melt %>% group_by(grupo) %>% dplyr::summarise(suma =  mean(value))

p3 <-
  ggplot(data = barras_grupo, aes(
    x = grupo,
    y = suma,
    fill = "steelblue"
  )) + geom_bar(stat = "identity",
                position = "dodge") +
  scale_fill_manual(values = c('steelblue')) + xlab("") + ylab("%") + theme(legend.position = "none")


barras_sexo <-
  X.melt %>% group_by(sexo) %>% dplyr::summarise(suma = mean(value))

p4 <-
  ggplot(data = barras_sexo, aes(
    x = sexo,
    y = suma,
    fill = "steelblue"
  )) + geom_bar(stat = "identity",
                position = "dodge") +
  scale_fill_manual(values = c('steelblue')) + xlab("") + ylab("%") + theme(legend.position = "none")
p5 <- ggarrange(p1, p2, p3, p4)



p6 <- annotate_figure(p5, top = "genero")

p6

ggsave(filename = file.path(directorio, "graficos_proporciones.jpg"),
       plot = p6)

### Realizamos el análisis univariante y graficamos
res <- analyze_data(X, grupo, obesidad,correccion = 3)
medianas <- as.matrix(res$medianas)
colnames(medianas)[1] <- "Obesidad"
p.values <- as.matrix(res$p_valores)[,-c(1, 3)]

p.values.mask <-
  (apply(as.matrix(p.values), 1, function(x)
    any(x < 0.05)))


jpeg(file = file.path(directorio, "univariante.jpg"))
corrplot::corrplot(
  medianas,
  method = "color",
  is.corr = F,
  tl.cex = 1,
  number.cex = 0.6,
  p.mat = p.values,
  sig.level = 0.05,
  insig = "label_sig",col =  ifelse(medianas <= 0, "blue", "red"),
  pch.cex = 1,
  pch.col = "white"
)
medianas_ <- medianas
dev.off()
### HACEMOS PCA


pcx <- prcomp((X), scale. = F)


scores <- pcx$x
loadings <- pcx$rotation

res <- analyze_data(scores, grupo, obesidad)
medianas <- as.matrix(res$medianas)
colnames(medianas)[1] <- "Obesidad"
p.values <- as.matrix(res$p_valores)[,-c(1, 3)]

p.values.mask <-
  (apply(as.matrix(p.values), 1, function(x)
    any(x < 0.05)))

p.values <- p.values[p.values.mask,]
medianas <- medianas[p.values.mask,]
jpeg(file = file.path(directorio, "asociacion_PCA.jpg"))
corrplot::corrplot(
  medianas,
  method = "color",
  is.corr = F,
  tl.cex = 1,
  number.cex = 0.6,
  p.mat = p.values,
  sig.level = 0.05,
  insig = "label_sig",
  bg = ifelse(medianas >= 0, "blue", "red")
)

dev.off()

principal_components <- data.frame(
  PC1 =  scores[,1],
  PC2 = scores[, 2],
  grupo = grupo,
  obesidad = obesidad,
  sexo = sexo
)
varianzas  <-  round(100*pcx$sdev ^ 2 / sum(pcx$sdev ^ 2),2)
p1 <-
  ggplot(principal_components, aes(PC1, PC2, color = grupo)) + geom_point() +
  facet_wrap( ~ obesidad) + ylab(paste("PC2", varianzas[2], "%")) + xlab("")
p2 <-
  ggplot(principal_components, aes(PC1, PC2, color = obesidad)) + geom_point() +
  facet_wrap( ~ grupo) + xlab(paste("PC1", varianzas[1], "%")) + ylab(paste("PC2", varianzas[2], "%"))
p3 <- ggarrange(p1, p2, ncol = 1, nrow = 2)


p4 <- annotate_figure(p3, top = 'PCA')
p4

ggsave(file.path(directorio, "PCA_scores.jpg"), p4)

top <- 10

p1 <-
  fviz_pca_var(pcx, choice = "var", select.var = list(contrib = top))
p2 <- fviz_screeplot(pcx, addlabels = TRUE, ylim = c(0, 15))

p3 <- ggarrange(p1, p2)
p3
ggsave(file.path(directorio, "PCA_scree_loadings.jpg"), p3)

dds <- p1$data
top.vars.interes <- rownames(dds)
objcor <-
  psych::corr.p(cor(cbind(X[, top.vars.interes], scores[, 1:top]), method = "spearman"),
                n = length(top.vars.interes),
                adjust = "BH")
cordata <- objcor$r
pdata <- objcor$p
cordata[c(ncol(cordata):(ncol(cordata) - length(top.vars.interes) + 1), ncol(cordata)),] <-
  NA
cordata[, 1:length(top.vars.interes)] <- NA
cordata <-
  matrix(cordata[!is.na(cordata)],
         ncol = length(top.vars.interes),
         byrow = T)
rownames(cordata) <- top.vars.interes
colnames(cordata) <- paste0("PC", 1:length(top.vars.interes))
pdata <- objcor$p
pdata[c(ncol(pdata):(ncol(pdata) - length(top.vars.interes) + 1), ncol(pdata)),] <-
  NA
pdata[, 1:length(top.vars.interes)] <- NA
pdata <-
  matrix(pdata[!is.na(pdata)], ncol = length(top.vars.interes), byrow = T)
rownames(pdata) <- top.vars.interes
colnames(pdata) <- paste0("PC", 1:length(top.vars.interes))

if (any(pdata < 0.05)) {
  p1 <-
    ggcorrplot::ggcorrplot(
      cordata,
      p.mat = 1 - pdata,
      hc.order = F,
      sig.level = 0.95,
      insig = "pch",
      pch.cex = 3
    )
  
} else{
  p1 <-
    ggcorrplot::ggcorrplot(
      cordata,
      p.mat = NULL,
      hc.order = F,
      sig.level = 0.95,
      insig = "pch",
      pch.cex = 3
    )
  
  
}
p1
ggsave(file.path(directorio, "correlacion_Vars_Pcs.jpg"), p1)

# Contributions of variables to PC1
p1 <- fviz_contrib(pcx,
                   choice = "var",
                   axes = 1,
                   top = 10)
# Contributions of variables to PC2
p2 <- fviz_contrib(pcx,
                   choice = "var",
                   axes = 3,
                   top = 10)

p3 <- ggarrange(p1, p2)
p3

ggsave(
  filename = file.path(directorio, "graficos_contribuciones_PCA.jpg"),
  plot = p3
)




correlacion  <- psych::corr.test(X,method = "spearman",adjust = "BH")
p.values <- correlacion$p
p.values[lower.tri(correlacion$p)] <- p.values[upper.tri(correlacion$p)]

diag(p.values) <- 1
jpeg(filename = file.path(directorio,"correlacion.jpeg"))
corrplot::corrplot(correlacion$r,tl.cex = 0.7,
                   p.mat = correlacion$p,
                   sig.level = 0.01,
                   insig = "label_sig")
dev.off
