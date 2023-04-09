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
datos <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")
set.seed(126581)

directorio <- "./scripts_R/metagenoma/resultados_filo_multivariante"

if (!dir.exists(directorio)) {
  dir.create(directorio)
}

mean_filo <- scale(t(mean_aldex(datos$comunes$microbiota$filo)))
filo.melt <- reshape2::melt(mean_filo)



ggplot(filo.melt, aes(value)) + geom_density()

obesidad <- datos$comunes$variables_in_bacteria$OBESE
grupo <- datos$comunes$variables_in_bacteria$GROUP
sexo <- datos$comunes$variables_in_bacteria$SEX
X <- mean_filo

resultado <-
  analyze_multivariate(t(mean_aldex(datos$comunes$microbiota$filo)), grupo, obesidad, sexo = sexo)



X.df <- as.data.frame(as.matrix(dist(X)))
X.df$obesidad <- obesidad
X.df$grupo <- grupo
X.df$sexo <- sexo
X.melt <- reshape2::melt(X.df)

all_grupos <-
  X.melt %>% group_by(obesidad, grupo) %>% summarise(sum = sum(exp(value)))
all_grupos$sum <- 100 * all_grupos$sum / sum(all_grupos$sum)

all_obesidad <-
  X.melt %>% group_by(obesidad) %>% summarise(sum = sum(exp(value)))
all_obesidad$sum <- 100 * all_obesidad$sum / sum(all_obesidad$sum)

all_grupo <-
  X.melt %>% group_by(grupo) %>% summarise(sum = sum(exp(value)))
all_grupo$sum <- 100 * all_grupo$sum / sum(all_grupo$sum)

resultados_prop <- bind_rows(all_obesidad, all_grupo, all_grupos)



write.csv(resultados_prop,
          file.path(directorio, "resultados_proporciones.csv"))

write.csv(resultado,
          file.path(directorio, "resultados_multivariantes.csv"))


barras_all <-
  X.melt %>% group_by(obesidad, grupo) %>% summarise(suma = mean(value))

p1 <-
  ggplot(data = barras_all, aes(x = grupo, y = suma, fill = obesidad)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = c('steelblue', 'red')) + xlab("") + ylab("%")

barras_obesidad <-
  X.melt %>% group_by(obesidad) %>% summarise(suma = mean(value))

p2 <-
  ggplot(data = barras_obesidad, aes(
    x = obesidad,
    y = suma,
    fill = "steelblue"
  )) + geom_bar(stat = "identity",
                position = "dodge") +
  scale_fill_manual(values = c('steelblue')) + xlab("") + ylab("%") + theme(legend.position = "none")

barras_grupo <-
  X.melt %>% group_by(grupo) %>% summarise(suma = mean(value))

p3 <-
  ggplot(data = barras_grupo, aes(
    x = grupo,
    y = suma,
    fill = "steelblue"
  )) + geom_bar(stat = "identity",
                position = "dodge") +
  scale_fill_manual(values = c('steelblue')) + xlab("") + ylab("%") + theme(legend.position = "none")


barras_sexo <-
  X.melt %>% group_by(sexo) %>% summarise(suma = mean(value))

p4 <-
  ggplot(data = barras_sexo, aes(
    x = sexo,
    y = suma,
    fill = "steelblue"
  )) + geom_bar(stat = "identity",
                position = "dodge") +
  scale_fill_manual(values = c('steelblue')) + xlab("") + ylab("%") + theme(legend.position = "none")
p5 <- ggarrange(p1, p2, p3, p4)



p6 <- annotate_figure(p5, top = "filo")

p6

ggsave(filename = file.path(directorio, "graficos_proporciones.jpg"),
       plot = p6)

pcx <- prcomp((X), scale. = F)
# Contributions of variables to PC1
p1 <- fviz_contrib(pcx,
                   choice = "var",
                   axes = 1,
                   top = 10)
# Contributions of variables to PC2
p2 <- fviz_contrib(pcx,
                   choice = "var",
                   axes = 2,
                   top = 10)

p3 <- ggarrange(p1, p2)
p3
ggsave(
  filename = file.path(directorio, "graficos_contribuciones_PCA.jpg"),
  plot = p3
)


p1 <- fviz_screeplot(pcx, addlabels = TRUE, ylim = c(0, 15))

resultados.univariantes <-
  analyze_data(filo, grupo = grupo, obesidad, correccion = 1)

nombres <-
  rownames(resultados.univariantes$p_valores)[which(resultados.univariantes$interpretacion !=
                                                      "")]
nombres <- nombres[complete.cases(nombres)]
p2 <-
  fviz_pca_var(pcx,
               col.var = "contrib",
               select.var = list(names = nombres)) + theme(legend.position = "none")
p2
p3 <- ggarrange(p1, p2)
p3
ggsave(
  filename = file.path(directorio, "graficos_contribuciones_scree_plotPCA.jpg"),
  plot = p3
)

principal_components <- data.frame(
  PC1 = pcx$x[, 1],
  PC2 = pcx$x[, 2],
  grupo = grupo,
  obesidad = obesidad,
  sexo = sexo
)
varianzas <- round(100 * pcx$sdev ^ 2 / sum(pcx$sdev ^ 2), 2)
p1 <-
  ggplot(principal_components, aes(PC1, PC2, color = sexo)) + geom_point() +
  facet_wrap(~ obesidad) + ylab(paste("PC2", varianzas[2], "%")) + xlab("")
p2 <-
  ggplot(principal_components, aes(PC1, PC2, color = obesidad)) + geom_point() +
  facet_wrap(~ grupo) + xlab(paste("PC1", varianzas[1], "%")) + ylab(paste("PC2", varianzas[2], "%"))
p3 <- ggarrange(p1, p2, ncol = 1, nrow = 2)


p4 <- annotate_figure(p3, top = 'PCA')
p4
ggsave(filename = file.path(directorio, "biplotPCA.jpg"),
       plot = p4)


## experimentar
set.seed(126581)

mean_filo <- scale(t(mean_aldex(datos$comunes$microbiota$filo)))

train <- sample(1:nrow(mean_filo), size = nrow(mean_filo) * 0.7)


train.data <- mean_filo[train,]
test.data <- mean_filo[-train,]

train.target <- obesidad[train]
train.test <- obesidad[-train]

metagenome.plsda <-
  mixOmics::plsda(train.data,
                  as.factor(as.character((train.target))),
                  ncomp = 8,
                  scale = F)

plsda.metagenome <-
  mixOmics::perf(
    metagenome.plsda,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # in

plsda.metagenome$choice.ncomp # what is the optimal value of components according to perf()

metagenome.plsda <-
  mixOmics::plsda(test.data, as.factor(((train.test))), ncomp = 1, scale =
                    F)


prediccion <-
  as.factor(predict(metagenome.plsda, test.data)$class$centroids.dist)


obese.caret <- caret::confusionMatrix(prediccion, train.test)
### grupo
train.target <- grupo[train]
train.test <- grupo[-train]

metagenome.plsda <-
  mixOmics::plsda(train.data,
                  as.factor(as.character((train.target))),
                  ncomp = 8,
                  scale = F)

plsda.metagenome <-
  mixOmics::perf(
    metagenome.plsda,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # in

plsda.metagenome$choice.ncomp # what is the optimal value of components according to perf()

metagenome.plsda <-
  mixOmics::plsda(test.data, (((train.test))), ncomp = 1, scale = F)


prediccion <-
  as.factor(predict(metagenome.plsda, test.data)$class$centroids.dist[, 1])

levels(prediccion) <- c("Female", "PCOS", "Male")

grupo.caret <- caret::confusionMatrix(prediccion, train.test)




### grupo
train.target <- sexo[train]
train.test <- sexo[-train]

metagenome.plsda <-
  mixOmics::plsda(train.data,
                  as.factor(as.character((train.target))),
                  ncomp = 8,
                  scale = F)

plsda.metagenome <-
  mixOmics::perf(
    metagenome.plsda,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # in

plsda.metagenome$choice.ncomp # what is the optimal value of components according to perf()

metagenome.plsda <-
  mixOmics::plsda(test.data, (((train.test))), ncomp = 1, scale = F)


prediccion <-
  as.factor(predict(metagenome.plsda, test.data)$class$centroids.dist[, 1])


sexo.caret <- caret::confusionMatrix(prediccion, train.test)


res.plsda <-
  data.frame(
    obesidad = obese.caret$overall,
    sexo = sexo.caret$overall,
    grupo = grupo.caret$overall
  )

train.target <- sexo[train]
train.test <- sexo[-train]
metagenome.plsda <-
  mixOmics::plsda(train.data,
                  as.factor(as.character((train.target))),
                  ncomp = 8,
                  scale = F)

plsda.metagenome <-
  mixOmics::perf(
    metagenome.plsda,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # in

ncopms <-
  plsda.metagenome$choice.ncomp[1, 2]# what is the optimal value of components according to perf()

metagenome.plsda <-
  mixOmics::plsda(mean_filo, sexo, ncomp = ncopms, scale = F)

plots_plsda <- data.frame(PC1 = metagenome.plsda$variates$X[, 1],
                          grupo = grupo,
                          obesidad = obesidad,
                          sexo)

ggplot(data = plots_plsda, aes(y = PC1, x = grupo, fill = obesidad)) + geom_boxplot() +
  ylab(paste("PC1", round(100 * metagenome.plsda$prop_expl_var$X, 2), "%")) +
  xlab("") + ggtitle("PLS-DA scores")


mod <- car::Anova(lm(plots_plsda$PC1 ~ grupo * obesidad))
pairwise.t.test(plots_plsda$PC1, grupo, "BH")

pairwise.t.test(plots_plsda$PC1, interaction(grupo, obesidad), "BH")


max.plsda <- (
  mixOmics::plotLoadings(
    metagenome.plsda,
    comp = 1,
    ndisplay = 10
    ,
    contrib = "max",
    plot = F,
    method = "median"
  )
)
max.plsda <- bind_cols(filos = rownames(max.plsda), max.plsda)
min.plsda <- (
  mixOmics::plotLoadings(
    metagenome.plsda,
    comp = 1,
    ndisplay = 10
    ,
    contrib = "min",
    plot = F,
    method = "median"
  )
)
min.plsda <- bind_cols(filos = rownames(min.plsda), min.plsda)


ggplot(max.plsda, aes(x = filos, y = importance, fill = GroupContrib)) +
  geom_bar(stat = "identity") + coord_flip() + 
  ylab("Loadings")+xlab("")




