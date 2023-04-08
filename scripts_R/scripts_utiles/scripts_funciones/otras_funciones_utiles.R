
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
funcion_plot <-
  function(X.,
           nombre_analisis.,
           obesidad_,
           grupo_,
           sexo_,
           interaccion_,
           directorio) {
    
    pval.obesidad <- funcion_multivariante2(
      X. ,
      categorica = obesidad_,
      categoricas = list(grupo_, obesidad_),
      nombre_analisis = nombre_analisis.
    )[[1]]
    
    pval.sexo <- funcion_multivariante2(
      X.,
      categorica = sexo_,
      categoricas = list(grupo_, obesidad_),
      nombre_analisis = nombre_analisis.
    )[[1]]
    
    pval.grupo <- funcion_multivariante2(
      X.,
      categorica = grupo_,
      categoricas = list(grupo_, obesidad_),
      nombre_analisis = nombre_analisis.
    )[[1]]
    
    pval.interaccion <- funcion_multivariante2(
      X.,
      categorica = grupo_,
      categoricas = list(grupo_, obesidad_),
      nombre_analisis = nombre_analisis.
    )[[2]]
    
    datos <- X.
    datos[is.na(datos)] <- 0
    pcx_ <- prcomp(datos)
    pcx <- pcx_$x
    varianzas <- 100 * round(pcx_$sdev ^ 2 / sum(pcx_$sdev ^ 2), 4)
    mds <- as.data.frame(cmdscale(dist(datos)))
    colnames(mds) <- c("PC1", "PC2")
    plotear <- bind_cols(
      PC1 = pcx[, 1],
      PC2 = pcx[, 2],
      grupo = grupo_,
      obesidad = obesidad_,
      sexo = sexo_
    )
    plotear2 <- bind_cols(
      PC1 = mds[, 1],
      PC2 = mds[, 2],
      grupo = grupo_,
      obesidad = obesidad_,
      sexo = sexo_
    )
    
    pvals <-
      paste0(
        "obesidad ",
        round(pval.obesidad, 4),
        "; grupo ",
        round(pval.grupo, 4),
        "; sexo" ,
        round(pval.sexo, 4),
        "; interaccion ",
        round(pval.interaccion, 4)
      )
    datos. <- datos
    datos.$grupo <- grupo_
    datos.$obesidad <- obesidad_
    datos.$sexo <- sexo_
    datos.$interaccion <- interaccion_
    datos.melt <- melt(datos.)
    ## PCA
    p1 <-
      ggplot(datos.melt, aes(value, color = grupo)) + geom_density()
    p2 <-
      ggplot(datos.melt, aes(value, color = obesidad)) + geom_density()
    p3 <-
      ggplot(datos.melt, aes(value, color = sexo)) + geom_density()
    p4 <-
      ggplot(datos.melt, aes(value, color = interaccion)) + geom_density()
    p_Dens <- ggarrange(p1, p2, p3, p4)
    
    ggsave(plot = p_Dens, filename = file.path(directorio, paste0("Densidad", nombre_analisis., ".jpg")))
    
    p1 <-
      ggplot(datos.melt, aes(y = value, color = grupo)) + geom_boxplot()
    p2 <-
      ggplot(datos.melt, aes(y = value, color = obesidad)) + geom_boxplot()
    p3 <-
      ggplot(datos.melt, aes(y = value, color = sexo)) + geom_boxplot()
    p4 <-
      ggplot(datos.melt, aes(y = value, color = interaccion)) + geom_boxplot()
    p_Box <- ggarrange(p1, p2, p3, p4)
    ggsave(plot = p_Box, filename = file.path(directorio, paste0("Boxplot", nombre_analisis., ".jpg")))
    
    
    p1 <-
      ggplot(plotear, aes(PC1, PC2, color = obesidad, shape = grupo)) + geom_point() +
      xlab(paste("PC1", varianzas[1], "%")) + ylab(paste("PC2", varianzas[2], "%"))
    p2 <-
      ggplot(plotear, aes(PC1, PC2, color = grupo, shape = obesidad)) + geom_point() +
      xlab(paste("PC1", varianzas[1], "%")) + ylab(paste("PC2", varianzas[2], "%"))
    p3 <-
      ggplot(plotear, aes(PC1, PC2, color = obesidad)) + geom_point() + xlab(paste("PC1", varianzas[1], "%")) +
      ylab(paste("PC2", varianzas[2], "%")) + facet_grid( ~ grupo)
    p4 <-
      ggplot(plotear, aes(PC1, PC2, color = grupo, shape = obesidad_)) + geom_point() +
      xlab(paste("PC1", varianzas[1], "%")) + ylab(paste("PC2", varianzas[2], "%")) +
      facet_grid( ~ obesidad)
    
    p5 <- ggarrange(p1, p2, p3, p4)
    
    p_PCA <- annotate_figure(p5, paste("PCA", pvals))
    
    ggsave(plot = p_PCA, filename = file.path(directorio, paste0("PCA", nombre_analisis., ".jpg")))
    
    
    p1 <-
      ggplot(plotear2, aes(PC1, PC2, color = obesidad, shape = grupo)) + geom_point()
    p2 <-
      ggplot(plotear2, aes(PC1, PC2, color = grupo, shape = obesidad)) + geom_point()
    p3 <-
      ggplot(plotear2, aes(PC1, PC2, color = obesidad)) + geom_point() + facet_grid( ~
                                                                                       grupo)
    p4 <-
      ggplot(plotear2, aes(PC1, PC2, color = grupo, shape = obesidad)) + geom_point() +
      facet_grid( ~ obesidad)
    
    p5 <- ggarrange(p1, p2, p3, p4)
    p_PCOA  <- annotate_figure(p5, paste("PCoA", pvals))
    
    ggsave(plot = p_PCOA, filename = file.path(directorio, paste0("PCoA", nombre_analisis., ".jpg")))
    
    
    p1 <-
      ggplot(plotear2, aes(PC1, PC2, color = obesidad, shape = sexo)) + geom_point()
    p2 <-
      ggplot(plotear2, aes(PC1, PC2, color = sexo, shape = obesidad)) + geom_point()
    p3 <-
      ggplot(plotear2, aes(PC1, PC2, color = sexo)) + geom_point() + facet_grid( ~
                                                                                   obesidad)
    p4 <-
      ggplot(plotear2, aes(PC1, PC2, color = obesidad, shape = obesidad)) + geom_point() +
      facet_grid( ~ sexo)
    
    p5 <- ggarrange(p1, p2, p3, p4)
    p_PCOA  <- annotate_figure(p5, paste("PCoA_sexo", pvals))
    
    ggsave(plot = p_PCOA, filename = file.path(directorio, paste0("PCoA_sexo", nombre_analisis., ".jpg")))
    
    p1 <-
      ggplot(plotear, aes(PC1, PC2, color = obesidad, shape = sexo)) + geom_point()
    p2 <-
      ggplot(plotear, aes(PC1, PC2, color = sexo, shape = obesidad)) + geom_point()
    p3 <-
      ggplot(plotear, aes(PC1, PC2, color = sexo)) + geom_point() + facet_grid( ~
                                                                                  obesidad)
    p4 <-
      ggplot(plotear, aes(PC1, PC2, color = obesidad, shape = obesidad)) + geom_point() +
      facet_grid( ~ sexo)
    
    p5 <- ggarrange(p1, p2, p3, p4)
    p_PCOA  <- annotate_figure(p5, paste("PCA_sexo", pvals))
    
    ggsave(plot = p_PCOA, filename = file.path(directorio, paste0("PCA_sexo", nombre_analisis., ".jpg")))
    
    
    
  }


mean_aldex <- function(X) {
  mc.all <- getMonteCarloInstances(X)
  mc.instances <- numMCInstances(X)
  row_mean <- function(y)
    apply(y, 1, mean)
  X <- data.frame(bind_rows(lapply(mc.all, row_mean)))
  rownames(X) <- names(mc.all)
  return(as.data.frame(t(X)))
}

add_classificacion_metaboloma <- function(X) {
  aromatic <- rownames(X)[c(2, 3, 4, 1, 6, 35, 36, 34)]
  other <-
    rownames(X)[c(21, 9, 16, 26, 22, 18, 13, 29, 30, 27)]
  aa <- rownames(X)[c(33, 15, 20, 19, 23)]
  carbo <-
    rownames(X)[!rownames(X) %in% c(aromatic, other, aa)]
  
  X$grupo_meta <-
    ifelse(rownames(X) %in% aromatic, "AROMATICO", NA)
  X$grupo_meta <-
    ifelse(rownames(X) %in% other,
           "OTHER",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(rownames(X) %in% aa,
           "DERIVADOS AA",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(
      rownames(X) %in% carbo,
      "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
      X$grupo_meta
    )
  
  return(X)
}



whichnot <- function(x, nivel)
  as.factor(as.character(x[x != nivel]))
whichnot.ob <-
  function(x, nivel, obese)
    as.factor(as.character(x[x != nivel & obesidad == obese]))
whichis <- function(x, y, nivel)
  as.factor(as.character(y[x == nivel]))

## vemos permanova
analyze_multivariate <- function(X,grupo,obesidad){
  
  pmanova.int <-
    adonis2(X ~ grupo * obesidad, by = "margin", method = "euclidean")$`Pr(>F)`[1]
  pmanova.grupo <-
    adonis2(X ~ grupo, by = "margin", method = "euclidean")$`Pr(>F)`[1]
  
  ## vemos grupo
  fm <- adonis2(X[grupo != "PCOS", ] ~ whichnot(grupo, "PCOS"),
                method = "euclidean", by = "margin")$`Pr(>F)`[1]
  fp <- adonis2(X[grupo != "Male", ] ~ whichnot(grupo, "Male"),
                method = "euclidean", by = "margin")$`Pr(>F)`[1]
  mp <- adonis2(X[grupo != "Female", ] ~ whichnot(grupo, "Female"),
                method = "euclidean",
                by = "margin")$`Pr(>F)`[1]
  ## vemos obesidad
  obe <-
    adonis2(X ~ obesidad, method = "euclidean", by = "margin")$`Pr(>F)`[1]
  ## vemos grupo dado la obesidad
  
  fm.o <-
    adonis2(X[grupo != "PCOS" &
                obesidad == "Obese", ] ~ whichnot.ob(grupo, "PCOS", "Obese"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  fp.o <-
    adonis2(X[grupo != "Male" &
                obesidad == "Obese", ] ~ whichnot.ob(grupo, "Male", "Obese"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  mp.o <-
    adonis2(X[grupo != "Female" &
                obesidad == "Obese", ] ~ whichnot.ob(grupo, "Female", "Obese"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  
  
  fm.no <-
    adonis2(X[grupo != "PCOS" &
                obesidad == "No Obese", ] ~ whichnot.ob(grupo, "PCOS", "No Obese"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  fp.no <-
    adonis2(X[grupo != "Male" &
                obesidad == "No Obese", ] ~ whichnot.ob(grupo, "Male", "No Obese"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  mp.no <-
    adonis2(X[grupo != "Female" &
                obesidad == "No Obese", ] ~ whichnot.ob(grupo, "Female", "No Obese"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  
  
  ## vemos la obesidad dado el grupo
  
  pcos <- adonis2(X[grupo == "PCOS", ] ~ whichis(grupo, obesidad, "PCOS"),
                  method = "euclidean",
                  by = "margin")$`Pr(>F)`[1]
  
  females <-
    adonis2(X[grupo == "Female", ] ~ whichis(grupo, obesidad, "Female"),
            method = "euclidean",
            by = "margin")$`Pr(>F)`[1]
  
  males <- adonis2(X[grupo == "Male", ] ~ whichis(grupo, obesidad, "Male"),
                   method = "euclidean",
                   by = "margin")$`Pr(>F)`[1]
  
  
  resultado <- c(
    obesidad = obe,
    grupo = pmanova.grupo,
    females_pcos = fp,
    females_males = fm,
    males_pcos = mp,
    interaccion = pmanova.int,
    females_pcos_obese = fp.o,
    females_males_obese = fm.o,
    males_pcos_obese = mp.o,
    females_pcos_Noobese = fp.no,
    females_males_Noobese = fm.no,
    males_pcos_Noobese = mp.no,
    pcos = pcos,
    females = females,
    males = males
  )
  
  resultado.df <- data.frame(p.values=resultado)
  resultado.df$p.adj <- c(resultado.df$p.values[1:2],
                          p.adjust(resultado.df$p.values[3:4],"BH"),
                          resultado.df$p.values[5:6],
                          p.adjust(resultado.df$p.values[7:length(resultado.df$p.values)],
                                   "BH"))
  resultado.df$significancia <- ifelse(resultado.df$p.values<0.05,resultado.df$p.values,NA)
  resultado.df$significancia.adj <- ifelse(resultado.df$p.adj<0.05,resultado.df$p.adj,NA)
  
  return(resultado.df)
}

