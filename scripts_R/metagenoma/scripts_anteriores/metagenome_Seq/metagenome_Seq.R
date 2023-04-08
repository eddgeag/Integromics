library(metagenomeSeq)
library(biomformat)
library(limma)
library(ggplot2)
library(reshape2)

datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")

contajes <-datos$comunes$microbiota$genero.abs
obesidad <- datos$comunes$variables_in_bacteria$OBESE
grupo <- datos$comunes$variables_in_bacteria$GROUP
sexo <- datos$comunes$variables_in_bacteria$SEX
pheno <- data.frame(obesidad=obesidad,grupo=grupo,sexo=sexo,row.names = rownames(contajes))
taxa <- data.frame(OTU=colnames(contajes),row.names = colnames(contajes))
contajes.t <- as.data.frame(t(contajes))

phenotypeData = AnnotatedDataFrame(pheno)
OTUdata = AnnotatedDataFrame(taxa)
obj = newMRexperiment(contajes.t,phenoData=phenotypeData,featureData=OTUdata)

obj_norm = filterData(obj, present = 30, depth = 1)

obj_norm = wrenchNorm(obj_norm, condition = obj$sexo)
mod <- model.matrix( ~ 1 + sexo, data = pheno)
# lungres1 = fitFeatureModel(obj_norm, mod)

mat <- as.data.frame(t(MRcounts(obj_norm, norm = TRUE, log = TRUE)))
mat$sexo <- sexo
mat.melt <- melt(mat)

ggplot(mat.melt,aes(value,color=sexo))+geom_density()

mat <- as.data.frame(t(MRcounts(obj_norm, norm = TRUE, log = TRUE)))

pcx <- prcomp(mat)

tmp <- data.frame(PC1=pcx$x[,1],PC2=pcx$x[,2],sexo=sexo)
ggplot(tmp,aes(PC1,PC2,color=sexo))+geom_point()


  