
directorio <- "./12-2-22/scripts_analysis_reunion/scripts_metagenoma_novoom/resultados/"
nombres <- list.files(directorio)
resultados <- lapply(list.files(directorio,full.names = T), read.csv)
names(resultados) <- nombres

directorio.datos <- "../../datos/preprocesado_05_02_23/novoom.rds"
datos <- readRDS(directorio.datos)

genero <- datos$comunes$microbiota$genero.abs

write.csv(colnames(genero),"~/Desktop/genero.csv")

resultados$obesidad.csv
