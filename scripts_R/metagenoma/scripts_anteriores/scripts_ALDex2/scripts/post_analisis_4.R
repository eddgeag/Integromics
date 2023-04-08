
csa <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")

write.csv(colnames(csa$comunes$microbiota$filo.abs),"~/Desktop/filo_tmp.csv")

nombres_archivos <- list.files("./12-2-22/scripts_analysis_reunion/scripts_metagenoma_novoom/resultados_filo/")
resultados <- lapply(list.files("./12-2-22/scripts_analysis_reunion/scripts_metagenoma_novoom/resultados_filo/",full.names = T), read.csv)
names(resultados) <- nombres_archivos

resultados$sexo.csv
resultados$obesidad.csv
resultados$diferencias_grupo.csv

resultados$PCOS.vs.Male.csv
resultados$PCOS.vs.Female.csv
resultados$Male.vs.Female.csv

resultados$NoObesidad_mujeres.csv
resultados$NoObesidad_sexo_control.csv
resultados$Noobesidad_hombresPCOS.csv

resultados$obesidad_hombresPCOS.csv
resultados$obesidad_mujeres_pcos.csv
resultados$obesidad_sexo_control.csv


resultados$obesidad_mujeres.csv
resultados$Obesidad_PCOS.csv
resultados$obesidad_hombres.csv


