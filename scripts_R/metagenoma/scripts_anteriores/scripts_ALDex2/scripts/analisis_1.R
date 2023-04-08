

library(ALDEx2)
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/manova_vanvalen.R")
source("./12-2-22/scripts_analysis_reunion/scripts_funciones/otras_funciones_utiles.R")

datos.obesidad <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_obesidad_.rds")

datos.grupo <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_grupo_.rds")

datos.sexo <- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_sex_.rds")

datos.interaccion<- readRDS("../../datos/resultados_analisis_datos_procesados_1000_novoom/resultados_analisis_datos_procesados_1000_novoom_interaccion_.rds")

obesidad <- datos.sexo$comunes$variables_in_bacteria$OBESE
grupo <- datos.sexo$comunes$variables_in_bacteria$GROUP
sexo <- datos.sexo$comunes$variables_in_bacteria$SEX

p.obesidad <- aldex.ttest(datos.sexo$comunes$microbiota$aldex$filo.aldex,verbose = T)
p.obesidad.filt <- aldex.ttest(datos.sexo$comunes$microbiota$aldex.filt$filo.aldex,verbose = T)
p.obesidad.genero <- aldex.ttest(datos.sexo$comunes$microbiota$aldex.filt$genero.aldex,verbose = T)
p.obesidad.genero.filt <- aldex.ttest(datos.sexo$comunes$microbiota$aldex$genero.aldex,verbose = T)




