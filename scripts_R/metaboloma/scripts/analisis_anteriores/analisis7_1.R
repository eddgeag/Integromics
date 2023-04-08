library(reshape2)
library(dplyr)
library(ggplot2)
library(limma)
library(car)
library(Hotelling)
library(factoextra)
library(FactoMineR)
library(ggpubr)

convert.factor <- function(x,nivel) as.factor(as.character(x[x!=nivel]))
convert.factor2 <- function(x,y,n1,n2) as.factor(as.character(x[x!=n1 & y==n2]))
convert.factor3 <- function(x,y,n1) as.factor(as.character(x[y==n1]))
sexo <- datos$totales$general_data$SEX
obesidad <- datos$totales$general_data$OBESE
grupo <- datos$totales$general_data$GROUP


metaboloma <- as.data.frame(t(datos$totales$metaboloma))

metaboloma.log <- scale(log(metaboloma),scale = F)

set.seed(123)
padonis <- adonis2(metaboloma.log ~ grupo*obesidad,method="euclidean")

### metaboloma grupo
met.fem.pcos <- metaboloma.log[grupo!="Male",]
met.fem.male <- metaboloma.log[grupo!="PCOS",]
met.pcos.male <- metaboloma.log[grupo!="Female",]

fem.pcos <- convert.factor(grupo,"Male")
fem.male <- convert.factor(grupo,"PCOS")
pcos.male <- convert.factor(grupo,"Female")

### metaboloma grupo|obesidad
met.fem.pcos.obese <- metaboloma.log[grupo!="Male" & obesidad=="Obese",]
met.fem.male.obese <- metaboloma.log[grupo!="PCOS" & obesidad=="Obese",]
met.pcos.male.obese <- metaboloma.log[grupo!="Female" & obesidad =="Obese",]

fem.pcos.obese <- convert.factor2(grupo,obesidad,"Male","Obese")
fem.male.obese <- convert.factor2(grupo,obesidad,"PCOS","Obese")
pcos.male.obese  <- convert.factor2(grupo,obesidad,"Female","Obese")

met.fem.pcos.Noobese <- metaboloma.log[grupo!="Male" & obesidad=="No Obese",]
met.fem.male.Noobese <- metaboloma.log[grupo!="PCOS" & obesidad=="No Obese",]
met.pcos.male.Noobese <- metaboloma.log[grupo!="Female" & obesidad =="No Obese",]


fem.pcos.Noobese <- convert.factor2(grupo,obesidad,"Male","No Obese")
fem.male.Noobese <- convert.factor2(grupo,obesidad,"PCOS","No Obese")
pcos.male.Noobese  <- convert.factor2(grupo,obesidad,"Female","No Obese")
### metaboloma obesidad|grupo

met.pcos <- metaboloma.log[grupo=="PCOS",]
met.females <- metaboloma.log[grupo=="Female",]
met.male <- metaboloma.log[grupo=="Male",]

pcos <- convert.factor3(obesidad,grupo,"PCOS")
female <- convert.factor3(obesidad,grupo,"Female")
male <- convert.factor3(obesidad,grupo,"Male")

### Hacemos el analisis del grupo

pmanova <-  function(x,y) adonis2(x ~ y,method = "euclidean")$`Pr(>F)`[1]
p.fem.pcos <- pmanova(met.fem.pcos,fem.pcos)
p.fem.male <- pmanova(met.fem.male,fem.male)
p.pcos.male <- pmanova(met.pcos.male,pcos.male)

p.adj_grupo <- p.adjust(c(p.fem.pcos=p.fem.pcos,
                          p.fem.male=p.fem.male,
                          p.pcos.male=p.pcos.male),"bonferroni")

## la interaccion

p.fem.male.obese <- pmanova(met.fem.male.obese,fem.male.obese)
p.fem.pcos.obese <- pmanova(met.fem.pcos.obese,fem.pcos.obese)
p.pcos.male.obese <- pmanova(met.pcos.male.obese,pcos.male.obese)

p.fem.male.Noobese <- pmanova(met.fem.male.Noobese,fem.male.Noobese)
p.fem.pcos.Noobese <- pmanova(met.fem.pcos.Noobese,fem.pcos.Noobese)
p.pcos.male.Noobese <- pmanova(met.pcos.male.Noobese,pcos.male.Noobese)

p.pcos <- pmanova(met.pcos,pcos)
p.female <- pmanova(met.females,female)
p.males <- pmanova(met.male,male)

p.adj_interaccion <- p.adjust(c(p.fem.male.obese=p.fem.male.obese,
                                p.fem.pcos.obese=p.fem.pcos.obese,
                                p.pcos.male.obese=p.pcos.male.obese,
                                p.fem.male.Noobese=p.fem.male.Noobese,
                                p.fem.pcos.Noobese=p.fem.pcos.Noobese,
                                p.pcos.male.Noobese=p.pcos.male.Noobese,
                                p.pcos=p.pcos,
                                p.female=p.female,
                                p.males=p.males),method = "bonferroni")



p.value.df <- data.frame(p.valor=c(obesidad=padonis$`Pr(>F)`[2],
                                   grupo=padonis$`Pr(>F)`[1],
                                   interaccion=padonis$`Pr(>F)`[3],
                                   p.fem.pcos=p.adj_grupo[1],
                                   p.fem.male=p.adj_grupo[2],
                                   p.pcos.male=p.adj_grupo[3],
                                   p.fem.male.obese=p.adj_interaccion[1],
                                   p.fem.pcos.obese=p.adj_interaccion[2],
                                   p.pcos.male.obese=p.adj_interaccion[3],
                                   p.fem.male.Noobese=p.adj_interaccion[4],
                                   p.fem.pcos.Noobese=p.adj_interaccion[5],
                                   p.pcos.male.Noobese=p.adj_interaccion[6],
                                   p.pcos=p.adj_interaccion[7],
                                   p.female=p.adj_interaccion[8],
                                   p.males=p.adj_interaccion[9]))
p.value.df$interpretacion <- ifelse(p.value.df$p.valor<0.05,"**","ns")
if(!dir.exists("./12-2-22/scripts_analysis_reunion/scripts_metaboloma2/resultados_multivariantes/")){
  dir.create("./12-2-22/scripts_analysis_reunion/scripts_metaboloma2/resultados_multivariantes/")
}
write.csv(p.value.df,"./12-2-22/scripts_analysis_reunion/scripts_metaboloma2/resultados_multivariantes/totales.csv")
