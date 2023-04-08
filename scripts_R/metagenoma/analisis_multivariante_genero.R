library(reshape2)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(ggrepel)
library(ALDEx2)
library(vegan)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
datos <- readRDS("../../datos/preprocesado_05_02_23/novoom.rds")
set.seed(126581)

directorio <- "./scripts_R/metagenoma/resultados_genero_multivariante"

if(!dir.exists(directorio)){
  dir.create(directorio)
}

genero <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
obesidad <- datos$comunes$variables_in_bacteria$OBESE
grupo <- datos$comunes$variables_in_bacteria$GROUP
sexo <- datos$comunes$variables_in_bacteria$SEX
X <- genero

resultado <- analyze_multivariate(X,grupo,obesidad)



X.df <- as.data.frame(X)
X.df$obesidad <- obesidad
X.df$grupo <- grupo
X.df$sexo <- sexo
X.melt <- reshape2::melt(X.df)

all_grupos <- X.melt %>% group_by(obesidad,grupo) %>% summarise(sum=sum(exp(value)))
all_grupos$sum <- 100* all_grupos$sum/sum(all_grupos$sum)

all_obesidad <- X.melt %>% group_by(obesidad) %>% summarise(sum=sum(exp(value)))
all_obesidad$sum <- 100* all_obesidad$sum/sum(all_obesidad$sum)

all_grupo <- X.melt %>% group_by(grupo) %>% summarise(sum=sum(exp(value)))
all_grupo$sum <- 100* all_grupo$sum/sum(all_grupo$sum)

resultados_prop <- bind_rows(all_obesidad,all_grupo,all_grupos)



write.csv(resultados_prop,file.path(directorio,"resultados_proporciones.csv"))

write.csv(resultado,file.path(directorio,"resultados_multivariantes.csv"))


barras_all <- X.melt %>% group_by(obesidad,grupo) %>% summarise(suma=100*exp((value))/sum(exp(value)))

p1 <- ggplot(data=barras_all,aes(x=grupo,y=suma,fill=obesidad))+geom_bar(stat="identity",
                                                                         position = "dodge")+scale_fill_manual(values=c('steelblue','red'))+xlab("")+ylab("%")

barras_obesidad<- X.melt %>% group_by(obesidad) %>% summarise(suma=100*exp((value))/sum(exp(value)))

p2 <- ggplot(data=barras_obesidad,aes(x=obesidad,y=suma,fill="steelblue"))+geom_bar(stat="identity",
                                                                                    position = "dodge")+scale_fill_manual(values=c('steelblue'))+xlab("")+ylab("%")+theme(legend.position = "none")

barras_grupo<- X.melt %>% group_by(grupo) %>% summarise(suma=100*exp((value))/sum(exp(value)))

p3 <- ggplot(data=barras_grupo,aes(x=grupo,y=suma,fill="steelblue"))+geom_bar(stat="identity",
                                                                              position = "dodge")+scale_fill_manual(values=c('steelblue'))+xlab("")+ylab("%")+theme(legend.position = "none")


barras_sexo<- X.melt %>% group_by(sexo) %>% summarise(suma=100*exp((value))/sum(exp(value)))

p4 <- ggplot(data=barras_sexo,aes(x=sexo,y=suma,fill="steelblue"))+geom_bar(stat="identity",
                                                                            position = "dodge")+scale_fill_manual(values=c('steelblue'))+xlab("")+ylab("%")+theme(legend.position = "none")
p5 <- ggarrange(p1,p2,p3,p4)



p6 <- annotate_figure(p5,top="genero")

ggsave(filename = file.path(directorio,"graficos_proporciones.jpg"),plot = p6)

pcx <- prcomp((genero),scale.=T)
# Contributions of variables to PC1
p1 <-fviz_contrib(pcx, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
p2 <- fviz_contrib(pcx, choice = "var", axes = 2, top = 10)

p3 <- ggarrange(p1,p2)

ggsave(filename = file.path(directorio,"graficos_contribuciones_PCA.jpg"),plot = p3)


p1 <- fviz_screeplot(pcx, addlabels = TRUE, ylim = c(0, 15))

resultados.univariantes <- analyze_data(genero,grupo = grupo,obesidad)
nombres <- rownames(resultados.univariantes$p_valores)[resultados.univariantes$p_valores<0.01]
nombres <- nombres[complete.cases(nombres)]
p2 <- fviz_pca_var(pcx,col.var = "contrib",select.var = list(names=nombres))+theme(legend.position = "none")

p3 <- ggarrange(p1,p2)

ggsave(filename = file.path(directorio,"graficos_contribuciones_scree_plotPCA.jpg"),plot = p3)

principal_components <- data.frame(PC1=pcx$x[,1],
                                   PC2=pcx$x[,2],
                                   grupo=grupo,
                                   obesidad=obesidad,
                                   sexo=sexo)
varianzas <- round(100*pcx$sdev^2/sum(pcx$sdev^2),2)
p1<-ggplot(principal_components,aes(PC1,PC2,color=sexo))+geom_point()+facet_wrap(~obesidad)+ylab(paste("PC2",varianzas[2],"%"))+xlab("")
p2 <- ggplot(principal_components,aes(PC1,PC2,color=obesidad))+geom_point()+facet_wrap(~grupo)+xlab(paste("PC1",varianzas[1],"%"))+ylab(paste("PC2",varianzas[2],"%"))
p3 <- ggarrange(p1,p2,ncol=1,nrow=2)


p4 <- annotate_figure(p3,top='PCA')
ggsave(filename = file.path(directorio,"biplotPCA.jpg"),plot = p4)


plsda <- mixOmics::plsda(genero,grupo)

principal_components_2 <- data.frame(PC1=plsda$variates$X[,1],
                                     PC2=plsda$variates$X[,2],
                                     grupo=grupo,
                                     obesidad=obesidad,
                                     sexo=sexo)

p1 <- ggplot(principal_components_2,aes(PC1,PC2,color=grupo))+facet_grid(~obesidad)+geom_point()+xlab(paste("PC1",round(100*plsda$prop_expl_var$X[1],2),"%"))+ylab(paste("PC2",round(100*plsda$prop_expl_var$X[2],2),"%"))+ggtitle("PLS-DA")


csa2 <- mixOmics::plotLoadings(plsda,ndisplay=16,title="Loadings",ncomps=c(1,2),size.name=0.9,contrib="max")
csa <- mixOmics::plotLoadings(plsda,ndisplay=16,title="Loadings",ncomps=c(1,2),size.name=0.9,contrib="min",plot=F)
csa$variables <- rownames(csa)
csa2$variables <- rownames(csa)

p1 <- ggplot(csa,aes(x=variables,y=importance,fill=GroupContrib))+geom_bar(stat="identity",
                                                                           position = "dodge")+coord_flip()+xlab("")

p2 <- ggplot(csa2,aes(x=variables,y=importance,fill=GroupContrib))+geom_bar(stat="identity",
                                                                            position = "dodge")+coord_flip()+xlab("")

p3 <- ggarrange(p1,p2,ncol=1,nrow=2)

ggsave(filename = file.path(directorio,"contribuciones_plsda.jpg"),plot=p3)

p4 <- mixOmics::plotVar(plsda,cutoff = 0.6,plot=F)
p4$nombres <- rownames(p4)
p5 <- ggplot(p4,aes(x,y,label=nombres))+geom_point()+ geom_text_repel()+xlab(paste("PC1",round(100*plsda$prop_expl_var$X[1],2),"%"))

p1 <- ggplot(principal_components_2,aes(PC1,PC2,color=grupo))+facet_grid(~obesidad)+geom_point()+xlab(paste("PC1",round(100*plsda$prop_expl_var$X[1],2),"%"))+ylab(paste("PC2",round(100*plsda$prop_expl_var$X[2],2),"%"))

p2 <- ggarrange(p1,p5,common.legend = T)

p3 <- annotate_figure(p2,top="PLS-DA")

ggsave(filename = file.path(directorio,"loadings_plsda.jpg"),plot=p3)

# 
