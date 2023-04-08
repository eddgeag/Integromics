
# BiocManager::install("mixOmics")
# BiocManager::install("BiocParallel")
library(mixOmics)
load("datos.RData")
variables_in_bacteria = salvar$variables
datos = salvar$datos
res1.pls.met_micro <- mixOmics::pls(datos$metabolome,datos$microbiome,ncomp=1)
(rho_met_micro <- cor(res1.pls.met_micro$variates$X,res1.pls.met_micro$variates$Y))

res1.pls.met_IP <- mixOmics::pls(datos$IP,datos$metabolome,ncomp = 1)
(rho_met_ip <- cor(res1.pls.met_IP$variates$X,res1.pls.met_IP$variates$Y))

res1.pls.met_clin <- mixOmics::pls(datos$metabolome,datos$clinical_data,ncomp=1)
(rho_pls_met_clin <- cor(res1.pls.met_clin$variates$X,res1.pls.met_clin$variates$Y))

res1.pls_micro_clin <- mixOmics::pls(datos$microbiome,datos$clinical_data,ncomp=1)
(rho_pls_micro_clin <- cor(res1.pls_micro_clin$variates$X,res1.pls_micro_clin$variates$Y))

res1.pls_micro_IP <- mixOmics::pls(datos$microbiome,datos$IP,ncomp = 1)
(rho_pls_micro_IP <- cor(res1.pls_micro_IP$variates$X,res1.pls_micro_IP$variates$Y))

res1.pls_IP_clin <- mixOmics::pls(datos$IP,datos$clinical_data,ncomp = 1)
(rho_pls_IP_clin <- cor(res1.pls_IP_clin$variates$X,res1.pls_IP_clin$variates$Y))

design =matrix(c(0,rho_pls_IP_clin,rho_pls_met_clin,rho_pls_micro_clin,
                 rho_pls_IP_clin,0,rho_met_ip,rho_pls_micro_IP,
                 rho_pls_met_clin,rho_met_ip,0,rho_met_micro,
                 rho_pls_micro_clin,rho_pls_micro_IP,rho_met_micro,0),ncol = 4,dimnames = list(names(datos),names(datos)))
interaccion <-with(variables_in_bacteria,interaction(GROUP,OBESE))
names(interaccion)<-variables_in_bacteria$Paciente
sgccda.res_total_1 = block.splsda(X = datos, Y=interaccion, ncomp = 5,
                                  design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run
perf.diablo_toltal_1 = mixOmics::perf(sgccda.res_total_1, validation = 'Mfold', folds = 5, nrepeat = 30)


perf.diablo_toltal_1$choice.ncomp$WeightedVote
ncomp = perf.diablo_toltal_1$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
list.keepX <- list(clinical_data=c(1:16),IP=c(1:5),metabolome=c(1:37),microbiome=c(1:83))
BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)
tune.splsda.srbct <- tune.block.splsda(X=datos, Y=interaccion, ncomp = ncomp, # we suggest to push ncomp a bit more, e.g. 4
                                  test.keepX = list.keepX,
                                 design=design,
                                 nrepeat = 50,
                                 BPPARAM=BPPARAM,
                                 fold=6)   # we suggest nrepeat = 50

save(tune.splsda.srbct,file="tune.RData")


