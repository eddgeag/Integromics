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

# set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run
perf.diablo_toltal_1 = mixOmics::perf(sgccda.res_total_1, validation = 'Mfold', folds = 3, nrepeat = 50)


perf.diablo_toltal_1$choice.ncomp$WeightedVote
ncomp1 = perf.diablo_toltal_1$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
list.keppX <- list(clinical_data=seq(1,ncol(datos$clinical_data),3),IP=c(1:ncol(datos$IP)),metabolome=c(seq(1,ncol(datos$metabolome),9)),
                   microbiome=c(seq(1,ncol(datos$microbiome),20)))
print(paste("numero de componentes",ncomp1))
print("aqui")
BPPARAM2 <- BiocParallel::MulticoreParam(workers = parallelly::availableCores())

tune.splsda.srbct <- tune.block.splsda(X=datos, Y=interaccion, ncomp = ncomp1, # we suggest to push ncomp a bit more, e.g. 4
                                       test.keepX=list.keppX,
                                       design=design,
                                       nrepeat=1,
                                       BPPARAM=BPPARAM2,
                                       fold=3,progressBar =T)   # we suggest nrepeat = 50


# Now tuning a new component given previous tuned keepX

already.tested.X = tune.splsda.srbct$choice.keepX
alred <- list(clinical_data=c(1,2,3,4),IP=c(2,3,4,5),metabolome=c(3,29,4,21),microbiome=c(4,1,44,3))
tune = tune.block.splsda(X = datos, Y = interaccion,nrepeat=1,
                         ncomp = ncomp1+1, test.keppX = list.keppX, design = design,
                         already.tested.X = already.tested.X,progressBar=T,fold=3)

ista <- list(tune$choice.ncomp,
             tune$choice.keepX)
save(tune,file="./lista.RData")
