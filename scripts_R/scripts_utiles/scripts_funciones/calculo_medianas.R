

colMedians <- function(x) apply(x,2,median)
compute_medians <- function(X,grupo,obesidad){
  
  obesos <- colMedians(X[obesidad=="Obese",])
  Noobesos <- colMedians(X[obesidad=="No Obese",])
  females <- colMedians(X[grupo=="Female",])
  males <- colMedians(X[grupo=="Male",])
  pcos <- colMedians(X[grupo=="PCOS",])
  females.o <- colMedians(X[grupo=="Female" & obesidad=="Obese",])
  males.o <- colMedians(X[grupo=="Male" & obesidad=="Obese",])
  pcos.o <- colMedians(X[grupo=="PCOS" & obesidad=="Obese",])
  females.noo <- colMedians(X[grupo=="Female" & obesidad=="No Obese",])
  males.noo <- colMedians(X[grupo=="Male" & obesidad=="No Obese",])
  pcos.noo <- colMedians(X[grupo=="PCOS" & obesidad=="No Obese",])
  
  medianas <- data.frame(obesos=obesos,
                         No.Obesos=Noobesos,
                         females=females,
                         males=males,
                         pcos=pcos,
                         females.obesas=females.o,
                         males.obesos=males.o,
                         pcos.obesas=pcos.o,
                         females.NoObesas=females.noo,
                         males.NoObesos=males.noo,
                         pcos.NoObesas=pcos.noo)
  
  return(medianas)

  
}
