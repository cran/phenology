#' .Lnegbin estimate a negative binomial likelihood.
#' @title The function ".Lnegbin"
#' @author Marc Girondot
#' @return Return the likelihood
#' @param x Set of parameters to be fitted
#' @param pt Transfer parameters
#' @description Function of the package phenology

.Lnegbin <- function(x, pt) {

#.phenology.env<- NULL
#rm(.phenology.env)

# pt=list(data=data, fixed=parametersfixed, incertitude=method_incertitude, zerocounts=zero_counts)

sum=0
# je mets tous les paramètres dans xpar
	xpar<-c(x, pt$fixed)


for(k in 1:length(pt$data)) {

# print(paste("pop", k))

datatot<-pt$data
data<-datatot[[k]]
# quel est le nom de la série en cours
nmser<-names(datatot[k])
# Prend en compte les 0 ou non 5/2/2012
zero<- pt$zerocounts[k]
deb<-ifelse(zero, 0, 1)
	
# je fais une fonction où j'envoie les paramètres et la série et il renvoie ceux à utiliser directement
	xparec<-.format_par(xpar, nmser)
#  if (xparec["MinE"]<0) print(xpar, xparec)
	th<-xparec["Theta"]

for (i in 1:dim(data)[1]) 
{
# est ce que j'ai une observation ?
# transformé en 0 ou non avec zero TRUE ou FALSE: 5/2/2012
       if ((data$nombre[i]!=0)||zero) {

       		if ((is.na(data$Date2[i]))) {
#######################
# on a une seule date #
#######################

# ici je devrais aussi faire le calcul en conditionnel
       			sumnbcount<-.daily_count(data$ordinal[i], xparec, print=FALSE)
# dans zero j'ai l'information si je prends les 0 ou non pour cette série
       			if (!zero) {
       				lnli2 <- -log(dnbinom(data$nombre[i], size=th, mu=sumnbcount, log=FALSE)/(1-dnbinom(0, size=th, mu=sumnbcount, log=FALSE)))
				} else {
       				lnli2 <- -dnbinom(data$nombre[i], size=th, mu=sumnbcount, log=TRUE)
 #      				if (is.na(lnli2)) print(xparec)
       			}

        	} else {
###################
# on a deux dates #
###################
            
# nombre de jours. On met des NA  
nbjour<-data$ordinal2[i]-data$ordinal[i]+1
nbcount<-rep(NA, nbjour)
nbcountrel<-nbcount

#countday prend les valeurs des jours possibles
for(countday in 1:nbjour) {
#je met le nombre théorique de ce jour dans nbcount
	nbcount[countday]<-.daily_count(countday+data$ordinal[i]-1,xparec, print=FALSE)
}

#je somme
sumnbcount<-sum(nbcount)

if (pt$incertitude==0 | pt$incertitude=="binomial") {

# je calcule pour chaque jour qu'elle est la fraction du nombre de ponte par rapport au total                
	nbcountrel<-nbcount/sumnbcount
	nbcountrel[nbcountrel==0]<-1E-9
	nbcountrel[nbcountrel==1]<-1-1E-9
	
	lnli2<-0
	
for(countday in 1:nbjour) {
	lnli<-0
	 for(countponte in deb:data$nombre[i]) {
	 	proba<-dbinom(countponte, size=data$nombre[i], prob=nbcountrel[countday], log=FALSE)
	 	if (!zero) {
       		dnb<-dnbinom(countponte, size=th, mu=nbcount[countday], log=FALSE)/(1-dnbinom(0, size=th, mu=nbcount[countday], log=FALSE))       			
		} else {
	 		dnb<-dnbinom(countponte, size=th, mu=nbcount[countday], log=FALSE)
       	}
	 	lnli<-lnli+proba*(dnb)
	}
	lnli2<-lnli2-log(lnli)
}

} else {

if (pt$incertitude==1 | pt$incertitude=="sum") {

# je suis sur la méthode 1 d'incertitude
if (!zero) {
	lnli2<--log(dnbinom(data$nombre[i], size=th, mu=sumnbcount, log=FALSE)/(1-dnbinom(0, size=th, mu=sumnbcount, log=FALSE)))
} else {
	lnli2<--dnbinom(data$nombre[i], size=th, mu=sumnbcount, log=TRUE)
}


} else {


# je suis sur la méthode 2 d'incertitude
	nbcountrel<-nbcount/sumnbcount
	nbcountrel[nbcountrel==0]<-1E-9
	nbcountrel[nbcountrel==1]<-1-1E-9

#dans nbcountrel j'ai les probabilités
#dans nbjour j'ai le nombre de jours
#dans data$nombre[i] j'ai le nombre de nids observés

#Indistinguishable Objects to Distinguishable Boxes
# http://2000clicks.com/mathhelp/CountingObjectsInBoxes.aspx
# The number of different ways to distribute n indistinguishable balls into
# k distinguishable boxes is C(n+k-1,k-1).
# For example, 5 balls into 3 boxes can be done in these C(7,2) = 21 ways:

N<-data$nombre[i]

# The number of different ways to distribute n indistinguishable balls into
# k distinguishable boxes is C(n+k-1,k-1).
# nb<-choose(N+nbjour-1,nbjour-1)=dim(tb)[1]
# divers<-matrix(rep(0, nbjour*nb), ncol=nbjour)

  # generate all possible positions of the boundaries
  xx <- combn(N+nbjour-1, nbjour-1)
  # compute the number of balls in each box
  a <- cbind(0, diag(nbjour)) - cbind(diag(nbjour), 0)
  tb<-t(a %*% rbind(0, xx, N+nbjour) - 1)


# je calcule déjà les dbnbinom pour toutes les solutions: 4/2/2012
# il ma faut un tableau de 1:nbjour et de 0:N
a <- try(matrix(rep(0, (N+1)*nbjour), nrow=N+1,ncol=nbjour), silent=TRUE)

if (class(a)=="try-error") {

	print("Too many incertitudes on the days. Use another incertitude method.")
	return(Inf)

}


for(ii in deb:N) {
	for(countday in 1:nbjour) {
		if (!zero) {
			a[ii+1, countday]<-log(dnbinom(ii, size=th, mu=nbcount[countday], log=FALSE)/(1-dnbinom(0, size=th, mu=nbcount[countday], log=FALSE)))
		} else {
			a[ii+1, countday]<-dnbinom(ii, size=th, mu=nbcount[countday], log=TRUE)
		}
	}
}


sump<-0
for(ii in 1:dim(tb)[1]) {
p<-dmultinom(tb[ii,1:nbjour], prob=nbcountrel, log=TRUE)
	for(countday in 1:nbjour) {
		p<-p+a[tb[ii,countday]+1, countday]
	}
sump <- sump+exp(p)
}
lnli2 <- -log(sump)


}

}

# fin du test if ((is.na(data$Date2[i]))) {
}

datatot[[k]]$LnL[i]<-lnli2
datatot[[k]]$Modeled[i]<-sumnbcount

# assign("data", datatot, envir=as.environment(.phenology.env))

sum<-sum+lnli2

# transformé en test de zero
}

# fin de la boucle des jours
}

# fin de la boucle des séries
}
return(sum)
# fin de la fonction
}
