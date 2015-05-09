#################################################################################
# Total number of nests for all RMUs being a constant
#################################################################################

.fonctionfitTcst<-function(x, fixed, data, index) {
  
  #___________________________________________________________
  # Je crée le vecteur avec les paramètres ajustés et les fixe
  #___________________________________________________________
  x<-c(x, fixed)
  
  dtaL_obs <- data[,index$mean]
  
  dtaL_SD <- data[,index$se]
  
  #	variance is mu + mu^2/size=V
  #  V - mu = mu^2/size
  #  size = (mu^2)/(V-mu)
  #	dtaL_SD<-(dtaL_obs^2)/(dtaL_SD^2-dtaL_obs)	
  #	dtaL_SD[is.na(dtaL_SD)]<-0
  dtaL_SD<-abs(dtaL_SD+x["Size"])
  
  nabeach <- colnames(data)[index$mean]
  #___________________________________________________________
  # Je crée le vecteur avec les proportions de chaque site
  #___________________________________________________________
  
  p <- c(abs(x[nabeach[nabeach %in% names(x)]]), 
         structure(100, .Names = nabeach[!(nabeach %in% names(x))]))
  pp<-p/sum(p)
  
  # L'ordre n'est pas forcément le bon
  # le bon ordre est dans names()
  # A revoir
  
  T<-abs(x["T"])
  
  dtaL_theo <- matrix(rep(T*pp, dim(data)[1]) , ncol=length(pp) , byrow = TRUE)
  colnames(dtaL_theo) <- names(pp)
  
  dtaL_theo <- dtaL_theo[, names(pp)]
  dtaL_obs <- dtaL_obs[, names(pp)]
  
  L<-sum(-dnorm(x=dtaL_obs[!is.na(dtaL_obs)], 
                mean=dtaL_theo[!is.na(dtaL_obs)],  
                sd=dtaL_SD[!is.na(dtaL_obs)], log=TRUE))
  
  
  
  #	if (!is.finite(L)) {print(x)}
  return(L)
}


#################################################################################
# r being constant
#################################################################################

.fonctionfitTexp<-function(x, fixed, data, index) {
  
  x<-c(x, fixed)
  
  dtaL_obs<-data[,index$mean]
  
  dtaL_SD<-data[,index$se]
  
  #	variance is mu + mu^2/size=V
  #  V - mu = mu^2/size
  #  size = (mu^2)/(V-mu)
  #	dtaL_SD<-(dtaL_obs^2)/(dtaL_SD^2-dtaL_obs)	
  #	dtaL_SD[is.na(dtaL_SD)]<-0
  dtaL_SD<-abs(dtaL_SD+x["Size"])
  
  nabeach <- colnames(data)[index$mean]
  #___________________________________________________________
  # Je crée le vecteur avec les proportions de chaque site
  #___________________________________________________________
  
  p <- c(abs(x[nabeach[nabeach %in% names(x)]]), 
         structure(100, .Names = nabeach[!(nabeach %in% names(x))]))
  pp<-p/sum(p)
  
  # L'ordre n'est pas forcément le bon
  # le bon ordre est dans names()
  # A revoir
  
  Tot <- abs(x["T"])*exp(x["r"]*(1:(dim(data)[1])))
  
  dtaL_theo <- matrix(rep(Tot, length(pp)) , ncol=length(pp) , byrow = FALSE)
  colnames(dtaL_theo) <- names(pp)
  
  for (j in 1:length(Tot))
    dtaL_theo[j,] <- dtaL_theo[j, ]*pp
  
#  dtaL_theo <- dtaL_theo[, names(pp)]
  dtaL_obs <- dtaL_obs[, names(pp)]
  
  L<-sum(-dnorm(x=dtaL_obs[!is.na(dtaL_obs)], 
                mean=dtaL_theo[!is.na(dtaL_obs)], 
                sd=dtaL_SD[!is.na(dtaL_obs)], log=TRUE))
  
  
  
  #	if (!is.finite(L)) {print(x)}
  return(L)
  
}


.fonctionfitYearspecific<-function(x, fixed, data, index) {
  
  x<-c(x, fixed)
  
  dtaL_obs<-data[,index$mean]
  
  dtaL_SD<-data[,index$se]
  
  #	variance is mu + mu^2/size=V
  #  V - mu = mu^2/size
  #  size = (mu^2)/(V-mu)
  #	dtaL_SD<-(dtaL_obs^2)/(dtaL_SD^2-dtaL_obs)	
  #	dtaL_SD[is.na(dtaL_SD)]<-0
  dtaL_SD<-abs(dtaL_SD+x["Size"])
  
  nabeach <- colnames(data)[index$mean]
  #___________________________________________________________
  # Je crée le vecteur avec les proportions de chaque site
  #___________________________________________________________
  
  p <- c(abs(x[nabeach[nabeach %in% names(x)]]), 
         structure(100, .Names = nabeach[!(nabeach %in% names(x))]))
  pp <- p/sum(p)
  
  Tot <- abs(x[paste("T", data[,index$year], sep="")])
  
  # je mets dans les colonnes les valeurs des années
  # il faut que je mette les proportions maintenant: pp
  
  dtaL_theo <- matrix(rep(Tot,length(pp)) , ncol=length(pp), byrow = FALSE)
  colnames(dtaL_theo) <- names(pp)
  
  # dtaL_theo <- dtaL_theo[, names(pp)]
  dtaL_obs <- dtaL_obs[, names(pp)]
  
  for(i in 1:length(Tot)) dtaL_theo[i, 1:length(pp)] <- dtaL_theo[i, 1:length(pp)]*pp
  
  L<-sum(-dnorm(x=dtaL_obs[!is.na(dtaL_obs)], 
                mean=dtaL_theo[!is.na(dtaL_obs)], 
                sd=dtaL_SD[!is.na(dtaL_obs)], log=TRUE))
  
  
  
  #	if (!is.finite(L)) {print(x)}
  return(L)
}


