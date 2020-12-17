#' Data generating process of MRVAR(n,p,S) 
#'
#' This function will generate data from an multi-regime staionary VAR(p) process and return a list containing data and parameters used in the VAR(p) process.
#'
#' @param n     : number of variables
#' @param S     : number of regimes 
#' @param p     : lag length, the lag length is the same across all regimes
#' @param T     : number of observations
#' @param SESVI : Index of the switching variable 
#' @param TH    : S-1 vector of threshold values 
#' @param Bo    : (n,n,p,S) array of the AR(p) process for each regime. If not given it will be generated.
#' @param Co    : (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an n vector for each regime, exog0 Co is a (n,k+1,S) array with first colume zero for each regime respectively,for exog1 Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao: (n,n,S)   array containing S covariance matrices of residuals
#' @param Uo    : residuals, if it is not NA it will be used as input to generate the VAR(p) for each regime.
#' @param SV    : exogeous switching variable
#' @param type  : type of the deterministic components type = ("none","const","exog0","exog1")   
#' @param X     : (T x k) matrix of exogeneous variables.
#' @param mu    : (n x S) matrix of the regime specific mean of the variables
#' @param Yo    : (p, n, S) array of initial values of the process
#' @param Do    : (T, n, S) array of extra exogenous components (not used with value zero) 
#' @param d     : lag of the self-exiting   

#'               (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
#' @return      A list containing the generated data, the parameters and the input exogeous variables. res = list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type)
#' 
#' @examples 
#' res_d = MRVARData4(n=2,p=2,T=300,S=2,type="exog1",SESVI=1,X=matrix(rnorm(2*300),300,2))
#' res_d = MRVARData4(n,p,T,S,SESVI,TH,Bo,Co,Sigmao,d=1,type="exog1",X=X)
#' res_d = MRVARData4(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao,d=1,type="exog1",X=X)
#' @export
MRVARData4 = function(n=2,p=2,T=100,S=1,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d=1) {
### T     : number of observations
### n     : number of variables
### p     : lag length
### S     : number of states of the underlying Markov Chain
###         (n,p,T,S,SESVI) are parameters which must be provided.
### SESVI : index of the switching variable in the endogeneous variables Y for the case of self-excited threshhold model.
### TH    : (S-1)-vector of threshold values
### Bo    : n x n x p x S array collecting the coefficients of VAR(p) in S different states 
### Sigmao: n x n x S array of the covariance matrix of the VAR(p) in S different states
### Uo    : an T x n x S  array of the temorally independent innovation processes
### SV    : exogeneous switching variable
### d     : lag of the self-exiting  
###         (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
###
### output:
### c("Y","Uo","Bo","Sigmao","TH","St","sv","SESVI","r_npo","check")
###
### Y     : the simulated data via of the MRVAR(p) 
### Uo    : T x n x S array of the simulated innovations in S different states
### Bo    : n x n x p x S array collecting the VAR(p) coefficients in S different states 
### Sigmao: n x n x S array collecting the covariance matrices of the simulated MRVAR(p) in S different states
### TH    : S-1 vector of threshold values
### St    : simulated time path of states
### sv    : the switching variable ( it is SW in the case of exogeneous switching and it is Y[,SESVI] in the case of self-excited threshold model.    
### r_npo : n x p matrix collecting the roots of the characteristic functions in L of the n dynamically independent processes.
### check : maximum of the data to check stationarity
###
### Remarks: The VAR(p) process is A transformed from n dynamically independent ar processes, such that the stationary is guaranteed. 
###          
###
  if (missing(Yo)) {
     Yo = NA
  }
  if (missing(Do)) {
     Do = NA
  }
  if (missing(d)) {
     d = NA
  }

  if (missing(TH)) {
      TH = NA
  }
  if (missing(Bo)) { Bo = NA }

  if (missing(Sigmao)) { Sigmao = NA  }

  if (missing(Uo))     { Uo = NA  }
  if (missing(type))   { type = NA }

  if (missing(Co))     { Co = NA  }  
  if (missing(mu))     { mu = NA  }
   
  if (missing(SV))     {SV = NA   }

  if (missing(X))      { X = NA  }
  if (missing(Yo))     { Yo = NA }

  if (missing(Do))     { Do = NA  }  
                       
  check = c(1:(S+1))*0
  if (anyNA(Yo)) {
     Yo = NA
  }
  if (anyNA(Do)) {
     Do = NA
  }

  if (anyNA(d)) {
     d = 1
  }

  P = max(p,d);

  if (anyNA(TH)) {
  TH = matrix(runif(S-1),S-1,1)
      TH = TH[order(TH)]
  }
  if (anyNA(Bo)) {
  	Bo = (1:(n*n*p*S))*0
      dim(Bo) = c(n,n,p,S)
      r_npo  = c(1:(n*p*S))*0
      dim(r_npo) = c(n,p,S)

      for (i in 1:S) {
      VARD = VARData(n,p,T) 
      Bo[,,,i]   = VARD$B
            r_npo[,,i] = VARD$r_np
            check[i]   = max(abs(VARD$Y))
      }      
      #dim(Bo) = c(n,n,p,S)
  }   else    {r_npo = NA    }



  if (anyNA(Sigmao)) {
  	Sigmao = (1:(n*n*S))*0
      dim(Sigmao) = c(n,n,S)
      for (i in 1:S) {
      	VARD = VARData(n,p,T) 
      	Sigmao[,,i] = VARD$Sigma
      }    
  }

  if (anyNA(Uo)) {
  	Uo = (1:(T*n*S))*0
      dim(Uo) = c(T,n,S)
      for (i in 1:S) Uo[,,i]= rnormSIGMA(T,as.matrix(Sigmao[,,i]))
  }
  
  Ct = Uo*0;

  if (anyNA(type)) {
  	type = "none"
      Co = NA
  }

  if (type=="exog0")  {
      m = dim(as.matrix(X))[2]
	CC = rnorm(n*(1+m)*S)
      dim(CC) = c(n,m+1,S); CC[,1,]=c(1:n)*0
  	if (anyNA(Co)) Co = CC  
      for (s in 1:S) Ct[,,s] = X%*%t(Co[,-1,s])      
  }

  if (type=="exog1")  {
      m = dim(as.matrix(X))[2]
	CC = rnorm(n*(1+m)*S)
      dim(CC) = c(n,m+1,S)
  	if (anyNA(Co)) Co = CC  
      for (s in 1:S) Ct[,,s] = cbind((1:T)/(1:T),X)%*%t(Co[,,s])      
  }   

  if (type == "const") {
  	if (anyNA(mu)) { mu = matrix(rnorm(n*S),n,S); dim(mu)=c(n,1,S)}
  	if (anyNA(Co)) { 
  		Co = mu
  		for (s in 1:S) { 
                  for (L in 1:p) Co[,1,s] = Co[,1,s] - Bo[,,L,s]%*%mu[,1,s]
                  Ct[,,s] = matrix(1,T,1)%*%t(Co[,1,s]) 
            	}
  	}  	else {
            mu = Co*NA

            for (s in 1:S) {
			H = diag(n)
			for (L in 1:p) H = H - Bo[,,L,s]
            	mu[,1,s] = solve(H)%*%Co[,1,s]
            }
      }
  }  
  

  St = (1:T)*0
  Y = as.matrix(Uo[,,1])*NA
  if (anyNA(Yo)) {
      Y[1:P,]=Uo[1:P,,1]
      Yo = Y[1:P,]      
  }   else {Y[1:P,] = Yo}
  
  if (anyNA(Do)) {Do = matrix(0,n,S)}  else  {Do = Do}   
  if (anyNA(SV)) {sv = Y[,SESVI]} else sv = SV
  
  

  for ( tt in (P+1):T )  { 
      if (anyNA(SV)) {sv = Y[,SESVI]} else {sv = SV}
      st = sv[tt-d];
      s  = sum(st>TH)+1
      St[tt] = s
      Y[tt,] = Uo[tt,,s]+Ct[tt,,s] + Do[,s]
      for (L in 1:p)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bo[,,L,s])
  }
  check[S+1] = max(abs(Y))
  resid = NA
  result=list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d)
  names(result) = c("Y","X","Uo","resid","Bo","Co","Sigmao","TH","St","sv","d","SESVI","r_npo","check","n","p","S","type","Yo","d")
  return(result)
}


#' Estimation of MRVAR(n,p,S) 
#'
#' This function estimate the unknown parameters of a specified MRVAR(n,p,S) model based on provided data.
#'
#' @param  res  :a list containing the components as the output of MRVARData4 including at least: n, p, type, Y, SESVI, TH, d, and optionally X. 
#' @return list(res,LH,vars_obj) res is a list as the input, but filled with estimated parameter values, LH is the estimated likelihood. vars_obj are a list contain S estimated VAR object in order to use the utilities of vars package. 
#' @examples 
#' res_d = MRVARData4(n=2,p=2,T=300,S=2,type="exog1",SESVI=1,X=matrix(rnorm(2*300),300,2))
#' res_d = MRVARData4(n,p,T,S,SESVI,TH,Bo,Co,Sigmao,d=1,type="exog1",X=X)
#' res_d = MRVARData4(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao,d=1,type="exog1",X=X)
#' RESS  = MRVARest1(res_d)
#' M1 <-   RESS$vars_obj[[1]]
#' plot(M1)
#' summary(M1)
#' IRF = irf(M1,boot=F)
#' plot(IRF)
#' @export
MRVARest1 = function(res) {
##
## This is a program used to estiamte MRVAR parameters under given order parameter of MRVAR
##
TH = res$TH
type = res$type
p     = res$p
S= res$S
n= res$n
Y= res$Y
X = res$X
T= dim(Y)[1]
SV    = res$SV
sv    = res$sv
SESVI = res$SESVI
b    = res$Bo*NA
cc    = res$Co*NA
sigma = res$Sigmao*NA

resid = (1:(T*n*S))*0; dim(resid) = c(T,n,S) 

Z = embed(Y,p+1)
if (type=="none")  {Z = embed(Y,p+1); cc = matrix(0,n,S);dim(cc)=c(n,1,S)}
if (type=="const")  Z = cbind(Z,(1:(T-p))/(1:(T-p)))
if (type=="exog0")  Z = cbind(Z,X[(p+1):T,])
if (type=="exog1")  Z = cbind(Z,((p+1):T)/((p+1):T),X[(p+1):T,])

m <- ncol(Z)
vars_obj = list()
St = (1:T)*0
if (is.null(SV)) {sv = Y[,SESVI]} else sv = SV 
for ( tt in (p+1):T )  { 
      svt = sv[tt-1];
      s = which.min(abs(TH-svt))
      if (svt>TH[s]) {s = s+1} else {s = s} 
      St[tt] = s
}
select = matrix(0,T-p,S)
select2 = select
for (s in 1:S) {
	select[,s] =  (1:(T-p))*(St[(p+1):T]==s)
      select2[,s] = ((p+1):T)*(St[(p+1):T]==s)
      if (!is.null(dim(as.matrix(Z[select[,s],1:n]))[1])) {
      	if (dim(as.matrix(Z[select[,s],1:n]))[1]-dim(as.matrix(Z[select[,s],(n+1):m]))[2]>0) {
      
      		LREG = lm(Z[select[,s],1:n]~0+Z[select[,s],(n+1):m])
     			bs = as.matrix(LREG$coefficients)[1:(n*p),]
      		if (type=="none")   cs = as.matrix(LREG$coefficients)[1,]*0
      		if (type=="const")  cs = as.matrix(LREG$coefficients)[m-n,]
      		if (type=="exog0")  cs = as.matrix(LREG$coefficients)[(n*p+1):(m-n),]
      		if (type=="exog1" ) cs = as.matrix(LREG$coefficients)[(n*p+1):(m-n),]
   
      		sigmas = t(LREG$residuals)%*%(LREG$residuals)/(dim(as.matrix(Z[select[,s],1:n]))[1]-dim(as.matrix(Z[select[,s],(n+1):m]))[2])
      		bs = t(bs); 
      		dim(bs) = c(n,n,p)
      		b[,,,s] = bs
      		if (type=="none")    cc[,1,s]      = (1:n)*0
      		if (type=="const" )  cc[,1,s]      = cs
      		if (type=="exog0") { cc[,-1,s]   = t(cs); cc[,1,s] = (1:n)*0 }
      		if (type=="exog1" )  cc[,,s]     = t(cs)


      		sigma[,,s] = sigmas
            	resid[select2[,s],,s] = LREG$residuals
      	}
      }
	#### the following code generate a list of S vars objects in order to use vars utilities
	Ys = as.matrix(Z[select[,s],1:n]); 
      Zs = as.matrix(Z[select[,s],(n+1):m]); 
      #Xs = Z[select[,s],(m-ncol(X)):m];

      YY = Y[1:(nrow(as.matrix(Z[select[,s],1:n]))+p),]; 
      if ((res$type=="exog0")|(res$type=="exog1")  ) XX = X[1:(nrow(Z[select[,s],1:n])+p),]

      if (n>1) {
		if (res$type=="none")    var_yy = vars::VAR(y=YY, p = p, type = "none")
		if (res$type=="const")   var_yy = vars::VAR(y=YY, p = p, type = "const")
		if (res$type=="exog0")   var_yy = vars::VAR(y=YY, p = p, type = "none",exogen = XX)
		if (res$type=="exog1")   var_yy = vars::VAR(y=YY, p = p, type = "const",exogen = XX)
     

		for ( i in 1:n) {
      		Names_of_coefficits = names(var_yy$varresult[[i]]$coefficients)
            	LREG = lm(Ys[,i]~0+Zs)
            	#str(var_yy$varresult)
   			var_yy$varresult[[i]]<-LREG
            	names(var_yy$varresult[[i]]$coefficients)[1:n]<-Names_of_coefficits
		}
		#str(var_yy$datamat)
		var_yy$datamat[,] <- Z[select[,s],]
	
      
      	if (p>1) { var_yy$y[,]       <- t(cbind(t(Ys[1:p,]),t(Ys)))         }
         	else  { var_yy$y[,]       <- t(cbind(as.matrix(Ys[1:p,]),t(Ys))) }
      }  else  { var_yy = lm(Ys~0+Zs) }
  
      vars_obj[[s]] = var_yy
}

   res$Bo      <- b
   res$Co      <- cc
   res$Sigmao  <- sigma
   res$St      <- St
   res$resid   <- resid
   LH = 0

      for (s in 1:S) { 

if (!is.na(log(det(as.matrix(sigma[,,s]))))) LH = LH + -dim(as.matrix(Z[select[,s],1:n]))[1]*log(det(as.matrix(sigma[,,s])))
}

result = list(res,LH,vars_obj)
names(result) = c("res","LH","vars_obj")
return(result)
}


#' Data generating process of MRVAR(n,p,S) 
#'
#' This function will generate data from an multi-regime staionary VAR(p) process and return a list containing data and parameters used in the VAR(p) process.
#'
#' @param n     : number of variables
#' @param S     : number of regimes 
#' @param p     : an S x 2 matrix. Each row contains the lag length for the corresponding state and the number of exogenous variables for the state.
#' @param T     : number of observations
#' @param SESVI : Index of the switching variable 
#' @param TH    : S-1 vector of threshold values 
#' @param Bo    : (n,n,p,S) array of the coeffiicent matrices of VAR(p) for each lag and each regime. If not given it will be generated.
#' @param Co    : (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an n vector for each regime, exog0 Co is a (n,k+1,S) array with first colume zero for each regime respectively,for exog1 Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao: (n,n,S)   array containing S covariance matrices of residuals, each for one regime.
#' @param Uo    : residuals, if it is not NA it will be used as input to generate the VAR(p) for each regime.
#' @param SV    : exogeous switching variable
#' @param type  : type of the deterministic components type = ("none","const","exog0","exog1")   
#' @param X     : (T x k x S) matrix of exogeneous variables for each state. The second dimension can be filled with zeros to take into account that the the exogenous variables are not identical in each state. 
#' @param mu    : (n x S) matrix of the regime specific mean of the variables
#' @param Yo    : (p, n, S) array of initial values of the process
#' @param Do    : (T, n, S) array of extra exogenous components (not used with value zero) 
#' @param d     : lag of the self-exiting switching variable.
#'   
#'               (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
#' @return      A object of MRVAR(n,p,S) containing the generated data, the parameters and the input of exogeous variables. res = list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type)
#' 
#' @examples 
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' res_e = MRVARest(res=res_d)
#' summary_MRVAR(res_e)
#' @export
MRVARData = function(n=2,p=matrix(2,2,2),T=100,S=1,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d=1) {
### T     : number of observations
### n     : number of variables
### p     : lag length, an S x 2 vector of lag length for each state and the number of exogenous variables for each state
### S     : number of states of the underlying Markov Chain
###         (n,p,T,S,SESVI) are parameters which must be provided.
### SESVI : index of the switching variable in the endogeneous variables Y for the case of self-excited threshhold model.
### TH    : (S-1)-vector of threshold values
### Bo    : n x n x p x S array collecting the coefficients of VAR(p) in S different states 
### Sigmao: n x n x S array of the covariance matrix of the VAR(p) in S different states
### Uo    : an T x n x S  array of the temorally independent innovation processes
### SV    : exogeneous switching variable
### d     : lag of the self-exiting  
###         (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
###
### output:
### c("Y","Uo","Bo","Sigmao","TH","St","sv","SESVI","r_npo","check")
###
### Y     : the simulated data via of the MRVAR(p) 
### Uo    : T x n x S array of the simulated innovations in S different states
### Bo    : n x n x p x S array collecting the VAR(p) coefficients in S different states 
### Sigmao: n x n x S array collecting the covariance matrices of the simulated MRVAR(p) in S different states
### TH    : S-1 vector of threshold values
### St    : simulated time path of states
### sv    : the switching variable ( it is SW in the case of exogeneous switching and it is Y[,SESVI] in the case of self-excited threshold model.    
### r_npo : n x p matrix collecting the roots of the characteristic functions in L of the n dynamically independent processes.
### check : maximum of the data to check stationarity
###
### Remarks: The VAR(p) process is A transformed from n dynamically independent ar processes, such that the stationary is guaranteed. 
###          
###
  if (missing(Yo)) {
     Yo = NA
  }
  if (missing(Do)) {
     Do = NA
  }
  if (missing(d)) {
     d = NA
  }

  if (missing(TH)) {
      TH = NA
  }
  if (missing(Bo)) { Bo = NA }

  if (missing(Sigmao)) { Sigmao = NA  }

  if (missing(Uo))     { Uo = NA  }
  if (missing(type))   { type = NA }

  if (missing(Co))     { Co = NA  }  
  if (missing(mu))     { mu = NA  }
   
  if (missing(SV))     {SV = NA   }

  if (missing(X))      { X = NA  }
  if (missing(Yo))     { Yo = NA }

  if (missing(Do))     { Do = NA  }  
                       
  check = c(1:(S+1))*0
  if (anyNA(Yo)) {
     Yo = NA
  }
  if (anyNA(Do)) {
     Do = NA
  }

  if (anyNA(d)) {
     d = 1
  }

  P = max(p[,1],d);
  Pmax = max(p[,1])
  if (anyNA(TH)) {
  TH = matrix(runif(S-1),S-1,1)
      TH = TH[order(TH)]
  }
  if (anyNA(Bo)) {
  	Bo = (1:(n*n*Pmax*S))*0
      dim(Bo) = c(n,n,Pmax,S)
      r_npo  = c(1:(n*Pmax*S))*0
      dim(r_npo) = c(n,Pmax,S)

      for (i in 1:S) {
      VARD = VARData(n,p[i],T) 
      Bo[,,1:p[i,1],i]   = VARD$B
            r_npo[,1:p[i,1],i] = VARD$r_np
            check[i]   = max(abs(VARD$Y))
      }      
      #dim(Bo) = c(n,n,p,S)
  }   else    {r_npo = NA    }



  if (anyNA(Sigmao)) {
  	Sigmao = (1:(n*n*S))*0
      dim(Sigmao) = c(n,n,S)
      for (i in 1:S) {
      	VARD = VARData(n,Pmax,T) 
      	Sigmao[,,i] = VARD$Sigma
      }    
  }

  if (anyNA(Uo)) {
  	Uo = (1:(T*n*S))*0
      dim(Uo) = c(T,n,S)
      for (i in 1:S) Uo[,,i]= rnormSIGMA(T,as.matrix(Sigmao[,,i]))
  }
  
  Ct = Uo*0;

  if (anyNA(type)) {
  	type = "none"
      Co = NA
  }

  if (type=="exog0")  {
      k = dim(X)[2]
	CC = rnorm(n*(1+k)*S)
      dim(CC) = c(n,k+1,S); CC[,1,]=c(1:n)*0
  	if (anyNA(Co)) Co = CC  
      for (s in 1:S) {
            if (p[s,2]<k) Co[,(p[s,2]+2):(k+1),s] = Co[,(p[s,2]+2):(k+1),s]*0
		Ct[,,s] = X[,,s]%*%t(Co[,-1,s])
      }      
  }

  if (type=="exog1")  {
      k = dim(X)[2]
	CC = rnorm(n*(1+k)*S)
      dim(CC) = c(n,k+1,S)
  	if (anyNA(Co)) Co = CC  
      for (s in 1:S) { 
            if (p[s,2]<k) Co[,(p[s,2]+2):(k+1),s] = Co[,(p[s,2]+2):(k+1),s]*0
		Ct[,,s] = cbind((1:T)/(1:T),X[,,s])%*%t(Co[,,s])      
      }
  }   

  if (type == "const") {
  	if (anyNA(mu)) { mu = matrix(rnorm(n*S),n,S); dim(mu)=c(n,1,S)}
  	if (anyNA(Co)) { 
  		Co = mu
  		for (s in 1:S) { 
                  for (L in 1:Pmax) Co[,1,s] = Co[,1,s] - Bo[,,L,s]%*%mu[,1,s]
                  Ct[,,s] = matrix(1,T,1)%*%t(Co[,1,s]) 
            	}
  	}  	else {
            mu = Co*NA

            for (s in 1:S) {
			H = diag(n)
			for (L in 1:Pmax) H = H - Bo[,,L,s]
            	mu[,1,s] = solve(H)%*%Co[,1,s]
            }
      }
  }  
  

  St = (1:T)*0
  Y = as.matrix(Uo[,,1])*NA
  if (anyNA(Yo)) {
      Y[1:P,]=Uo[1:P,,1]
      Yo = Y[1:P,]      
  }   else {Y[1:P,] = Yo}
  
  if (anyNA(Do)) {Do = matrix(0,n,S)}  else  {Do = Do}   
  if (anyNA(SV)) {sv = Y[,SESVI]} else sv = SV
  
  

  for ( tt in (P+1):T )  { 
      if (anyNA(SV)) {sv = Y[,SESVI]} else {sv = SV}
      st = sv[tt-d];
      s  = sum(st>TH)+1
      St[tt] = s
      Y[tt,] = Uo[tt,,s]+Ct[tt,,s] + Do[,s]
      for (L in 1:Pmax)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bo[,,L,s])
  }
  check[S+1] = max(abs(Y))
  resid = NA
  result=list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d)
  names(result) = c("Y","X","Uo","resid","Bo","Co","Sigmao","TH","St","sv","d","SESVI","r_npo","check","n","p","S","type","Yo","d")
  return(result)
}


#' Estimation of MRVAR(n,p,S) 
#'
#' This function estimates the unknown parameters of a specified MRVAR(n,p,S) based on provided data.
#'
#' @param  res  :a object of MRVAR which is an output of MRVARData including at least: n, S, p, type, Y, SESVI, TH, d, and optionally X. 
#' @return list(res, LH_AIC, LH_BIC, LH_P, LH_N, LH_NN,LLH_AIC, LLH_BIC, ORBIC, ORAIC, vars_obj): res is an object of MRVAR containing the estimated parameter values, LH is the estimated likelihood. vars_obj is a list contain S estimated VAR object in order to use the utilities of vars package. 
#' @examples 
#'
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' res_e = MRVARest(res=res_d)
#' res_e$LH_AIC
#' res_e$LH_BIC
#' res_e$ORBIC
#' res_e$ORAIC
#' res_e$LH_P
#' res_e$LH_N
#' summary_MRVAR(res_e)
#' RESS = res_e
#' #### IRF without regime migration
#' IRF = irf_MRVAR_NM1(RESS,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM1_CB2(RESS,nstep=20,comb=NA,irf="gen",runs=100,conf=c(0.05,0.90))
#' 
#' par(mfrow=c(4,2))
#' plott(IRF_CB[,,,,1],1,1)
#' plott(IRF_CB[,,,,1],1,2)
#' plott(IRF_CB[,,,,1],2,1)
#' plott(IRF_CB[,,,,1],2,2)
#' plott(IRF_CB[,,,,2],1,1)
#' plott(IRF_CB[,,,,2],1,2)
#' plott(IRF_CB[,,,,2],2,1)
#' plott(IRF_CB[,,,,2],2,2)
#'
#' @export
MRVARest = function(res) {
##
## This is a program used to estiamte MRVAR parameters under given order parameter of MRVAR
    TH = res$TH
    type = res$type
    p = res$p
    S = res$S
    n = res$n
    Y = res$Y
    X = res$X
    T = dim(Y)[1]
    SV = res$SV
    sv = res$sv
    d = res$d
    SESVI = res$SESVI
    P = max(p[, 1], d)
    Pmax = max(p[, 1])
    ms = (1:S) * 0
    LH_AIC = 0
    LH_BIC = 0
    LLH_AIC = 0
    LLH_BIC = 0
    LH_P = 0
    LH_N = 0
    LH_NN = 0
    b = res$Bo * NA
    cc = res$Co * NA
    sigma = res$Sigmao * NA
    resid = (1:(T * n * S)) * 0
    dim(resid) = c(T, n, S)
    Z = embed(Y, Pmax + 1)
    if (!anyNA(X)) {
        XX = JointX(X)
        if (type == "none") {
            ZO = embed(Y, Pmax + 1)
            cc = matrix(0, n, S)
            dim(cc) = c(n, 1, S)
        }
        if (type == "const") 
            ZO = cbind(Z, (1:(T - Pmax))/(1:(T - Pmax)))
        if (type == "exog0") 
            ZO = cbind(Z, XX[(Pmax + 1):T, ])
        if (type == "exog1") 
            ZO = cbind(Z, ((Pmax + 1):T)/((Pmax + 1):T), XX[(Pmax + 
                1):T, ])
        mso = dim(ZO)[2]
        LREGOR = lm(ZO[, 1:n] ~ 0 + ZO[, (n + 1):mso])
    }
    if (!anyNA(X)) {
        XX = JointX(X)
        Co = matrix(1, n, dim(XX)[2] + 1)
        if (type == "exog0") 
            Co[, 1] = 0
        VARDataHelp = VARData(n = n, p = Pmax, T = T, type = type, 
            Co = Co, X = XX)
    }
    else {
        if (type == "none") 
            Co = matrix(0, n, 1)
        if (type == "const") 
            Co = matrix(1, n, 1)
        VARDataHelp = VARData(n = n, p = Pmax, T = T, type = type, 
            Co = Co)
    }
    VARDataHelp$Y = res$Y
    ORBIC = VARest(res = VARDataHelp)$BIC
    ORAIC = VARest(VARDataHelp)$AIC
    m <- ncol(Z)
    mm = c(1:S) * 0
    vars_obj = list()
    St = (1:T) * 0
    if (is.null(SV)) {
        sv = Y[, SESVI]
    }
    else sv = SV
    for (tt in (P + 1):T) {
        svt = sv[tt - d]
        s = which.min(abs(TH - svt))
        if (svt > TH[s]) {
            s = s + 1
        }
        else {
            s = s
        }
        St[tt] = s
    }
    select = matrix(0, T - P, S)
    select2 = select
    for (s in 1:S) {
        select[, s] = (1:(T - P)) * (St[(P + 1):T] == s)
        select2[, s] = ((P + 1):T) * (St[(P + 1):T] == s)
        Z = embed(Y, P + 1)[, 1:((p[s, 1] + 1) * n)]
        if (type == "none") {
            Z = Z
            cc = matrix(0, n, S)
            dim(cc) = c(n, 1, S)
        }
        if (type == "const") 
            Z = cbind(Z, (1:(T - P))/(1:(T - P)))
        if (type == "exog0") 
            Z = cbind(Z, X[(P + 1):T, 1:p[s, 2], s])
        if (type == "exog1") 
            Z = cbind(Z, ((P + 1):T)/((P + 1):T), X[(P + 1):T, 
                1:p[s, 2], s])
        ms[s] = dim(Z)[2]
        if (!is.null(dim(as.matrix(Z[select[, s], 1:n]))[1])) {
            if (dim(as.matrix(Z[select[, s], 1:n]))[1] - dim(as.matrix(Z[select[, 
                s], (n + 1):((p[s, 1] + 1) * n)]))[2] > 0) {
                LREG = lm(Z[select[, s], 1:n] ~ 0 + Z[select[, 
                  s], (n + 1):ms[s]])
                bs = as.matrix(LREG$coefficients)[1:(n * p[s, 
                  1]), ]
                if (type == "none") 
                  cs = as.matrix(LREG$coefficients)[1, ] * 0
                if (type == "const") 
                  cs = as.matrix(LREG$coefficients)[ms[s] - n, 
                    ]
                if (type == "exog0") 
                  cs = as.matrix(LREG$coefficients)[(n * p[s, 
                    1] + 1):(ms[s] - n), ]
                if (type == "exog1") 
                  cs = as.matrix(LREG$coefficients)[(n * p[s, 
                    1] + 1):(ms[s] - n), ]
                sigmas = t(LREG$residuals) %*% (LREG$residuals)/(dim(as.matrix(Z[select[, 
                  s], 1:n]))[1] - dim(as.matrix(Z[select[, s], 
                  (n + 1):ms[s]]))[2])
                bs = t(bs)
                dim(bs) = c(n, n, p[s, 1])
                b[, , 1:p[s, 1], s] = bs
                if (type == "none") 
                  cc[, 1, s] = (1:n) * 0
                if (type == "const") 
                  cc[, 1, s] = cs
                if (type == "exog0") {
                  cc[, c(2:(1 + p[s, 2])), s] = t(cs)
                  cc[, 1, s] = (1:n) * 0
                }
                if (type == "exog1") 
                  cc[, c(1:(1 + p[s, 2])), s] = t(cs)
                sigma[, , s] = sigmas
                resid[select2[, s], , s] = LREG$residuals
            }
        }
        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LH_AIC = LH_AIC  + 2 * n * (dim(Z[,(n + 1):ms[s]])[2] + (n + 1)/2) + dim(as.matrix(Z[select[, s], 1:n]))[1] * log(det(as.matrix(sigma[, , s]))) + dim(as.matrix(Z[select[,s], 1:n]))[1] * n * (1 + log(2 * pi)) 


        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LH_BIC = LH_BIC + log(dim(as.matrix(Z[select[,s], 1:n]))[1]) * n * (dim(Z[, (n + 1):ms[s]])[2] + (n + 1)/2) + dim(as.matrix(Z[select[, s], 1:n]))[1] * log(det(as.matrix(sigma[, , s]))) + dim(as.matrix(Z[select[, s], 1:n]))[1] * n * (1 + log(2 * pi)) 
 

        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LH_P =  LH_P +     dim(as.matrix(Z[select[,s], 1:n]))[1] * log(det(as.matrix(sigma[, , s]))) +  dim(as.matrix(Z[select[,s], 1:n]))[1] * n * (1 + log(2 * pi))

        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LH_N = LH_N  +    2 * n * (dim(Z[,(n + 1):ms[s]])[2] + (n + 1)/2) 

        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LH_NN = LH_NN +  log(dim(as.matrix(Z[select[,s], 1:n]))[1]) * (n * (dim(Z[, (n + 1):ms[s]])[2] + (n + 1)/2)) 
   

        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LLH_AIC = LH_P + LH_N

        if (!is.na(log(det(as.matrix(sigma[, , s]))))) 
            LLH_BIC = LH_P + LH_NN 

       
        Ys = as.matrix(Z[select[, s], 1:n])
        Zs = as.matrix(Z[select[, s], (n + 1):ms[s]])
        YY = Y[1:(nrow(as.matrix(Z[select[, s], 1:n])) + p[s, 
            1]), ]
        if ((res$type == "exog0") | (res$type == "exog1")) 
            XX = X[1:(nrow(Z[select[, s], 1:n]) + p[s, 1]), 1:p[s, 
                2], s]
        if (n > 1) {
            if (res$type == "none") 
                var_yy = vars::VAR(y = YY, p = p[s, 1], type = "none")
            if (res$type == "const") 
                var_yy = vars::VAR(y = YY, p = p[s, 1], type = "const")
            if (res$type == "exog0") 
                var_yy = vars::VAR(y = YY, p = p[s, 1], type = "none", 
                  exogen = XX)
            if (res$type == "exog1") 
                var_yy = vars::VAR(y = YY, p = p[s, 1], type = "const", 
                  exogen = XX)
            for (i in 1:n) {
                Names_of_coefficits = names(var_yy$varresult[[i]]$coefficients)
                LREG = lm(Ys[, i] ~ 0 + Zs)
                var_yy$varresult[[i]] <- LREG
                names(var_yy$varresult[[i]]$coefficients) <- Names_of_coefficits
            }
            var_yy$datamat[, ] <- Z[select[, s], ]
            if (p[s, 1] > 1) {
                var_yy$y[, ] <- t(cbind(t(Ys[1:p[s, 1], ]), t(Ys)))
            }
            else {
                var_yy$y[, ] <- t(cbind(as.matrix(Ys[1:p[s, 1], 
                  ]), t(Ys)))
            }
        }
        else {
            var_yy = lm(Ys ~ 0 + Zs)
        }
        vars_obj[[s]] = var_yy
    }
    res$Bo <- b
    res$Co <- cc
    res$Sigmao <- sigma
    res$St <- St
    res$resid <- resid
    result = list(res, LH_AIC, LH_BIC, LH_P, LH_N, LH_NN,LLH_AIC, LLH_BIC, ORBIC, ORAIC, vars_obj)
    names(result) = c("res", "LH_AIC", "LH_BIC", "LH_P", "LH_N","LH_NN", "LLH_AIC", "LLH_BIC","ORBIC", "ORAIC", "vars_obj")
    return(result)
}


#' @export
JointX = function(X) {
### X: the T x Nx x S array collecting the exogeneous variables in MRVAR with state dependent exog varibales
	S  = dim(X)[3]
      Nx = dim(X)[2]
      ss = 0
      sidx = 0
	for (s in 1:S) {
        if (sum(colSums(X[,,s] != 0) !=0)>ss ) {
           ss = sum(colSums(X[,,s] != 0) !=0)
           sidx = s
      }  
      }
      XX = X[,1:ss,sidx]
      return(XX)
}




#' Generalized impulse response functions of MRVAR(n,p,S) with regime migrations 
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S).
#'  
#' ### for a given shock vector SHCK.
#"
#' ### GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK)) 
#'
#'#############################################################################################
#' @param  RES   : a list containing the components as the output of MRVARest1.
#' @param  shock : an n vector containing the shocks as impulse.
#' @param  R     : the number runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  state : the state from which we start the shocks
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric boostrap. 
#' @return a n x n x nstep+1 matrix of impulse response functions. The rows represent response the columns represent impulses. 
#' @examples 
#'
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' res_e = MRVARest(res=res_d)
#' RF3 = girf_MRVAR_RM(RES=res_e,shock=c(1,1),R = 800,nstep=20,state=1,resid_method="parametric")
#' RF4 = girf_MRVAR_RM_CB(RES=res_e, shock=c(1,1), R=800, nstep=20, state=1, resid_method = "parametric", conf_level=c(0.05,0.95), N=300)
#' x11()
#' par(mfrow=c(2,2))
#' plot(RF3[1,1,],type="l")
#' plot(RF3[1,2,],type="l")
#' plot(RF3[2,1,],type="l")
#' plot(RF3[2,2,],type="l")
#' dim(RF4)
#' x11()
#' par(mfrow=c(2,2))
#' plott(RF4[,,,],1,1)
#' plott(RF4[,,,],1,2)
#' plott(RF4[,,,],2,1)
#' plott(RF4[,,,],2,2)
#'
#' @export
girf_MRVAR_RM = function(RES,shock,R,nstep,state,resid_method) {
####  this function generate the impulse response funciton of MRVAR with migration
####
####                  	GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK)) 
####   
####  
#### RES is the output of MRVARest 
      n = RES$res$n
      p = RES$res$p
      S = RES$res$S
      SESVI= RES$res$SESVI; if (is.null(SESVI)) SESVI = NA
      TH= RES$res$TH;    if (is.null(TH))       TH = NA
      Bo= RES$res$Bo
      Co= RES$res$Co
      Sigmao= RES$res$Sigmao
	d= RES$res$d
      type= RES$res$type
      P           = max(p,d);
      X= RES$res$X;  if (!anyNA(X))  X = X[1:(P+nstep+1),]
      Yo= RES$res$Yo; if (is.null(Yo)) Yo = NA
      Uo          = RES$res$Uo
      res = RES$res
      vars_obj = RES$vars_obj
      TT          = dim(RES$res$Y)[1]

	IRF = list()
      YR          = list()
      YS          = list()
      residR   <-  list()
      residS   <-  residR
      IRFHLPR  = (1:(n*(P+nstep+1)*R))*0;   dim(IRFHLPR) = c(n,P+nstep+1,R) 
      MeanR       = (1:(n*(P+nstep+1)))*0;     dim(MeanR)   = c(n,P+nstep+1)
	IRFHLPS  = (1:(n*n*(P+nstep+1)*R))*0; dim(IRFHLPS) = c(n,n,P+nstep+1,R)
      MeanS       = (1:(n*n*(P+nstep+1)))*0;   dim(MeanS)   = c(n,n,P+nstep+1)

      #R = 400
      #shock       = (1:n)/(1:n)
      DIMresid    = dim(Uo)
      if (length(DIMresid) == 3)  residI = Uo[1:(P+1+nstep),,]*0 
      if (length(DIMresid) == 2)  residI = Uo[1:(P+1+nstep),]*0
      shock_mat = matrix(0,n,n)
      for (i in 1:n) {   
          shock_mat[i,i]  = shock[i]; 
          shock_mat[-i,i] = rnormSIGMA_cond(T=1,sigma=Sigmao[,,state],I=i, v=shock[i],switch=0)
      }
    	for (i in 1:R) { 
      if ( resid_method=="resid" ) {
	if (length(DIMresid) == 2)  { 
	residI  = Uo[NoutofT(P+nstep+1,DIMresid[1]),]; 
            }
            if (length(DIMresid) == 3)  {
            for (j in 1:DIMresid[3]) {
		residI[,,j]  = Uo[P+NoutofT(nstep+1,DIMresid[1]),,j];                        
                  }
            } 
      }
      if ( resid_method=="parametric" ) {
	if (length(DIMresid) == 2)  {
		residI   = rnormSIGMA(P+nstep+1,Sigmao)
            residI[1:P,]= residI[1:P,]*0 
	}
            if (length(DIMresid) == 3)  {
		for (j in 1:DIMresid[3]) {
                        residI[,,j] = rnormSIGMA(P+nstep+1,Sigmao[,,j])
                        residI[1:P,,j]= residI[1:P,,j]*0
                  }
            }
      }
      residR[[i]] = residI
      residS[[i]] = residI 
	Yo = RES$res$Y[1:P,]*0 
	YR[[i]] = MRVARData(n=n,p=p,T=(P+nstep+1),S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo=residR[[i]],type=type,X=X,Yo=Yo,d=d)
	IRFHLPR[,,i] = t(YR[[i]]$Y)

	for (k in 1:DIMresid[2]) {
            if (length(DIMresid) == 2)  residS[[i]][P+1,]  = t(shock_mat[,k])  
		if (length(DIMresid) == 3)  residS[[i]][P+1,,] = t(as.matrix((1:DIMresid[3])/(1:DIMresid[3]))%*%shock_mat[,k])
                  
            YS[[i]] = MRVARData(n=n,p=p,T=(P+nstep+1),S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo=residS[[i]],d=d,type=type,X=X,Yo=Yo)
                  IRFHLPS[,k,,i] = t(YS[[i]]$Y)

    }
    }  # end of R loop
    ############## integrating out the disturbances
      for (i in 1:n ) {
          MeanR[i,] = rowMeans(IRFHLPR[i,,],dims=1)
    	for (j in 1:n) {
                 MeanS[i,j,]      = rowMeans(IRFHLPS[i,j,,],dims=1)
                 MeanS[i,j,]      = MeanS[i,j,] - MeanR[i,] 
      }
      }
   return(MeanS)
}


#' @export
NoutofT = function(N,T) {
  unique(round(runif(3*N)*(T-1)))[1:N]+1
}


#' Regime specific impulse response functions of MRVAR(n,p,S) with confidence intervals
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S) for orhtogoalized shocks as well as for generalized irf. 
#' Using the estimated Bo[,,,s] and Sigma[,,s] matrices of the estimated MRVAR, this function calculated the regime speicfic impulse response functions.
#' @param B (n,n,p) of the sth regime. 
#' @param sigma the covariance matrix of the s-th regime
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a matrix of (n,n,nstep) as the IRF colunms respesenting the impulse rows the responses. 
#' @examples 
#' var_d = VARData(n=2,p=2,T=200)
#' IRF_G = impulsdtrf_g(B=var_d$B,sigma=var_d$Sigma,nstep=11,comb=NA,irf="gen1")
#' plot(IRF_G[1,1,])

#' res_d = MRVARData4(n=2,p=2,T=300,S=2,type="exog1",SESVI=1,X=matrix(rnorm(2*300),300,2))
#' RESS  = MRVARest1(res_d)
#' RF3 = girf_MRVAR_RM(RES=RESS,shock=c(1,1),R = 1800,nstep=10,state=1,resid_method="parametric")
#' IRF = impulsdtrf_g(B=RESS$res$Bo[,,,1],sigma=RESS$res$Sigmao[,,1],nstep=11,comb=NA,irf="gen1")
#' plot(IRF[1,1,])
#' lines(RF3[1,1,],col="green")
#'
#'
#' @export
impulsdtrf_g=function(B,sigma,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1"))
### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
###                 impulse responses" it should have the property that crossprod(smat)=sigma,
###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
###                 way to get such a smat is to set smat=chol(sigma).  To get the smat
###                 corresponding to a different ordering, use
###                 smat=chol(P %*% Sigma %*% t(P)) %*% P, where P is a permutation matrix.
###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
###                 positive powers of L in A), W=A_0^(-1)'.
###                 In general, though, it is not required that W be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
### response:       nvar(response) x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
  {
    ##-----debug--------
    ##browser()
    ##------------------
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    if ( irf == "chol" )    { smat  = chol(sigma)};
    if ( irf == "gen"  )    { smat  = t(sigma%*%diag(1/(sqrt(diag(sigma))))) };
    if ( irf == "chol1")    { smat  = chol(sigma)%*%diag(1/diag(chol(sigma)))};
    if ( irf == "gen1" )    { smat  = t(sigma%*%diag(1/(diag(sigma)))) };
    if ( irf == "concert" ) { smat  = t((sigma%*%diag(1/(diag(sigma))))%*%comb) };

    if ( irf == "comb") { 
               DD = diag(t(comb)%*%sigma%*%comb);
                           for (i in 1:length(DD)) { if (DD[i]>0) DD[i]=sqrt(1/DD[i]) }
                           DD = diag(DD)
                           smat  = t(sigma%*%comb%*%DD)

         #smat  = t(sigma%*%comb*as.numeric((1/(sqrt(t(comb[,1])%*%sigma%*%comb[,1])))))
 
                   };
    if ( irf == "comb1") { 
               DD = diag(t(comb)%*%sigma%*%comb);
                           for (i in 1:length(DD)) { if (DD[i]>0) DD[i]=(1/DD[i]) }
                           DD = diag(DD)
                           smat  = t(sigma%*%comb%*%DD)

         #smat  = t(sigma%*%comb*as.numeric((1/(sqrt(t(comb[,1])%*%sigma%*%comb[,1])))))
 
                   };


    if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
    response <- array(0,dim=c(neq,nvar,nstep));
    response[,,1] <- t(smat);
    for (it in 2:nstep)
      {
        for (ilag in 1:min(lags,it-1))
          response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
      }
    dimnames(response) <- list(dimnames(B)[[2]],dimnames(smat)[[1]],as.character(0:(nstep-1)))
    ## dimnames(response)[2] <- dimnames(smat)[1]
    ## dimnames(response)[1] <- dimnames(B)[2]
    return(response)
  }

#' Impulse response functions of MRVAR(n,p,S) with regime migrations and confidence intervals
#'
#' This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) for given a shock vector.
#'                   GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK)) 
#' It also calculated the boostraped confidence intervals.
#' @param  RES   : a list containing the components as the output of MRVARest1.
#' @param  shock : an n vector containing the shocks as impulse.
#' @param  R     : the number runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  state : the state from which we start the shocks
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric boostrap. 
#' @param  conf_level : a vecter contain the level of confidences 
#' @N            : the number of runs in ontaining the bootstraped confidence intervals.
#' @return a n x n x nstep+1 x 3 containing of impulse response functions with lower and upper confidence bonds. The rows represent response the columns represent impulses. 
#' @examples 

#'
##' @export
girf_MRVAR_RM_CB = function (RES, shock, R, nstep, state, resid_method = "parametric", conf_level, N) 
{
    n = RES$res$n
    p = RES$res$p
    S = RES$res$S
    SESVI = RES$res$SESVI
    TH = RES$res$TH
    Bo = RES$res$Bo
    Co = RES$res$Co
    Sigmao = RES$res$Sigmao
    d = RES$res$d
    type = RES$res$type
    X = RES$res$X
    T = dim(RES$res$Y)[1]
    P = max(p, d)
    GIRF = (1:(n * n * (P + nstep + 1) * N)) * 0
    dim(GIRF) = c(n, n, P + nstep + 1, N)
    GIRFBd = (1:(n * n * (P + nstep + 1) * (length(conf_level) + 1))) * 0
    dim(GIRFBd) = c(n, n, P + nstep + 1, length(conf_level) + 1)
    GIRFBd[, , , 1] = girf_MRVAR_RM(RES, shock = shock, R = R, nstep, state=state, resid_method)
    for (i in 1:N) {
        res_d = MRVARData(n = n, p = p, T = T, S = S, SESVI = SESVI, TH = TH, Bo = Bo, Co = Co, Sigmao, type = type, X = X, d = d)
        RESS = MRVARest(res_d)
        RF3 = girf_MRVAR_RM(RES = RESS, shock = shock, R = R, nstep, state=state, resid_method)
        GIRF[, , , i] = RF3
    }
    GIRF[,,,1] = GIRFBd[, , , 1] 
    for (tt in 1:(P + nstep + 1)) {
        for (i in 1:n) {
            for (j in 1:n) {
                GIRFBd[i, j, tt, -1] = quantile(GIRF[i, j, tt, ], conf_level)
            }
        }
    }
    return(GIRFBd)
}



#' Summary of MRVAR estimation results 
#'
#' This function sumerizes the estimation results of MRVAR in its regime specific VAR model
#'
#' @param  res  :a list of the output of MRVARest1
#' @examples 
#'
#' res_d = MRVARData4(n=2,p=2,T=300,S=2,SESVI=1)
#' RESS  = MRVARest1(res_d)
#' summary_MRVAR(RESS)
#'
#' @export
summary_MRVAR = function(res) {
 var_obj = res$vars_obj
  n = length(var_obj)
  for (i in 1:n) {
    print(summary(var_obj[[i]]))
  }
}




#' Initial U1 for DATA generating IRF
#' @export
U1=function(B,sigma,irf=c("gen","chol","chol1","gen1","comb1"))
  {
    ##-----debug--------
    ##browser()
    ##------------------
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    if ( irf == "chol" )    { smat  = chol(sigma)};
    if ( irf == "gen"  )    { smat  = t(sigma%*%diag(1/(sqrt(diag(sigma))))) };
    if ( irf == "chol1")    { smat  = chol(sigma)%*%diag(1/diag(chol(sigma)))};
    if ( irf == "gen1" )    { smat  = t(sigma%*%diag(1/(diag(sigma)))) };
    if ( irf == "concert" ) { smat  = t((sigma%*%diag(1/(diag(sigma))))%*%comb) };
    if ( irf == "comb") { 
               DD = diag(t(comb)%*%sigma%*%comb);
                           for (i in 1:length(DD)) { if (DD[i]>0) DD[i]=sqrt(1/DD[i]) }
                           DD = diag(DD)
                           smat  = t(sigma%*%comb%*%DD)

         #smat  = t(sigma%*%comb*as.numeric((1/(sqrt(t(comb[,1])%*%sigma%*%comb[,1])))))
 
                   };
    if ( irf == "comb1") { 
               DD = diag(t(comb)%*%sigma%*%comb);
                           for (i in 1:length(DD)) { if (DD[i]>0) DD[i]=(1/DD[i]) }
                           DD = diag(DD)
                           smat  = t(sigma%*%comb%*%DD)

         #smat  = t(sigma%*%comb*as.numeric((1/(sqrt(t(comb[,1])%*%sigma%*%comb[,1])))))
 
                   };


    if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
    return(smat)
  }


#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S), using the Bo[,,,s] and Sigma[,,s] matrices of the estimated MRVAR.
#' @param res an MRVAR object as an output of MRVARest 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the weights used in GVAR models. For MRVAR its value is NA. 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol for Cholezky decomposition, chol1 for Cholezky decomposition with unit impulse, comb1 for concerted action with unit impulse.
#' @return an (n,n,nstep) array of IRF with colunms respesenting the impulse and rows the responses. 
#' @examples 
#'
#'
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1)
#' max(abs(res_d$Y))
#' res_e  = MRVARest(res_d)
#' summary_MRVAR(res_e)
#' 
#' IRF    = irf_MRVAR_NM1(res_e,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM1_CB2(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.90))
#' plott(IRF_CB[,,,,1],1,2)
#' plott(IRF_CB[,,,,2],1,2)
#' GIRF_migration    = girf_MRVAR_RM(RES=RESS,shock=c(1,1),R = 1800,nstep=10,state=1,resid_method="parametric")
#' GIRF_migration_cb = girf_MRVAR_RM_CB(RES=RESS,shock=c(1,1),R = 1800,nstep=10,state=1,resid_method="parametric",conf_level=c(0.1,0.9), N=200)
#' plott(GIRF_migration_cb,1,2)
#' @export
irf_MRVAR_NM1 = function(res,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1")) {
      p      = RESS$res$p
      Bo         = RESS$res$Bo
        Sigmao = RESS$res$Sigmao
      neq        = dim(Bo)[1]
        nvar     = dim(Bo)[2] 
         S   = dim(Bo)[4]   
      response <- array(0,dim=c(neq,nvar,nstep,S));
      BBo    = list()
      for (s in 1:S) { BBo[[s]] = Bo[,,1:p[s,1],s] ; if ( length(dim(BBo[[s]]))<3 ) dim(BBo[[s]]) = c(dim(BBo[[s]]),1)  }
      for (s in 1:S) response[,,,s] <- impulsdtrf_g(BBo[[s]],Sigmao[,,s],nstep,comb,irf=irf)
        return(response)
}
irf_MRVAR_NM = function(RESS,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1")) {
      Bo  = RESS$res$Bo
Sigmao = RESS$res$Sigmao
      neq  = dim(Bo)[1]
nvar = dim(Bo)[2] 
         S   = dim(Bo)[4]   
      response <- array(0,dim=c(neq,nvar,nstep,S));
      for (s in 1:S) response[,,,s] <- impulsdtrf_g(Bo[,,,s],Sigmao[,,s],nstep,comb,irf=irf)
return(response)
}


#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S) together with confidence bands.
#' Using the estimated Bo[,,,s] and Sigma[,,s] matrices of the estimated MRVAR, this function calculated the regime speicfic impulse response functions.
#' @param RESS a list of estimated MRVAR as output of MRVARest1 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a matrix of (n,n,nstep) as the IRF colunms respesenting the impulse and rows the responses. 
#' @examples 
#'
#' res_d = MRVARData4(n=2,p=5,T=300,S=2,SESVI=1)
#' RESS  = MRVARest1(res_d)
#' summary_MRVAR(RESS)
#' RESS$res$Bo
#' res_d$Bo
#' plot(irf(RESS$vars_obj[[1]]))  # only chol
#' #### IRF without regime migration 
#' IRF    = irf_MRVAR_NM(RESS,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM_CB(RESS,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.90))
#' plott(IRF_CB,1,2)
#'
#'
#' @export
irf_MRVAR_NM_CB = function(RESS,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95)) {
	n 	= RESS$res$n
      p 	= RESS$res$p
      T 	= dim(RESS$res$Y)[1]
      S 	= RESS$res$S
      SESVI = RESS$res$SESVI
      TH	= RESS$res$TH
	Bo	= RESS$res$Bo
      Co	= RESS$res$Co
      Sigmao= RESS$res$Sigmao
      SV 	= RESS$res$SV;   if (is.null(SV)) SV = NA 
      type 	= RESS$res$type
      X    	= RESS$res$X;   if (is.null(X))  X = NA  
      mu   	= RESS$res$mu;  if (is.null(mu)) mu = NA 
      Yo 	= RESS$res$Yo;  if (is.null(Yo)) Yo = NA 
      Do 	= RESS$res$Do;  if (is.null(Do)) Do = NA 
      Uo_run= RESS$res$Uo
   	neq 	 = dim(Bo)[1]
	nvar	 = dim(Bo)[2] 
      
      response <- array(0,dim=c(neq,nvar,nstep,length(conf)+1,S))
      response[,,,1,] <- irf_MRVAR_NM(RESS,nstep,comb,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs,S))
      for (i in 1:runs) {
      	for (s in 1:S) Uo_run[,,s] = rnormSIGMA(T,Sigmao[,,s]) 
            res_run = MRVARData4(n,p,T,S,SESVI,TH,Bo,Co,Sigmao,Uo_run,SV,type,X,mu,Yo,Do,d=1)
            res_e   = MRVARest1(res_run)
 		responseR[,,,i,] <- irf_MRVAR_NM(res_e,nstep,comb,irf)
	}

      for (tt in 1:(nstep) )    {
	for (i in 1:neq)           {
		for (j in 1:nvar)     {
			for (s in 1:S) response[i,j,tt,-1,s] = quantile(responseR[i,j,tt,,s], conf)
			
	} 
	}
	}
	return(response)
} 

#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S) together with confidence bands.
#' Using the estimated Bo[,,,s] and Sigma[,,s] matrices of the estimated MRVAR, this function calculated the regime speicfic impulse response functions.
#' @param RESS a list of estimated MRVAR as output of MRVARest1 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a matrix of (n,n,nstep) as the IRF colunms respesenting the impulse and rows the responses. 
#' @examples 
#'
#'
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' res_e = MRVARest(res=res_d)
#' summary_MRVAR(res_e)
#' RESS = res_e
#' #### IRF without regime migration
#' IRF = irf_MRVAR_NM1(RESS,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM1_CB2(RESS,nstep=20,comb=NA,irf="gen",runs=100,conf=c(0.05,0.90))
#' par(mfrow=c(4,2))
#' plott(IRF_CB[,,,,1],1,1)
#' plott(IRF_CB[,,,,1],1,2)
#' plott(IRF_CB[,,,,1],2,1)
#' plott(IRF_CB[,,,,1],2,2)
#' plott(IRF_CB[,,,,2],1,1)
#' plott(IRF_CB[,,,,2],1,2)
#' plott(IRF_CB[,,,,2],2,1)
#' plott(IRF_CB[,,,,2],2,2)
#'
#'
#' @export
irf_MRVAR_NM_CB1 = function(RESS,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95)) {
	n 	= RESS$res$n
      p 	= RESS$res$p
      T 	= dim(RESS$res$Y)[1]
      S 	= RESS$res$S
      SESVI = RESS$res$SESVI
      TH	= RESS$res$TH
	Bo	= RESS$res$Bo
      Co	= RESS$res$Co
      Sigmao= RESS$res$Sigmao
      SV 	= RESS$res$SV;   if (is.null(SV)) SV = NA 
      type 	= RESS$res$type
      X    	= RESS$res$X;   if (is.null(X))  X = NA  
      mu   	= RESS$res$mu;  if (is.null(mu)) mu = NA 
      Yo 	= RESS$res$Yo;  if (is.null(Yo)) Yo = NA 
      Do 	= RESS$res$Do;  if (is.null(Do)) Do = NA 
      Uo_run= RESS$res$Uo
   	neq 	 = dim(Bo)[1]
	nvar	 = dim(Bo)[2] 
      U_run  = Uo_run*0 
      response <- array(0,dim=c(neq,nvar,nstep,length(conf)+1,S))
      response[,,,1,] <- irf_MRVAR_NM1(RESS,nstep,comb,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs,S))
      for (i in 1:runs) {
      	for (s in 1:S) U_run[,,s] = rnormSIGMA(T,Sigmao[,,s]) 
            res_run = MRVARData(n,p,T,S,SESVI,TH,Bo,Co,Sigmao,U_run,SV,type,X,mu,Yo,Do,d=1)
            max(res_run$Y)
            res_e   = MRVARest(res_run)
 		responseR[,,,i,] <- irf_MRVAR_NM1(res_e,nstep,comb,irf)
            #if ( anyNA(responseR) ) break
	}

      for (tt in 1:(nstep) )    {
	for (i in 1:neq)           {
		for (j in 1:nvar)     {
			for (s in 1:S) response[i,j,tt,-1,s] = quantile(responseR[i,j,tt,,s], conf)
			
	} 
	}
	}
	return(response)
} 


#' Regime specific impulse response functions of MRVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRVAR(n,p,S) together with confidence bands, using the Bo[,,,s] and Sigma[,,s] matrices of the estimated MRVAR.
#' @param RESS a MRVAR object which is an output of MRVARest. 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol for Cholezky decomposition, chol1 for Cholezky decomposition with unit impulse, comb1 is used for global VAR.
#' @return an (n,n,nstep,3,2) array containing the IRF of the two regimes. The IRF colunms respesent the impulse and rows the responses. 
#' @examples 
#'
#'
#' p = matrix(c(2,1,0,0),2,2)
#' res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
#' max(res_d$Y)
#' res_e = MRVARest(res=res_d)
#' summary_MRVAR(res_e)
#' RESS = res_e
#' #### IRF without regime migration
#' IRF = irf_MRVAR_NM1(RESS,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_MRVAR_NM1_CB2(RESS,nstep=20,comb=NA,irf="gen",runs=100,conf=c(0.05,0.90))
#' par(mfrow=c(4,2))
#' plott(IRF_CB[,,,,1],1,1)
#' plott(IRF_CB[,,,,1],1,2)
#' plott(IRF_CB[,,,,1],2,1)
#' plott(IRF_CB[,,,,1],2,2)
#' plott(IRF_CB[,,,,2],1,1)
#' plott(IRF_CB[,,,,2],1,2)
#' plott(IRF_CB[,,,,2],2,1)
#' plott(IRF_CB[,,,,2],2,2)
#'
#' @export
irf_MRVAR_NM1_CB2 = function (RESS, nstep, comb, irf = c("gen", "chol", "chol1", 
    "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95)) 
{
    n = RESS$res$n
    p = RESS$res$p
    T = dim(RESS$res$Y)[1]
    S = RESS$res$S
    SESVI = RESS$res$SESVI
    TH = RESS$res$TH
    Bo = RESS$res$Bo
    Co = RESS$res$Co
    Sigmao = RESS$res$Sigmao
    SV = RESS$res$SV
    if (is.null(SV)) 
        SV = NA
    type = RESS$res$type
    X = RESS$res$X
    if (is.null(X)) 
        X = NA
    mu = RESS$res$mu
    if (is.null(mu)) 
        mu = NA
    Yo = RESS$res$Yo
    if (is.null(Yo)) 
        Yo = NA
    Do = RESS$res$Do
    if (is.null(Do)) 
        Do = NA
    Uo_run = RESS$res$Uo
    neq = dim(Bo)[1]
    nvar = dim(Bo)[2]
    U_run = Uo_run * 0
    response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 
        1, S))
    response[, , , 1, ] <- irf_MRVAR_NM1(RESS, nstep, comb, irf)
    responseR <- array(0, dim = c(neq, nvar, nstep, runs, S))
    for (i in 1:runs) {
        for (s in 1:S) {
            U_run[, , s] = rnormSIGMA(T, Sigmao[, , s])
            B = Bo[, , 1:p[s, 1], s]
            if (length(dim(B)) == 2) 
                dim(B) = c(dim(B), 1)
            res_run = VARData(n = n, p = p[s, 1], T = T, B = B, 
                Co = 0 * (1:n), Sigma = Sigmao[, , s], U = U_run[, 
                  , s], type = "none")
            res_e = VARest(res_run)
            responseR[, , , i, s] <- impulsdtrf_g(res_e$B, res_e$Sigma, 
                nstep, comb, irf)
        }
    }
    for (tt in 1:(nstep)) {
        for (i in 1:neq) {
            for (j in 1:nvar) {
                for (s in 1:S) response[i, j, tt, -1, s] = quantile(responseR[i, 
                  j, tt, , s], conf)
            }
        }
    }
    return(response)
}




#' Calculation of information criteria for a MRVAR model.  
#' 
#' @param  res  : an object generated from MRVARData or estimated from MRVARest.
#' @param  L_V  : a 2 components vector containing the maxima of the domestic lag and the foreign lag respectively. 
#' @param  TH_V  : a vector containing the possible thrshold values over which the model selection criteria valuee will be calculated respectively. 

#' @return      : A matrix with different lag specifications and threshold values as well as the corresponding model selection criterion values.
#' @examples 
#' 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)*2
#' Sigma[,,2] = diag(4)*2
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' 
#' res_d = MRVARData(n=4,p=p,T=8000,S=2,SESVI=1,Sigmao=Sigma,type="none")
#' max(res_d$Y)
#' res_e = MRVARest(res=res_d)
#' TH_v = c(0,0.0)
#' L_v = c(5,5)
#' Sel = MRVAR_Select(res=res_e,L_V=L_v,TH_V=TH_v)
#' Sel[which.min(Sel[,5]),]
#' @export
MRVAR_Select = function(res,L_V=L_v,TH_V=TH_v) {
res_dd = res$res
p      = res_dd$p
n 	 = res_dd$n
T    	 = dim(res_dd$Y)[1]
SESVI  = res_dd$SESVI
type   = res_dd$type
X      = res_dd$X
p[,1]  = max(L_V)
S      = res_dd$S
res_dd = MRVARData(n=n,p=p,T=T,S=S,SESVI=SESVI,type=type,X=X)
res_dd$Y= res$res$Y

Criteria = matrix(0,L_V[1]*L_V[2]*length(TH_v),13)
idx = 0

for (l_d in 1: L_V[1] )             {
   for (l_f in 1:L_V[2] )           {
      for (l_th in 1:length(TH_v))  {
      idx = idx + 1
      res_dd$p[1,1] = l_d
      res_dd$p[2,1] = l_f
      res_dd$TH     = TH_V[l_th]
	res_ss = MRVARest(res_dd)
      TT = res_ss$LH_N + res_ss$LH_P
	Criteria[idx,] = c(l_d,l_f,TH_V[l_th],res_ss$LH_AIC,res_ss$LH_BIC,res_ss$LH_P,res_ss$LH_N,res_ss$LH_NN,res_ss$ORAIC,res_ss$ORBIC,TT,res_ss$LLH_AIC,res_ss$LLH_BIC)
      colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LH_P","LH_N","LH_NN","ORAIC","ORBIC","TT","LLH_AIC","LLH_BIC") 
    }
}
}
return(Criteria)
}



#' @export
MRVAR_Select_Summary = function(Sel) {
  result = list()
  result[[1]] = Sel[which.min(Sel[,5]),c(1,2,3,5)]
  result[[2]] = Sel[which.min(Sel[,4]),c(1,2,3,4)]
  return(result)
}