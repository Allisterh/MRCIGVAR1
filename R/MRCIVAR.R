#' Data generating process of MRCIVAR(n,p,S) 
#'
#' This function will generate data from an multi-regime cointegrated VAR(p) process with identical cointegration relations and identical adjustment speeds and return an object of MRCIVAR containing data and parameters used in the multi regime cointegrated VAR(p) process.
#'
#' @param n     : number of variables
#' @param S     : number of regimes 
#' @param p     : an S x 2 matrix. Each row contains the lag length for the corresponding state and the number of exogenous variables for the state.
#' @param T     : number of observations
#' @param SESVI : Index of the switching variable, switching sv_t = Y[tt-1,SESVI]>Y[tt-2,SESVI] 
#' @param TH    : S-1 vector of threshold values 
#' @param Bo    : (n,n,p,S) array of the coeffiicent matrices of VAR(p) for each lag and each regime. If not given it will be generated.
#' @param Co    : (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an n vector for each regime, exog0 Co is a (n,k+1,S) array with first colume zero for each regime respectively,for exog1 Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao: (n,n,S)   array containing S covariance matrices of residuals, each for one regime.
#' @param Uo    : residuals, if it is not NA it will be used as input to generate the CIVAR(p) for each regime.
#' @param SV    : exogeous switching variable
#' @param type  : type of the deterministic components type = ("none","const","exog0","exog1")   
#' @param X     : (T x k x S) matrix of exogeneous variables for each state. The second dimension can be filled with zeros to take into account that the the exogenous variables are not identical in each state. 
#' @param mu    : (n x S) matrix of the regime specific mean of the variables
#' @param Yo    : (p, n, S) array of initial values of the process
#' @param Do    : (T, n, S) array of extra exogenous components (not used with value zero) 
#' @param d     : lag of the self-exiting switching variable  
#' @param r     : the number of unit roots in the system. n-r is then the cointegration rank.    
#' 
#'               (TH,Bo,Sigmao,Uo) if not provided, they will be generated randomly.
#' @return      A MRCIVAR object containing the generated data, the parameters and the input exogeous variables. res = list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d,r,crk)
#' @examples 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' res_d$p
#' res_d = MRCIVARData(n=4,p=p,T=300,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="none",r=1)
#' res_e = MRCIVARest(res=res_d)
#' sum(abs(res_d$Bo - res_e[[1]]$Bo))
#' sum(abs(res_d$Co - res_e[[1]]$Co))
#' res_d$p
#' res_d = MRCIVARData(n=4,p=p,T=300,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
#' res_e = MRCIVARest(res=res_d)
#' sum(abs(res_d$Bo - res_e[[1]]$Bo))
#' sum(abs(res_d$Co - res_e[[1]]$Co))
#' @export
MRCIVARData = function(n=2,p=matrix(2,2,2),T=100,S=2,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d=1,r=1) {
### n     : number of variables
### p     : lag length, an S x 2 vector of lag length for each state and the number of exogenous variables for each state
### T     : number of observations
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
  Pmin = min(p[,1])
  if (anyNA(TH)) {
      TH = matrix(runif(S-1),S-1,1)
      TH = TH[order(TH)]
  }
  if (anyNA(Bo)) {
  	Bo = (1:(n*n*Pmax*S))*0
      dim(Bo) = c(n,n,Pmax,S)
      r_npo = t(matrix( rep((2:(n+1)),Pmax),Pmax,n))
      for (i in 1:r ) r_npo[i,1] = 1
	res_d  = VARData(n,Pmax,T=100,r_np=r_npo)
		#plot(ts(res_d$Y))
	B1 = res_d$B
	CIB1 = B2CIB(B1)[[1]]
		#STAT(B1)
	repeat {
		res_d2  = VARData(n,Pmin,T=100)
      	B2 = B1*0
		B2[,,1:Pmin] = res_d2$B
		CIB2 = B2CIB(B2)[[1]]
		#CIBB = B1*0
		#CIBB[,,1] = CIB1[,,1]
		#for (i in 2:Pmin ) CIBB[,,i] = CIB2[,,i] 
		#CBB = CIB3B(CIBB)
		#STAT(CBB)
            CIB2[,,1]  = CIB1[,,1]
            CBB = CIB3B(CIB2)
		if (abs(STAT(CBB)[2])<1) break
	}
      
      if (p[1,1] > p[2,1]) { 
      	Bo[,,,1] = B1 ; Bo[,,,2] = CBB      
      }     else           {
            Bo[,,,1] = CBB ; Bo[,,,2] = B1    
      }

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

  if (anyNA(type)|type=="none") {
  	type = "none"
      Co = (1:(S*n))*0
      dim(Co) = c(n,1,S)
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
      if (anyNA(SV)) {sv = Y[tt-1,SESVI]>Y[tt-2,SESVI]} else {sv = SV[tt]}
      st = sv;
      s  = sum(st>TH)+1
      St[tt] = s
      Y[tt,] = Uo[tt,,s]+Ct[tt,,s] + Do[,s]
      for (L in 1:Pmax)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bo[,,L,s])
  }
  check[S+1] = max(abs(Y))
  resid = NA
  crk   = n - r
  result=list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d,r,crk)
  names(result) = c("Y","X","Uo","resid","Bo","Co","Sigmao","TH","St","sv","d","SESVI","r_npo","check","n","p","S","type","Yo","d","r","crk")
  return(result)
}




#' Data generating process of MRCIVAR(n,p,S) 
#'
#' This function will generate data from an multi-regime cointegrated VAR(p) process with identical cointegration relations and different adjustment speeds and return an object of MRCIVAR containing data and parameters used in the multi regime cointegrated VAR(p) process.
#'
#' @param n     : number of variables
#' @param S     : number of regimes 
#' @param p     : an S x 2 matrix. Each row contains the lag length for the corresponding state and the number of exogenous variables for the state.
#' @param T     : number of observations
#' @param SESVI : Index of the switching variable, switching sv_t = Y[tt-1,SESVI]>Y[tt-2,SESVI] 
#' @param TH    : S-1 vector of threshold values 
#' @param Bo    : (n,n,p,S) array of the coeffiicent matrices of VAR(p) for each lag and each regime. If not given it will be generated.
#' @param Co    : (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an n vector for each regime, exog0 Co is a (n,k+1,S) array with first colume zero for each regime respectively,for exog1 Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao: (n,n,S)   array containing S covariance matrices of residuals
#' @param Uo    : residuals, if it is not NA it will be used as input to generate the VAR(p) for each regime.
#' @param SV    : exogeous switching variable
#' @param type  : type of the deterministic components type = ("none","const","exog0","exog1")   
#' @param X     : (T x k x S) matrix of exogeneous variables for each state. The second dimension can be filled with zeros to take into account that the the exogenous variables are not identical in each state. 
#' @param mu    : (n x S) matrix of the regime specific mean of the variables
#' @param Yo    : (p, n, S) array of initial values of the process
#' @param Do    : (T, n, S) array of extra exogenous components (not used with value zero) 
#' @param d     : lag of the self-exiting   
#' @param r     : the number of unit roots in the system. n-r is then the cointegration rank.    
#' 
#'               (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
#' @return      A MRCIVAR object containing the generated data, the parameters and the input exogeous variables. res = list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d,r,crk)
#' @examples 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' res_d$p
#' res_d = MRCIVARDatam(n=4,p=p,T=300,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="none",r=1)
#' res_e = MRCIVARestm(res=res_d)
#' sum(abs(res_d$Bo - res_e[[1]]$Bo))
#' sum(abs(res_d$Co - res_e[[1]]$Co))
#' res_d$p
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y))
#' res_e$res$tst[[1]]
#' STAT(res_e$res$B[,,,1])
#' STAT(res_e$res$B[,,,2])
#'
#' @export
MRCIVARDatam = function(n=2,p=matrix(2,2,2),T=100,S=2,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d=1,r=1) {
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
  Pmin = min(p[,1])
  if (anyNA(TH)) {
      TH = matrix(runif(S-1),S-1,1)
      TH = TH[order(TH)]
  }
  if (anyNA(Bo)) {
  Bo = (1:(n*n*Pmax*S))*0
      dim(Bo) = c(n,n,Pmax,S)
      r_npo = t(matrix( rep((2:(n+1)),Pmax),Pmax,n))
      for (i in 1:r ) r_npo[i,1] = 1
res_d  = VARData(n,Pmax,T=100,r_np=r_npo)
#plot(ts(res_d$Y))
B1 = res_d$B
CIB1 = B2CIB(B1)[[1]]
#STAT(B1)
repeat {
res_d2  = VARData(n,Pmin,T=100)
      B2 = B1*0
B2[,,1:Pmin] = res_d2$B
CIB2 = B2CIB(B2)[[1]]
#CIBB = B1*0
#CIBB[,,1] = CIB1[,,1]
#for (i in 2:Pmin ) CIBB[,,i] = CIB2[,,i] 
#CBB = CIB3B(CIBB)
#STAT(CBB)
            values = eigen(CIB1[,,1])$values
            VECTOR = eigen(CIB1[,,1])$vectors

            if ((n-r) ==1 )   CIB2[,,1]  = VECTOR[,(r+1):n]%*%as.matrix(-runif(n-r))%*%solve(VECTOR)[(r+1):n,] 
            if ((n-r) > 1 )   CIB2[,,1]  = VECTOR[,(r+1):n]%*%(diag(-runif(n-r)))%*%solve(VECTOR)[(r+1):n,]   
            #if ((n-r) ==1 )   CIB2[,,1]  = VECTOR[,1:(n-r)]%*%as.matrix(-runif(n-r))%*%solve(VECTOR)[(r+1):n,] 
            #if ((n-r) > 1 )   CIB2[,,1]  = VECTOR[,1:(n-r)]%*%(diag(-runif(n-r)))%*%solve(VECTOR)[(r+1):n,]   


            CBB = CIB3B(CIB2)

if (max(abs(STAT(CBB)))==1) break
}
      
      if (p[1,1] > p[2,1]) { 
      Bo[,,,1] = B1 ; Bo[,,,2] = CBB      
      }     else           {
            Bo[,,,1] = CBB ; Bo[,,,2] = B1    
      }

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

  if (anyNA(type)|type=="none") {
  type = "none"
      Co = (1:(S*n))*0
      dim(Co) = c(n,1,S)
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
  }  else {
            mu = Co*NA

            for (s in 1:S) {
H = diag(n)
for (L in 1:Pmax) H = H - Bo[,,L,s]
            #mu[,1,s] = solve(H)%*%Co[,1,s]
            mu = NA
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
      if (anyNA(SV)) {sv = Y[tt-1,SESVI]>Y[tt-2,SESVI]} else {sv = SV[tt]}
      st = sv;
      s  = sum(st>TH)+1
      St[tt] = s
      Y[tt,] = Uo[tt,,s]+Ct[tt,,s] + Do[,s]
      for (L in 1:Pmax)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bo[,,L,s])
  }
  check[S+1] = max(abs(Y))
  resid = NA
  crk   = n - r
  result=list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d,r,crk)
  names(result) = c("Y","X","Uo","resid","Bo","Co","Sigmao","TH","St","sv","d","SESVI","r_npo","check","n","p","S","type","Yo","d","r","crk")
  return(result)
}



#' Data generating process of MRCIVAR(n,p,S) 
#'
#' This function will generate data from an multi-regime cointegrated VAR(p) process with identical cointegration relations and different adjustment speeds and return a list containing data and parameters used in the multi regime cointegrated VAR(p) process. The function assumes a underlying processes consisting of a set of I(1) processes and two sets of I(0) processes for the two different regimes.  The data generated are obtained through a linear transformation of the underlying processes.  
#' @param n     : number of variables
#' @param S     : number of regimes 
#' @param p     : an S x 2 matrix. Each row contains the lag length for the corresponding state and the number of exogenous variables for the state.
#' @param T     : number of observations
#' @param SESVI : Index of the switching variable, switching sv_t = Y[tt-1,SESVI]>Y[tt-2,SESVI] 
#' @param TH    : S-1 vector of threshold values 
#' @param Bo    : (n,n,p,S) array of the coeffiicent matrices of VAR(p) for each lag and each regime. If not given it will be generated.
#' @param Co    : (n,k+1,S) array of the coefficients of the deterministic components. For type="none" Co = O*(1:n,1:S), for "const" Co is an n vector for each regime, exog0 Co is a (n,k+1,S) array with first colume zero for each regime respectively,for exog1 Co is an (n,1+k, S) array without zero restrictions.
#' @param Sigmao: (n,n,S)   array containing S covariance matrices of residuals
#' @param Uo    : residuals, if it is not NA it will be used as input to generate the VAR(p) for each regime.
#' @param SV    : exogeous switching variable
#' @param type  : type of the deterministic components type = ("none","const","exog0","exog1")   
#' @param X     : (T x k x S) matrix of exogeneous variables for each state. The second dimension can be filled with zeros to take into account that the the exogenous variables are not identical in each state. 
#' @param mu    : (n x S) matrix of the regime specific mean of the variables
#' @param Yo    : (p, n, S) array of initial values of the process
#' @param Do    : (T, n, S) array of extra exogenous components (not used with value zero) 
#' @param d     : lag of the self-exiting   
#' @param r     : the number of unit roots in the system. n-r is then the cointegration rank.    
#' 
#'               (TH,Bo,Sigmao,Uo,SV) if not provided, they will be generated randomly.
#' @return      A list containing the generated data, the parameters and the input exogeous variables. res = list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d,r,crk)
#' @examples 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' res_d$p
#' res_d = MRCIVARDatam1(n=4,p=p,T=300,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="none",r=2)
#' res_e = MRCIVARestm(res=res_d)
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y))
#' res_e$res$tst[[1]]
#' STAT(res_d$B[,,,1])
#' STAT(res_d$B[,,,2])
#' @export
MRCIVARDatam1 = function(n=2,p=matrix(2,2,2),T=100,S=2,SESVI,TH,Bo,Co,Sigmao,Uo,SV,type,X,mu,Yo,Do,d=1,r=1) {
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
### r     : number of unit roots
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
### r_npo : n x p x S array collecting the roots of the characteristic functions in L of the n dynamically independent processes.
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
  Pmin = min(p[,1])
  if (anyNA(TH)) {
      TH = matrix(runif(S-1),S-1,1)
      TH = TH[order(TH)]
  }
  r_npo = NA



 if (anyNA(Bo)) {
 	### unit root process with Pmin lag regime
 	r_np = c(1:(r*Pmin))/c(1:(r*Pmin))*1.5
      dim(r_np) = c(r,Pmin)
      r_np[,1]  = 1     
      VAR_IONE  = VARData(n=r,p=Pmin,T=100,r_np=r_np)
      VAR_IONE$B
      ### stationary process with lag p
      r_np = c(1:((n-r)*p[1,1]))/c(1:((n-r)*p[1,1]))*2
      dim(r_np)   = c(n-r,p[1,1])
      VAR_IZeroI  = VARData(n-r, p[1,1], T=100,r_np=r_np)
      plot(ts(VAR_IZeroI$Y))
      VAR_IZeroI$B
      r_np = c(1:((n-r)*p[2,1]))/c(1:((n-r)*p[2,1]))*1.5
      dim(r_np)   = c(n-r,p[2,1])
      VAR_IZeroII  = VARData(n-r, p[2,1], T=100,r_np=r_np)
      
      ### B is the transformation matrix to mix the block diagonal the I(1) and I(0) process
      B = diag(n)*0.5+matrix(rnorm(n*n),n,n)
  
       Bo = (1:(n*n*Pmax*S))*0
       dim(Bo) = c(n,n,Pmax,S)
 
       for (j in 1: p[1,1]) {
               
               Dblock    = matrix(0,n,n)
               Dblock[1:(n-r),1:(n-r)]     = VAR_IZeroI$B[,,j]
               if (dim(VAR_IONE$B)[3]>=j) Dblock[(n-r+1):n,(n-r+1):n] = VAR_IONE$B[,,j]
               Bo[,,j,1] = B%*%Dblock%*%solve(B)       
	}

      for (j in 1: p[2,1]) {
               Dblock    = matrix(0,n,n)
               Dblock[1:(n-r),1:(n-r)]   = VAR_IZeroII$B[,,j]
               if (dim(VAR_IONE$B)[3]>=j) Dblock[(n-r+1):n,(n-r+1):n] = VAR_IONE$B[,,j]
               Bo[,,j,2] = B%*%Dblock%*%solve(B)       
	}
      
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

  if (anyNA(type)|type=="none") {
  type = "none"
      Co = (1:(S*n))*0
      dim(Co) = c(n,1,S)
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
  }  else {
            mu = Co*NA

            for (s in 1:S) {
H = diag(n)
for (L in 1:Pmax) H = H - Bo[,,L,s]
            #mu[,1,s] = solve(H)%*%Co[,1,s]
            mu = NA
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
      if (anyNA(SV)) {sv = Y[tt-1,SESVI]>Y[tt-2,SESVI]} else {sv = SV[tt]}
      st = sv;
      s  = sum(st>TH)+1
      St[tt] = s
      Y[tt,] = Uo[tt,,s]+Ct[tt,,s] + Do[,s]
      for (L in 1:Pmax)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bo[,,L,s])
  }
  check[S+1] = max(abs(Y))
  resid = NA
  crk   = n - r
  result=list(Y,X,Uo,resid,Bo,Co,Sigmao,TH,St,sv,d,SESVI,r_npo,check,n,p,S,type,Yo,d,r,crk)
  names(result) = c("Y","X","Uo","resid","Bo","Co","Sigmao","TH","St","sv","d","SESVI","r_npo","check","n","p","S","type","Yo","d","r","crk")
  return(result)
}






#' Estimation of MRCIVARm(n,p,S) 
#'
#' This function estimate the unknown parameters with different adjustment speed of a specified MRVAR(n,p,S) model based on provided data.
#'
#' @param  res  : a list containing the components as the output of MRVARData4 including at least: n, p, type, Y, SESVI, TH, d, and optionally X. 
#' @return 	: c("res","LH_AIC","LH_BIC","LH_P","LH_N","LHH_AIC","LHH_BIC","ORBIC","ORAIC","ORLH","ORLHN","vars_obj","St","TT","TTT","LH_TN","LH_TN1","n2LHP") res is a list as the input.
#' @examples 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' res_d$p
#' res_d = MRCIVARDatam1(n=4,p=p,T=300,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="none",r=2)
#' res_e = MRCIVARestm(res=res_d)
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y))
#' res_e$res$tst[[1]]
#' STAT(res_d$B[,,,1])
#' STAT(res_d$B[,,,2])
#' summary(res_e$vars_obj)
#'
#' @export
MRCIVARestm = function(res) {
##
## This is a program used to estiamte MRVAR parameters under given order parameter of MRVAR
##
TH = res$TH
type = res$type
p     = res$p   # S x 2 vector of lag lengths for each state, the second column is the the number of exogenous variables
S = res$S
n = res$n
Y = res$Y
X = res$X
T = dim(Y)[1]
SV    = res$SV
sv    = res$sv
d     = res$d
SESVI = res$SESVI
St    = 2-res$St
P     = max(p[,1],d)
Pmax  = max(p[,1])
Co    = res$Co
r     = res$r
crk   = res$crk
ms    = (1:S)*0
LH_AIC = 0
LH_BIC = 0
LH_P   = 0
LH_N   = 0

b     = res$Bo*NA
cc    = res$Co*NA
sigma = res$Sigmao*NA

resid  = (1:(T*n*S))*0; dim(resid) = c(T,n,S) 
Tresid = resid[,,1]
#### estimation of one regime CIVAR
    ORCIVARD   = CIVARData(n,p=Pmax,T=T,Co=Co[,,1],type=type)   
    ORCIVARD$Y = Y
    ORCIVAR_e  = CIVARest(res=ORCIVARD)

#### estimation of multi regime CIVAR
    if (type=="none")  Model = "I"
    if (type=="const") Model = "III"
    tst <- MRCVECMestm(y=Y,x=0,s=St,model = Model, type = "eigen",P = p , crk = n-r , q = 0.95)


    Tresid[(T-dim(tst[[2]]$residuals)[1]+1):T,]=tst[[2]]$residuals
    
    Sigma_one = t(Tresid)%*%Tresid/(T-dim(tst[[2]][[1]])[1])
    resid[,,1] = Tresid*St
    resid[,,2] = Tresid*(1-St)
    sigma[,,1] = t(resid[,,1])%*%resid[,,1]/(sum(St)-1*(n*(p[2,1]-1)+crk)) 
    sigma[,,2] = t(resid[,,2])%*%resid[,,2]/(sum(1-St)-1*(n*(p[2,1]-1)+crk)) 
    CIVAREST = VECM2VARm(param=tst$VECM1[[1]],beta=tst$betaS,p=c(crk,p[1,1]-1,p[2,1]-1),s=1)
    Bo = CIVAREST[[1]]
    if (type=="const")  Co = CIVAREST[[3]] 

    #LH_P     = -(T*n/2)*log(2*pi) -(T*n/2) + (T/2)*log(det(solve(Sigma_one)))
    LH_P      = -(T*n/2)*log(2*pi) -(T*n/2) + (sum(St))/2*log(det(solve(sigma[,,1])))+(sum(1-St))/2*log(det(solve(sigma[,,2])))

    LH_P      = -(T*n/2)*log(2*pi) -(T*n/2) + (sum(St))/2*log(det(solve(sigma[,,1])))+(sum(1-St))/2*log(det(solve(sigma[,,2])))



    LH_AIC   =           2*(n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2) - 2*LH_P
    LHH_AIC  =           2*(n*(n*(p[1,1]-1+crk/2)+n*(n+1)/2)+n*(n*(p[2,1]-1+crk/2)+n*(n+1)/2))  - 2*LH_P

    LH_BIC   =      log(T)*(n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2) - 2*LH_P
    LHH_BIC  =            log(sum(St))*n*(n*(p[1,1]-1+crk/2)+n*(n+1)/2)+ log(sum(1-St))*n*(n*(p[2,1]-1+crk/2)+n*(n+1)/2)  - 2*LH_P
    TT       =     (p[1,1]-1)*n+(p[2,1]-1)*n + crk - dim(tst[[2]][[1]])[1]
    TTT      =      dim(tst[[2]][[1]])[1]

    LH_N     =     2*n*(dim(tst[[2]][[1]])[1])+n*(n+1)
    LH_TN    =     log(sum(St))*n*(n*(p[1,1]-1+crk/2)+n*(n+1)/2)+ log(sum(1-St))*n*(n*(p[2,1]-1+crk/2)+n*(n+1)/2)
    LH_TN1   =     log(T)*(n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2)
    lm_obj   = tst$VECM1


    ORBIC    =  ORCIVAR_e$BIC 
    ORAIC    =  ORCIVAR_e$AIC
    ORLH     =  ORCIVAR_e$LH
    ORLHN    =  ORCIVAR_e$LH_N   
    res$Bo      <- Bo
    res$Co      <- Co
    res$Sigmao  <- sigma
    res$St      <- 2-St
    res$resid   <- resid
    res$tst     <- tst
result = list(res,LH_AIC,LH_BIC,LH_P,LH_N,LHH_AIC,LHH_BIC,ORBIC,ORAIC,ORLH,ORLHN,lm_obj,St,TT,TTT,LH_TN, LH_TN1,-2*LH_P)
names(result) = c("res","LH_AIC","LH_BIC","LH_P","LH_N","LHH_AIC","LHH_BIC","ORBIC","ORAIC","ORLH","ORLHN","vars_obj","St","TT","TTT","LH_TN","LH_TN1","n2LHP")
return(result)
}


#' Estimation of MRCIVAR(n,p,S) 
#'
#' This function estimates the unknown parameters with identical adjustment speed of a specified MRVAR(n,p,S) model based on provided data.
#
#'
#' @param  res  : an MRCIVAR object as the output of MRCIVARData including at least: n, S, p, type, Y, SESVI, TH, d, and optionally X. 
#' @return 	: list("res","LH_AIC","LH_BIC","LH_P","LH_N","ORBIC","ORAIC","vars_obj") res is an object of MRCIGVAR.
#' @examples 
#'
#'
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' 
#' res_d = MRCIVARData(n=4,p=p,T=261,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=2)
#' 
#' res_e = MRCIVARest(res=res_d)
#' res_e[[1]]$tst[[1]]
#' summary(res_e$vars_obj)
#' res_em = MRCIVARestm(res=res_d)
#' res_em[[1]]$tst[[1]]
#' summary(res_em$vars_obj)
#'
#' @export
MRCIVARest = function(res) {
##
## This is a program used to estiamte MRVAR parameters under given order parameter of MRVAR
##
TH = res$TH
type = res$type
p     = res$p   # S x 2 vector of lag lengths for each state, the second column is the the number of exogenous variables
S= res$S
n= res$n
Y= res$Y
X = res$X
T= dim(Y)[1]
SV    = res$SV
sv    = res$sv
d     = res$d
SESVI = res$SESVI
St    = 2-res$St
P     = max(p[,1],d)
Pmax  = max(p[,1])
Co    = res$Co
r     = res$r
crk   = res$crk
ms    = (1:S)*0
LH_AIC = 0
LH_BIC = 0
LH_P   = 0
LH_N   = 0

b     = res$Bo*NA
cc    = res$Co*NA
sigma = res$Sigmao*NA

resid  = (1:(T*n*S))*0; dim(resid) = c(T,n,S) 
Tresid = resid[,,1]
#### estimation of one regime CIVAR
    ORCIVARD   = CIVARData(n,p=Pmax,T=T,Co=Co[,,1],type=type)   
    ORCIVARD$Y = Y
    ORCIVAR_e  = CIVARest(res=ORCIVARD)

#### estimation of multi regime CIVAR
    if (type=="none")  Model = "I"
    if (type=="const") Model = "III"
    tst <- MRCVECMest2(y=Y,x=0,s=St,model = Model, type = "eigen",P = p , crk = n-r , q = 0.95)


    Tresid[(T-dim(tst[[2]]$residuals)[1]+1):T,]=tst[[2]]$residuals
    Sigma_one = t(Tresid)%*%Tresid/(T-dim(tst[[2]][[1]])[1])
    resid[,,1] = Tresid*St
    resid[,,2] = Tresid*(1-St)
    sigma[,,1] = t(resid[,,1])%*%resid[,,1]/(sum(St)-(crk+p[1,1])*n) 
    sigma[,,2] = t(resid[,,2])%*%resid[,,2]/(sum(1-St)-(crk+p[2,1])*n) 
    CIVAREST = VECM2VAR(param=tst[[2]][[1]],beta=tst$beta,p=c(crk,p[1,1]-1,p[2,1]-1),s=1)
    Bo = CIVAREST[[1]]
    if (type=="const")  Co = CIVAREST[[3]] 

    LH_P     = -(T*n/2)*log(2*pi) -(T*n/2) +(T/2)*log(det(solve(Sigma_one)))
    LH_AIC   = 2*n*(dim(tst[[2]][[1]])[1])+n*(n+1) - 2*LH_P
    LH_BIC   = log(T)*(n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2) - 2*LH_P
    LH_N     = (n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2)
    lm_obj   = tst[[2]]


    ORBIC    =  ORCIVAR_e$BIC 
    ORAIC    =  ORCIVAR_e$AIC
   
  

    res$Bo      <- Bo
    res$Co      <- Co
    res$Sigmao  <- sigma
    res$St      <- St
    res$resid   <- resid
    res$tst     <- tst
result = list(res,LH_AIC,LH_BIC,LH_P,LH_N,ORBIC,ORAIC,lm_obj)
names(result) = c("res","LH_AIC","LH_BIC","LH_P","LH_N","ORBIC","ORAIC","vars_obj")
return(result)
}





#' Estimation of a two regime conditional vector error correction process.
#'
#' This function estimates the unknown parameters of a conditional VECM based on provided data.
#'
#' @param y	: data matrix of the endogenous variable 
#' @param x	: x data matrix of the conditioning variables. If x is missing, it will estimate a unconditional VECM.
#' @param s     : the series of state variable of value (0,1) for the two regime case. Missing s implies single regime VECM.  
#' @param model : It assumes one of the values in c("I","II","III","IV","V") corresponding to the five cases discussed in JH book. 
#' @param type  : c("eigen", "trace") corresponding to the trace test or the max eigenvalue test
#' @param crk   : cointegration rank
#' @param P     : 2 x 2 matrix containing the lag of the cointegrated VAR process. The first row is the lag of the first regime and the number of exogenous variables. If the second row is zero, this is a one regime VECM.  
#' @return      : a list containing: the result of JH test, VECM in regression format, lambda, beta, PI, GAMMA, model, P, s 
#' @examples 
#' @export
MRCVECMestm = function (y,x,s,model = c("I","II","III","IV","V"),type = c("eigen", "trace"), constant = TRUE, ret = c("statistic", 
    "test"), ctable = c("A3", "A1", "A5"), crk=2, P =matrix(2,2,2), q = 0.95,Dxflag=0) 
{
    y <- as.matrix(y)
    model <- match.arg(model)
    type <- match.arg(type)
    ret <- match.arg(ret)
    ctable <- match.arg(ctable)
    if (q != 0.9 && q != 0.95 && q != 0.99) {
        print("please correct significance level")
        break()
    }

    S = 2
    p <- as.integer(max(P))
    pmin <- as.integer(min(P))
    N1 <- ncol(as.matrix(y))
    NN1 <-ncol(as.matrix(x)) 
    N  <- ncol(as.matrix(y)) 
    n = N
    if (N1<crk) {
        print("y's dimension must be larger than crk")
        break()
    }
    if (missing(s)) {
      s = NA
    }
    if (missing(Dxflag)) {
      Dxflag = 0
    }

    if (!anyNA(s)) {
	St  = s[ (p+1):length(s)] 
      NSt = 1-s[(p+1):length(s)]
    } 
     
 
    
    if (!sum(abs(x))==0)  {
      z  <-  cbind(y,x)    
      Zy <-  embed(diff(y),p)
      Zx <-  embed(diff(x),p)
      ## taking out the non-appearing lags    

      if (P[1,1]<p) {
          Aa = P[1,1]*N1+(1:((p-P[1,1])*N1))    	
          Zy1 = Zy[,-Aa]
      }   else  Zy1 = Zy

      if (P[1,2]<p) {
          Bb = P[1,2]*NN1+(1:((p-P[1,2])*NN1))    	
          Zx1 = Zx[,-Bb]
      }   else Zx1 = Zx

      if (P[2,1]<p) {
          Aa = P[2,1]*N1+(1:((p-P[2,1])*N1))    	
          Zy2 = Zy[,-Aa]
      }   else Zy2 = Zy
      if (P[2,2]<p) {
          Bb = P[2,2]*NN1+(1:((p-P[2,2])*NN1))    	
          Zx2 = Zx[,-Bb]
      }   else Zx2 = Zx

      Z = cbind( Zy1,Zx1); 
      Z_2 = cbind(Zy2,Zx2)
      if (Dxflag == 0) {
		Z2    <-   Z[, -c(1:N1,P[1,1]*N1+1:NN1)] 
		ZS_2  <- Z_2[, -c(1:N1,P[2,1]*N1+1:NN1)] 
		}   else  {
      	Z2    <-   Z[, -c(1:N1)] 
      	ZS_2  <- Z_2[, -c(1:N1)] 
	}   
   
    }  else  {   
      z = y
      Z = embed(diff(y),p)
	if (Dxflag == 0) {
		Z2    <-   Z[, -c(1:N1)]
            ZS_2  <- Z2
            Z2    <- Z2[,1:((P[1,1]-1)*n)]
            ZS_2  <- ZS_2[,1:((P[2,1]-1)*n)]
	  	}   else  {
      	Z2    <-   Z[, -c(1:N1)]
            Zs_2  <- Z2 
            Z2    <- Z2[,1:((P[1,1]-1)*n)]
            ZS_2  <- ZS_2[,1:((P[2,1]-1)*n)]
      }     	

    }

    M1 <- ncol(as.matrix(z))
    T <-  nrow(as.matrix(z))

    MM1 = ncol(Z);
    Y0 <- Z[, c(1:N1)]                # \Delta Y
    #X0 <- Z[,c((N1+1):M1)]           # \Delta X          
    Z1 <- z[-T, ][p:(T - 1), ]        #      Z_1

												    # \Delta Z_    
 
    if (!anyNA(s))  Z2 = cbind(St*Z2,NSt*ZS_2)                                          # two regimes with common CI space
      
    T1 <- nrow(as.matrix(Y0))
    lT = (1:T1)/(1:T1)
    Trend <- matrix(1:T1, T1, 1)
    if ( model == "I" ) {
    Y0 = Y0
    Z1 = Z1
    Z2 = Z2
    } 
    if ( model == "II" ) {
    Y0 = Y0
    Z1 = cbind(lT,Z1)
    Z2 = Z2
    } 
    if ( model == "III" ) {
    Y0 = Y0
    Z1 = Z1
    if (!anyNA(s))  Z2 = cbind(Z2,St,NSt)   else  Z2 = cbind(Z2,lT) 
    } 
    if ( model == "IV" ) {
    Y0 = Y0
    Z1 = cbind(Trend,Z1)
    Z2 = cbind(Z2,lT)
    } 
    if ( model == "V" ) {
    Y0 = Y0
    Z1 = cbind(Z1)
    Z2 = cbind(Z2,lT,Trend)
    } 

        M00 <- crossprod(Y0)/T1
        M11 <- crossprod(Z1)/T1
        M22 <- crossprod(Z2)/T1
        M01 <- crossprod(Y0, Z1)/T1
        M02 <- crossprod(Y0, Z2)/T1
        M10 <- crossprod(Z1, Y0)/T1
        M20 <- crossprod(Z2, Y0)/T1
        M12 <- crossprod(Z1, Z2)/T1
        M21 <- crossprod(Z2, Z1)/T1
        M22inv <- solve(M22)
        R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
        R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))
    
    S00 <- crossprod(R0)/T1
    S01 <- crossprod(R0, R1)/T1
    S10 <- crossprod(R1, R0)/T1
    S11 <- crossprod(R1)/T1
    Ctemp <- chol(S11, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv <- solve(S00)
    valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))
 
    #valeigen <- eigen(S11inv %*% S10 %*% S00inv %*% S01 )
    lambda <- valeigen$values
    e      <- valeigen$vector
    
    V <- t(Cinv) %*% e     
    Vorg <- V
    V <- sapply(1:M1, function(j) V[, j]/V[1, j])
    W <- S01 %*% V %*% solve(t(V) %*% S11 %*% V)
    PI <- S01 %*% solve(S11)
    DELTA <- S00 - S01 %*% V %*% solve(t(V) %*% S11 %*% V) %*% t(V) %*% S10
    GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
    beta   <- as.matrix(V[,1:crk])
    beta0  <- as.matrix(V[,1:crk])
    if (crk>0) {
       
       betaS = nlm(f, as.vector(beta0), beta,Z1,St,NSt,Y0,Z2)$estimate
       dim(betaS) = dim(beta)
       for (i in 1:ncol(as.matrix(betaS))) betaS[,i] = betaS[,i]/(betaS[1,i])
       CI = Z1%*%betaS
       estimation<-lm(Y0~0+CI+Z2) 
       CI1 = CI*St
       CI2 = CI*NSt
       VECM1<-lm(Y0~0+CI1+CI2+Z2) 
    }
############ ADD in interation to find beta


############

    if (crk == 0 ) {
       CI = 0
       estimation<-lm(Y0~0+Z2) 
    } 
    E = -T1 * log(1 - lambda)
    E = E[1:N]
    ##### Estimation of VECM
  
    #####
    ##### transform to level model
    #####
    ##### 
    resultsvecm <- summary(estimation)
    #if ( model == "I" )   { Tab<-read.table("tab1.txt", header = FALSE) }
    #if ( model == "II" )  { Tab<-read.table("tab2.txt", header = FALSE) }
    #if ( model == "III" ) { Tab<-read.table("tab3.txt", header = FALSE) }
    #if ( model == "IV" )  { Tab<-read.table("tab4.txt", header = FALSE) }
    #if ( model == "V" )   { Tab<-read.table("tab5.txt", header = FALSE) }
	if ( model == "I" )   { Tab<- Tab1 }
	if ( model == "II" )  { Tab<- Tab2 }
	if ( model == "III" ) { Tab<- Tab3 }
	if ( model == "IV" )  { Tab<- Tab4 }
	if ( model == "V" )   { Tab<- Tab5 }



    
    b = c(1:12)
    for (i in 1:12)   { b[i] = Tab[2*(i-1)+1,2+M1-N1] }
    a = c(1:12)
    for (i in 1:12)   { a[i] = Tab[2*i,2+M1-N1] }
    
    if (type == "eigen") {
        critical_vals = b
        M = matrix(0, N, 1)
        j = 1
        rank = 0
        while (j <= N && E[j] > critical_vals[N + 1 - j]) {
            M[j, ] = M[j, ] + 1
            j = j + 1
            rank = rank + 1
        }
        if (ret == "test") {
            return(M)
        }
        erg <- cbind(E, critical_vals[N:1])
        colnames(erg) <- c("teststatistic", "critical_value")

        if ( N > 1 ) {
           rownames(erg) <- c("crk <= 0 |", paste("crk <= ", 1:(N - 1),  " |", sep = ""))
        }
        if ( N == 1 ) {
           rownames(erg) <- c("crk <= 0 |" )
        }

        coint_rank <- paste("Johansen-Test (with maximum-eigenvalue-teststatistic) indicates", 
            rank, "cointegrating equation(s) at the", 1 - q, 
            "level")
        #print(coint_rank)
        result <- new.env()
        result$erg = erg
        result$estimation = estimation
        result$VECM1      = VECM1
        result$lambda     = E
        result$z          = z
        result$Z2         = Z2
        result$beta       = beta

	  result$betaS      = betaS
        result$PI         = PI
        result$GAMMA      = GAMMA
        result$model      = model
        result$P          = P
        result$NN1        = NN1
        result$s          = s
        rr<-as.list(result)
        return(rr)
    }
    else {
        type = "trace"
        critical_vals = a
        stat = matrix(0, N, 1)
        for (i in 1:N) {
            sum = 0
            for (j in i:N) {
                sum = sum + E[j]
            }
            stat[i] = sum
        }
        M = matrix(0, N, 1)
        j = 1
        rank = 0
        while (stat[j] > critical_vals[N + 1 - j] && j <= N) {
            M[j, ] = M[j, ] + 1
            j = j + 1
            rank = rank + 1
        }
        if (ret == "test") {
            return(M)
        }
        erg <- cbind(stat, critical_vals[N:1])
        colnames(erg) <- c("teststatistic", "critical_value")

        if ( N > 1 ) {
           rownames(erg) <- c("crk <= 0 |", paste("crk <= ", 1:(N - 1),  " |", sep = ""))
        }
        if ( N == 1 ) {
           rownames(erg) <- c("crk <= 0 |")
        }
 
            coint_rank <- paste("Johansen-Test (with trace-teststatistic) indicates", 
            rank, "cointegrating equation(s) at the", 1 - q, 
            "level")
        #print(coint_rank)
        result <- new.env()
        result$erg = erg
        result$estimation = estimation
	  result$VECM1       = VECM1
        result$lambda     = E
        result$z          = z
        result$Z2         = Z2
        result$beta       = beta

	  result$betaS      = betaS
        result$PI         = PI
        result$GAMMA      = GAMMA
        result$model      = model
        result$P          = P
        result$NN1        = NN1
        result$s          = s
        rr<-as.list(result)
        return(rr)
    }
}




#' Regime specific impulse response functions of MRCIVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions with confidence bands, using Bo[,,,s] and Sigma[,,s] matrices of the estimated MRCIVAR, this function calculated the regime speicfic impulse response functions.
#'
#' @param res_e an object of MRCIVAR as output of MRVARestm 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol Cholezky decomposition, chol1 Cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @param conf a vector of the tail probabilities of the confidence interval. 
#' @return a list of impulse reaponse function in two regimes and the bootstrap parameters. 
#' @examples 
#'
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#' 
#' res_d = MRCIVARDatam(n=4,p=p,T=261,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
#' res_e = MRCIVARestm(res=res_d)
#' ## check the roots of the characteristic equation
#' STAT(res_d$Bo[,,,1])
#' STAT(res_d$Bo[,,,2])
#' plot(ts(res_d$Y))
#' ## Johasen test of the cointegration rank
#' res_e$res$tst[[1]]
#' ## Regression format of the VECM   
#' #summary(res_e$res$tst$VECM1)
#' ## within regime impulse response function
#' IRF   = irf_MRCIVARm(res_e,nstep=20,irf="gen1",runs=400,comb=NA,G=NA,conf=c(0.05,0.95))
#' ## For the impusle response function of regime 2: 
#' IRF_CB = IRF[[2]] 
#' ## For the impusle response function of regime 1: 
#' IRF_CB = IRF[[1]]    
#' x11()
#' par(mfrow=c(2,1))
#' plott(IRF_CB,1,1) 
#' plott(IRF_CB,2,1)
#' @export
irf_MRCIVARm = function(res_e,nstep=20,irf=c("gen","gen1"),runs=100,comb=NA,G=NA,conf=c(0.05,0.95)) {
	IRF1   = irf_B_sigma(B=res_e[[1]]$Bo[,,,1],sigma=res_e[[1]]$Sigma[,,1],nstep=nstep,comb=NA,irf=irf,G=NA)
	IRF2   = irf_B_sigma(B=res_e[[1]]$Bo[,,,2],sigma=res_e[[1]]$Sigma[,,2],nstep=nstep,comb=NA,irf=irf,G=NA)
	n 	 = res_e[[1]]$n
      p 	 = res_e[[1]]$p
      T 	 = dim(res_e[[1]]$Y)[1]
      S 	 = res_e[[1]]$S
      SV 	 = res_e[[1]]$SV
      SESVI  = res_e[[1]]$SESVI
      TH 	 = res_e[[1]]$TH
      Sigmao = res_e[[1]]$Sigmao
      Bo     = res_e[[1]]$Bo;
      Co     = res_e[[1]]$Co;
      type   = res_e[[1]]$type
      r      = res_e[[1]]$r
      Uo     = res_e[[1]]$Uo
      Y      = res_e[[1]]$Y
      res_d = MRCIVARDatam(n,p,T,S,SESVI,TH,Bo,Co,Sigmao,Uo = Uo,type,r)
      neq   = res_e[[1]]$n
      nvar  = res_e[[1]]$n
      
      response1 <- array(0,dim=c(neq,nvar,nstep,length(conf)+1))
      response2 <- array(0,dim=c(neq,nvar,nstep,length(conf)+1))


      response1[,,,1] <- IRF1
      response2[,,,1] <- IRF2

      responseR <- array(0,dim=c(neq,nvar,nstep,runs,S))
      BoColect = array(0,c(dim(Bo),runs))
      UoColect = array(0,c(dim(Uo),runs))
      CoColect = array(0,c(dim(Co),runs))
      YColect  = array(0,c(dim(Y),runs))
      Uo_run = array(0,c(T,n,S))
      
      for (i in 1:runs) {
            for (s in 1:S)    Uo_run[,,s] = rnormSIGMA(T,Sigmao[,,s])
            if (length(p)>1)	{  
               		res_run    = MRCIVARDatam(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Bo=Bo,Co=Co,Sigmao=Sigmao,Uo = Uo_run,type=type,r=r)

			  	res_e_run  = MRCIVARestm(res=res_run)
		}
		IRF1   = irf_B_sigma(B=res_e_run[[1]]$Bo[,,,1],sigma=res_e_run[[1]]$Sigma[,,1],nstep=nstep,comb=NA,irf=irf,G=NA)
		IRF2   = irf_B_sigma(B=res_e_run[[1]]$Bo[,,,2],sigma=res_e_run[[1]]$Sigma[,,2],nstep=nstep,comb=NA,irf=irf,G=NA)

            responseR[,,,i,1]  <- IRF1
            responseR[,,,i,2]  <- IRF2
            BoColect[,,,,i]    <- res_e_run[[1]]$Bo
  		CoColect[,,,i]     <- res_e_run[[1]]$Co

            UoColect[,,,i]     <-  res_run$Uo 
 		YColect[,,i]       <- res_run$Y                   
	}
      
      for (tt in 1:(nstep) )      {
	for (i in 1:neq)            {
		for (j in 1:nvar)     {
			response1[i,j,tt,-1] = quantile(responseR[i,j,tt,,1], conf)
			response2[i,j,tt,-1] = quantile(responseR[i,j,tt,,2], conf)

			
	} 
	}
	}
	return(list(response1,response2,BoColect,UoColect,YColect))
}



#' Calculation of the model selection criterion values for MRCIVAR models  
#' 
#' @param  res  : an object generated from CIGVARDatam or estimated from CIGVARestm.
#' @param  L_V  : a 2 components vector containing the maxima of lags in the two regimes. 
#' @param  TH_V  : a vector containing possible threshold values for  the selection. 

#' @return      : A matrix with different lag specifications and the model selection criteria.
#' @examples 
#' 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#'
#' res_d = MRCIVARData(n=4,p=p,T=2610,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
#' res_e = MRCIVARest(res=res_d)
#' 
#' TH_v = c(0,0.1)
#' L_v = c(6,6)
#' plot(ts(res_d$Y))
#' 
#' Sel = MRCIVAR_Select(res=res_e,L_V=L_v,TH_V=TH_v)
#' Sel[which.min(Sel[,5]),]
#' Sel[which.min(Sel[,4]),]
#' MRVAR_Select_Summary(Sel)
#' @export
MRCIVAR_Select = function(res=res_e,L_V=L_v,TH_V=TH_v) {
res_dd = res$res
p      = res_dd$p
n 	 = res_dd$n
T    	 = dim(res_dd$Y)[1]
SESVI  = res_dd$SESVI
type   = res_dd$type
TH     = res_dd$TH 
X      = res_dd$X
p[,1]  = max(L_V)
S      = res_dd$S
r      = res_dd$r
#res_dd = MRCIVARData(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Sigmao=NA,type=type,r=r)
        
Criteria = matrix(0,(L_V[1]-1)*(L_V[2]-1)*length(TH_v),9)
idx = 0

for (l_d in 2: L_V[1] )             {
   for (l_f in 2:L_V[2] )           {
      for (l_th in 1:length(TH_v))  {
      idx = idx + 1
      res_dd$p[1,1] = l_d
      res_dd$p[2,1] = l_f
      res_dd$TH     = TH_V[l_th]
	
	#res_dd = MRCIVARData(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Sigmao=NA,type=type,r=r)
	#res_dd$Y= res$res$Y
      
	# for MRCIVARestm
      #res_s         = MRCIVARestm(res=res_dd)
      #Criteria[idx,] = c (l_d,l_f,TH_V[l_th],res_s$LH_AIC,res_s$LH_BIC,res_s$LHH_AIC,res_s$LHH_BIC,res_s$ORAIC,res_s$ORBIC,res_s$TTT,res_s$LH_TN,res_s$LH_TN1,res_s$n2LHP)
      #colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LHH_AIC","LHH_BIC","ORAIC","ORBIC","TTT","LH_TN","LH_TN1","n2LHP")
      # for MECIVARest
      res_s         = MRCIVARest(res=res_dd)
	Criteria[idx,] = c (l_d,l_f,TH_V[l_th],res_s$LH_AIC,res_s$LH_BIC,res_s$LH_P,res_s$LH_N,res_s$ORAIC,res_s$ORBIC)
      colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LH_P","LH_N","ORAIC","ORBIC")
    }
}
}
return(Criteria)
}


#' @export
VECM2VARm = function (param, beta, p = c(1, 2, 2, 2, 2), s = NA) 
{
    m = dim(param)[2]
    VECB = t(param)
    if (anyNA(s)) {
        if (!p[3] == 0) {
            B = (1:(m * m * (p[2] + 1))) * 0
            dim(B) = c(m, m, (p[2] + 1))
            A = (1:(m * m * (p[3] + 1))) * 0
            dim(A) = c(m, m, (p[3] + 1))
            AB1 = VECB[, 1:p[1]] %*% t(beta)
		AB2 = VECB[, (1+p[1]):(p[1]+p[1])] %*% t(beta)
            B[, , 1] = AB[, 1:m]
            A[, , 1] = AB[, (m + 1):(2 * m)]
            for (i in 2:(p[2] + 1)) B[, , i] = VECB[, p[1]*2  + ((i - 2) * m + 1):((i - 2) * m + m)]
            for (i in 2:(p[3] + 1)) A[, , i] = VECB[, (p[1]*2 + p[2] * m) + ((i - 2) * m + 1):((i - 2) * m + m)]
            B = CIB3B(B)
            A = CIA2A(A)
        }
        else {
            B = (1:(m * m * (p[2] + 1))) * 0
            dim(B) = c(m, m, (p[2] + 1))
            AB = VECB[, 1:p[1]] %*% t(beta)
            B[, , 1] = AB[, 1:m]
            for (i in 2:(p[2] + 1)) B[, , i] = VECB[, p[1] + 
                ((i - 2) * m + 1):((i - 2) * m + m)]
            B = CIB3B(B)
            A = NA
        }
        if (dim(param)[1] > p[1] + (p[2] + p[3]) * m) 
            C = as.matrix(t(param)[, p[1] + (p[2] + p[3]) * m + 
                1])
        else C = NA
    }
    else {
        if (length(p) == 5) {   # still not cleaned
            P = max(p[-1])
            B = (1:(m * m * (P + 1) * 2)) * 0
            dim(B) = c(m, m, (P + 1), 2)
            A = B
            AB = VECB[, 1:p[1]] %*% t(beta)
            B[, , 1, 1] = AB[, 1:m]
            A[, , 1, 1] = AB[, (m + 1):(2 * m)]
            B[, , 1, 2] = AB[, 1:m]
            A[, , 1, 2] = AB[, (m + 1):(2 * m)]
            for (i in 2:(p[2] + 1)) B[, , i, 1] = VECB[, p[1] + 
                ((i - 2) * m + 1):((i - 2) * m + m)]
            for (i in 2:(p[3] + 1)) A[, , i, 1] = VECB[, (p[1] + 
                p[2] * m) + ((i - 2) * m + 1):((i - 2) * m + 
                m)]
            for (i in 2:(p[4] + 1)) B[, , i, 2] = VECB[, (p[1] + 
                (p[2] + p[3]) * m) + ((i - 2) * m + 1):((i - 
                2) * m + m)]
            for (i in 2:(p[5] + 1)) A[, , i, 2] = VECB[, (p[1] + 
                (p[2] + p[3] + p[4]) * m) + ((i - 2) * m + 1):((i - 
                2) * m + m)]
            B[, , 1:(p[2] + 1), 1] = CIB3B(B[, , 1:(p[2] + 1), 
                1])
            B[, , 1:(p[4] + 1), 2] = CIB3B(B[, , 1:(p[4] + 1), 
                2])
            A[, , 1:(p[3] + 1), 1] = CIA2A(A[, , 1:(p[3] + 1), 
                1])
            A[, , 1:(p[5] + 1), 2] = CIA2A(A[, , 1:(p[5] + 1), 
                2])
            if (dim(param)[1] > p[1] + (sum(p[-1])) * m) {
                C = (1:(m * 2)) * 0
                dim(C) = c(m, 1, 2)
                C[, 1, 1] = as.matrix(t(param)[, p[1]*2 + sum(p[-1]) * 
                  m + 1])
                C[, 1, 2] = as.matrix(t(param)[, p[1]*2 + sum(p[-1]) * 
                  m + 2])
            }
            else {
                C = NA
            }
        }
        if (length(p) == 3) {
            P = max(p[-1])
            B = (1:(m * m * (P + 1) * 2)) * 0
            dim(B) = c(m, m, (P + 1), 2)
   		AB1 = VECB[, 1:p[1]] %*% t(beta)
		AB2 = VECB[, (1+p[1]):(p[1]+p[1])] %*% t(beta)

            B[, , 1, 1] = AB1[, 1:m]
            B[, , 1, 2] = AB2[, 1:m]
            for (i in 2:(p[2] + 1)) B[, , i, 1] = VECB[, 2*p[1] + 
                ((i - 2) * m + 1):((i - 2) * m + m)]
            for (i in 2:(p[3] + 1)) B[, , i, 2] = VECB[, (2*p[1] + 
                p[2] * m) + ((i - 2) * m + 1):((i - 2) * m + 
                m)]
            B[, , 1:(p[2] + 1), 1] = CIB3B(B[, , 1:(p[2] + 1), 
                1])
            B[, , 1:(p[3] + 1), 2] = CIB3B(B[, , 1:(p[3] + 1), 
                2])
            A = NA
            if (dim(param)[1] > p[1]*2 + (sum(p[-1])) * m) {
                C = (1:(m * 2)) * 0
                dim(C) = c(m, 1, 2)
                C[, 1, 1] = as.matrix(t(param)[, p[1]*2 + sum(p[-1]) * 
                  m + 1])
                C[, 1, 2] = as.matrix(t(param)[, p[1]*2 + sum(p[-1]) * 
                  m + 2])
            }
            else {
                C = NA
            }
        }
    }
    return(list(B, A, C))
}





#' Calculation of the model selection values for MRVAR models  
#' 
#' @param  res  : an object generated from CIGVARData or estimated from CIGVARest.
#' @param  L_V  : a 2 components vector containing the maxima of the domestic lag and the foreign lag respectively. 
#' @param  TH_V  : a vector containing possible threshold values for the selection. 

#' @return      : A matrix with different lag specifications and the model selection criteria.
#' @examples 
#' 
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2) 
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(3,2)
#'
#' res_d = MRCIVARDatam(n=4,p=p,T=2610,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
#' res_e = MRCIVARestm(res=res_d)
#' 
#' TH_v = c(0,0.1)
#' L_v = c(6,6)
#' plot(ts(res_d$Y))
#' 
#' Selm = MRCIVAR_Selectm(res=res_e,L_V=L_v,TH_V=TH_v)
#' Selm[which.min(Selm[,5]),]
#' Selm[which.min(Selm[,4]),]
#' MRVAR_Select_Summary(Selm)
#'
#' @export
MRCIVAR_Selectm = function(res=res_e,L_V=L_v,TH_V=TH_v) {
res_dd = res$res
p      = res_dd$p
n      = res_dd$n
T      = dim(res_dd$Y)[1]
SESVI  = res_dd$SESVI
type   = res_dd$type
TH     = res_dd$TH 
X      = res_dd$X
p[,1]  = max(L_V)
S      = res_dd$S
r      = res_dd$r
#res_dd = MRCIVARData(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Sigmao=NA,type=type,r=r)
        
Criteria = matrix(0,(L_V[1]-1)*(L_V[2]-1)*length(TH_v),9)
idx = 0

for (l_d in 2: L_V[1] )             {
   for (l_f in 2:L_V[2] )           {
      for (l_th in 1:length(TH_v))  {
      idx = idx + 1
      res_dd$p[1,1] = l_d
      res_dd$p[2,1] = l_f
      res_dd$TH     = TH_V[l_th]
	
	#res_dd = MRCIVARData(n=n,p=p,T=T,S=S,SESVI=SESVI,TH=TH,Sigmao=NA,type=type,r=r)
	#res_dd$Y= res$res$Y
      
	# for MRCIVARestm
      res_s         = MRCIVARestm(res=res_dd)
      Criteria[idx,] = c (l_d,l_f,TH_V[l_th],res_s$LH_AIC,res_s$LH_BIC,res_s$LHH_AIC,res_s$LHH_BIC,res_s$ORAIC,res_s$ORBIC)
      colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LHH_AIC","LHH_BIC","ORAIC","ORBIC")
      # for MECIVARest
      #res_s         = MRCIVARest(res=res_dd)
	#Criteria[idx,] = c (l_d,l_f,TH_V[l_th],res_s$LH_AIC,res_s$LH_BIC,res_s$LH_P,res_s$LH_N,res_s$ORAIC,res_s$ORBIC)
      #colnames(Criteria) = c("Lag_regime1","Lag_regime2" ,"threshold","AIC","BIC","LH_P","LH_N","ORAIC","ORBIC")
    }
}
}
return(Criteria)
}



#' @export
f <- function(x,beta,Z1,St,NSt,Y0,Z2) {
      dim(x) = dim(beta)
	CI = Z1%*%x
      CI1 = CI*St
      CI2 = CI*NSt
      residuals <-lm(Y0~0+CI1+CI2+Z2)$residuals 
      SSR = sum(diag(t(residuals)%*%residuals))
      return(SSR)
} 




