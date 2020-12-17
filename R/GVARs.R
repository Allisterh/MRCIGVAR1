#### GVAR 
#### n number of countries
#### m number of variables
#### p number of lags
#### w weiting matrix
#### 


#' Data generating process of GVAR(n,m,p) 
#'
#' This function will generate data from a staionary global autoregressive process GVAR(p) and return an GVAR object containing data and parameters used in the GVAR(p) process.
#'
#' @param n     : number of countries/units
#' @param m     : number of variables in each unit
#' @param p     : lag length, the lag length is assumed to be identical across all countries for both the domestic and foreign variables  
#' @param T     : number of observations
#'                (m,n,p,T) are parameters which must be provided. 
#' @param W     : n x n weighting matrix. w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' r_npo 	    : n x m x p array collecting the roots of the characteristic functions of the country VAR
#' Ao		    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
#' Bo    	    : m x m x p x n array collecting the n country VAR(p) coefficients.  Bo are coefficents of stationary domestic VAR(p). 
#' Co           : m x (k+1) x n array collecting the coeeficients of the deterministic components of the n countries.
#' Uo		    : an T x mn matrix of the temporally independent innovation processes
#' Sigmao	    : mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
#'  		      (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly.
#'
#' type  	    : deterministic component "const" and "none" are two options 
#' mu    	    : if type = "const" mu has the same dimension as Co. is an muV is nm vector of the means of the time series in the system
#'
#' @param Sigma : The covariance matrix of the n dynamically independent processes
#' @param A     : an n x n full rank matrix of transformation to generate correlated VAR(p) from the n independent AR(p)  
#' @param B     : (n,n,p) array of the AR(p) process. If B is not given, it will be calculated out of r_np and A.
#' @param Co    : (m,k+1,n) vector of intercept of the VAR(p) process. for type="none" Co = O*(1:n), for const Co is an n vector, exog0 Co is a (n,k) matrix for exog1 Co is an (n,1+k) matrix  Depending on type it will be zeros for none,  
#' @param type  : deterministic component "none", "const" "exog0" and "exog1" are four options
#' @param U     : residuals, if it is not NA it will be used as input to generate the VAR(p)
#' @param X     : (T x k x n) matrix of exogeneous variables. For each country there is a (T x k) exogenous variable. 
#' @param Yo    : (p x n) matrix of initial values of the process
#' @param mu    : n vector of the expected mean of the VAR(p) process
#' @return        A list containing the generated data, the parameters and the input exogeous variables. res = list("Y","Uo","G","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
#'
#' @param Y     : T x nm simulated data via of the GVAR(m,n,p) 
#' @param Uo    : T x mn array of the simulated innovations of the GVAR(m,n,p) 
#' @param G     : mn x mn x p  array of the GVAR(m,n,p) coefficients. G is contructed from Bo, Ao and W:  
#' @param         off-diagonal          G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L]  = Ao[,,L,i]*W[i,j]
#' @param         diagonal              G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]  = Bo[,,L,i] 
#' @param Sigmao: mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
#' @param r_npo : m x p x n matrix collecting the roots of the characteristic functions in L of the n dynamically independent domestic VAR(p)s.
#' @param Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients 
#' @param Bo    : m x m x p X n array collecting the n country VAR(p) coefficients. Bo are coefficents of stationary domestic VAR(p). 
#' @param W     : n x n weiting matrix, w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' @param check : maximum of the data to check stationarity
#' @examples 
#' res = GVARData0(m=2,n=4,p=3,T=200)
#' res = GVARData0(m=2,n=4,p=3,T=200,type="const")
#' plot(res$Y[,1],type="l")
#' lines(res$Y[,2],col="red")
#' lines(res$Y[,3],col="blue")
#' lines(res$Y[,4],col="green")
#'
#' mu=c(50,50,50,50,50,50,50,50)
#' dim(mu) = c(2,4)
#' Co = c(50,50,50,50,50,50,50,50)
#' dim(Co) = c(2,4)
#' res = GVARData0(m=2,n=4,p=3,T=200,type="const",mu=mu)
#' res = GVARData0(m=2,n=4,p=3,T=200,type="const",Co=Co)
#' plot(res$Y[,1],type="l")
#' lines(res$Y[,2],col="red")
#' lines(res$Y[,3],col="blue")
#' lines(res$Y[,4],col="green")
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X=cbind(X1,X2)
#' Co = c(1:(2*3*4))
#' dim(Co) = c(2,3,4)
#' res = GVARData0(m=3,n=4,p=3,T=200,type="exog",Co=Co,X=X)
#' @export
GVARData0 = function(m,n,p,T,W=NA,r_npo=NA,Ao=NA,Bo=NA,Go=NA,Co=NA,Uo=NA,Sigmao=NA,type=NA,X=NA,mu=NA) {
### m     : number of variables in a country
### n     : number of countries
### p     : lag length
### T     : number of observations
###         (m,n,p,T) are parameters which must be provided.
### W     : n x n weiting matrix, w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
### r_npo : n x m x p array collecting the roots of the characteristic functions of the country VAR
### Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
### Bo    : m x m x p x n array collecting the n country VAR(p) coefficients.  Bo are coefficents of stationary domestic VAR(p). 
### Uo    : an T x mn matrix of the temporally independent innovation processes
### Sigmao: mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
###         (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly.
###
### type  : deterministic component "const" and "none" are two options 
### mu    : if type = "const" mu is an nm vector of the means of the time series in the system
### output:
###  c("Y","Uo","G","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
###
### Y     : T x nm simulated data via of the GVAR(m,n,p) 
### Uo    : T x mn array of the simulated innovations of the GVAR(m,n,p) 
### G     : mn x mn x p  array of the GVAR(m,n,p) coefficients. B is contructed from Bo, Ao and W:  
###         off-diagonal          B[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L]  = Ao[,,L,i]*W[i,j]
###         diagonal              B[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]  = Bo[,,L,i] 
### Sigmao: mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
### r_npo : m x p x n matrix collecting the roots of the characteristic functions in L of the n dynamically independent domestic VAR(p)s.
### Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients 
### Bo    : m x m x p x n array collecting the n country VAR(p) coefficients. Bo are coefficents of stationary domestic VAR(p). 
### W     : n x n weiting matrix, w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
### m     : number of variables in a country
### n     : number of countries
### check : maximum of the data to check stationarity
###
### 
###
### Remarks: The n domestic VAR(p)s are A transformed from n dynamically independent ar processes, such that the stationary is guaranteed. 
###          But the coefficients of domestic VAR(p)s are restricted the A tranformation. 
###                        
  if (missing(Bo)) 	{Bo=NA}
  if (missing(Sigmao)) 	{Sigmao=NA}
  if (missing(Uo)) 	{Uo=NA}
  if (missing(W)) 	{W=NA}
  if (missing(Ao)) 	{Ao=NA}
  if (missing(Go)) 	{Go=NA}
  if (missing(type)) 	{type=NA}
  if (missing(Co)) 	{Co=NA}
  if (missing(mu)) 	{mu=NA}


  if (anyNA(Bo)) {
  	Bo = (1:(m*m*p*n))*0
      dim(Bo) = c(m,m,p,n)
      r_npo  = c(1:(m*p*n))*0
      dim(r_npo) = c(m,p,n)

      for (i in 1:n) {
      	VARD = VARData(m,p,T) 
      	Bo[,,,i]   = VARD$B
            r_npo[,,i] = VARD$r_np
      }      
  }

  if (anyNA(Sigmao)) {
  	Sigmao = matrix(0,m*n,m*n)
      VARD = VARData(m*n,p,T) 
      Sigmao = VARD$Sigma
  }

  if (anyNA(Uo)) {
      Uo= rnormSIGMA(T,Sigmao)
  }
  
  
  if (anyNA(W)) {
  	W = matrix(((1:(n*n))/(1:(n*n))*1/(n-1)),n,n)
      for (i in 1:n) W[i,i]= 0
  }

  if (anyNA(Ao)) {
  	Ao = (1:(m*m*p*n))*0
      dim(Ao) = c(m,m,p,n)
      for (i in 1:n) {
		for (L in 1:p) Ao[,,L,i] = matrix((runif(m*m)-0.5)/10,m,m)
	}
  }


  if (anyNA(type)) {
  		type = "none"
  }
  
  if (type=="none") {
  		Co = matrix(0,m,n); dim(Co) = c(m,1,n) 
            mu = matrix(0,m,n); dim(mu) = c(m,1,n)
  }
  
  G = (1:(n*m*n*m*p))*0
  dim(G) = c(n*m,n*m,p)
  for (i in 1:n) {
      for (j in 1:n) { 
         for (L in 1:p)    G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L] = Ao[,,L,i]*W[i,j]
      }
  }
  for (i in 1:n) {
	for (L in 1:p)       G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L] = Bo[,,L,i]  
  }
  if (!anyNA(Go)) G = Go

  Ct = Uo*0
  if (type=="const") { 
  	if (anyNA(mu))   {  mu = matrix(rnorm(n*m),m,n); dim(mu) = c(m,1,n)}
	if (anyNA(Co))   {
            Co = mu
            muV = as.vector(mu)
            CoV = muV
 		for (L in 1:p) CoV = CoV - G[,,L]%*%muV
      }   else   {
            H = diag(n*m)
		for (L in 1:p) H = H - G[,,L]
            CoV = as.vector(Co)
            muV = solve(H)%*%CoV
            mu  = muV
            dim(mu) = c(m,1,n)
      }
      Ct = matrix(1,T,1)%*%t(CoV)   
  }

  if (type=="exog0" ) 	{
            Co[,1,] = (1:m)*0
            k = dim(X)[2]
            DMCo = dim(Co[,-1,])
            if (!(DMCo[2]==dim(X)[2])|!(DMCo[1]==m)|!(DMCo[3]==n)) {print("dimension problem"); return("dimension")}
            CoV = matrix(0,dim(X)[2],m*n);
            for (i in 1:dim(X)[2])   CoV[i,] = as.vector(Co[,1+i,])
            Ct = X%*%CoV;         #
		mu = NA 
  }   else 	{ 
 
  	if (type=="exog1" ) 	{
            k = dim(X)[2]
            DMCo = dim(Co[,-1,])
            if (!(DMCo[2]==dim(X)[2])|!(DMCo[1]==m)|!(DMCo[3]==n)) {print("dimension problem"); return("dimension")}
            CoV = matrix(0,dim(X)[2]+1,m*n);
            for (i in 1:(1+dim(X)[2]))   CoV[i,] = as.vector(Co[,i,])
 
            Ct = cbind(matrix(1,dim(X)[1],1),X)%*%CoV;       #
		mu = NA 
  	}   else 			{ 
		X  = NA  
  	}
  }
  Y = Uo + Ct  
  for ( t in (p+1):T )  {
      for (L in 1:p)    Y[t,] = Y[t,] + Y[t-L,]%*%t(G[,,L])
  }

  check = max(abs(Y))
  result=list(Y,X,Uo,G,C,Sigmao,r_npo,Ao,Bo,Co,W,m,n,p,mu,check,type)
  names(result) = c("Y","X","Uo","G","C","Sigmao","r_npo","Ao","Bo","Co","W","m","n","p","mu","check","type")
  return(result)
}

#' Estimation of GVAR(m,n,p): GVARest(res) 
#'
#' This function estimate the unknown parameters of a specified GVAR(m,n,p) model based on provided data.
#'
#' @param  res  :a list containing the components as the output of GVARData0 including at least: m,n,p,type,Y and optionally X. 
#' @return res  :a list as the input, but filled with estimated parameter values, AIC and LH
#' @examples 
#' res = GVARData0(m=2,n=4,p=3,T=200)
#' GVARest0(res)
#' res = GVARData0(m=2,n=4,p=3,T=200,type="const")
#' GVARest0(res)
#' plot(res$Y[,1],type="l")
#' lines(res$Y[,2],col="red")
#' lines(res$Y[,3],col="blue")
#' lines(res$Y[,4],col="green")
#'
#' mu=c(50,50,50,50,50,50,50,50);dim(mu) = c(2,4)
#' Co = c(50,50,50,50,50,50,50,50);dim(Co) = c(2,4)
#' res = GVARData0(m=2,n=4,p=3,T=200,type="const",mu=mu)
#' GVARest0(res)
#'
#' res = GVARData0(m=2,n=4,p=3,T=200,type="const",Co=Co)
#' GVARest0(res)
#' plot(res$Y[,1],type="l")
#' lines(res$Y[,2],col="red")
#' lines(res$Y[,3],col="blue")
#' lines(res$Y[,4],col="green")
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X=cbind(X1,X2)
#' Co = c(1:(2*3*4))
#' dim(Co) = c(2,3,4)
#' res = GVARData0(m=3,n=4,p=3,T=200,type="exog",Co=Co,X=X)
#' GVARest0(res)
#' @export
GVARest0 = function(res) {
m	= res$m
n	= res$n
p	= res$p
Y	= res$Y
X     = res$X
W 	= res$W
type  = res$type
Bo	= res$Bo
Ao	= res$Ao
Co    = res$Co
resid = Y*0
VAR_domestic = list()


T	= dim(Y)[1]
FY	= Y%*%t(W%x%diag(m))

Bo    = Bo = (1:(m*m*p*n))*0
dim(Bo) = c(m,m,p,n)
Ao 	= Bo*0

if (type == "none")  C   = matrix(0,m*n,1)*0 
if (type == "const") C   = matrix(0,m*n,1)*0 
if (type == "exog0") C   = matrix(0,m*n,dim(X)[2]+1)*0 
if (type == "exog1") C   = matrix(0,m*n,dim(X)[2]+1)*0 

res_VAR = VARData(n=m,p=p,T=T) ### note the different use of (m n) in GVAR an n in VAR 

for (i in 1:n) {
	Z         	=  matrix(0,T,p*m)
      FYp	    	=  embed(FY[,(m*(i-1)+1):(i*m)],(p+1))
	Z[(p+1):T,] =  FYp[,(m+1):(p*m+m)]
      if (type=="none")   Z = Z
	if (type=="const")  Z = Z
      if (type=="exog0" ) Z = cbind(Z,X) 	              	
      if (type=="exog1" ) Z = cbind(Z,X) 	              	
	### put data into a VARData object in order to replace the data and order parameters
	res_VAR$Y 		=  Y[,(m*(i-1)+1):(i*m)]
	res_VAR$X		=  Z
	if ((type=="none" )|(type=="exog0")) {res_VAR$type 	=  "exog0"  }
	if ((type=="const")|(type=="exog1")) {res_VAR$type 	=  "exog1"  }

	RR			=  VARest(res_VAR)
      VAR_domestic[[i]] =  RR$varp
	Bo[,,,i]		=  RR$B
      if (type=="none")  Ai	=  RR$Co[,-1]
      if (type=="const") Ai	=  RR$Co[,-1]
      if (type=="exog0") Ai	=  RR$Co[,c(2:(m*p+1))]
      if (type=="exog1") Ai	=  RR$Co[,c(2:(m*p+1))]

	     #Ai		=  t(Ai)
       dim(Ai)		=  c(m,m,p)
	Ao[,,,i]		=  Ai
	if (type=="none")  { C[(m*(i-1)+1):(i*m),] = RR$Co[,1]; Co[,,i]=RR$Co[,1] }
	if (type=="const") { C[(m*(i-1)+1):(i*m),] = RR$Co[,1]; Co[,,i]=RR$Co[,1] }

	if (type=="exog0" ) { C[(m*(i-1)+1):(i*m),] = RR$Co[,c(1,(m*p+2):(m*p+1+dim(X)[2]))];
				    Co[,,i]=RR$Co[,c(1,(m*p+2):(m*p+1+dim(X)[2]))]
      }
	if (type=="exog1" ) { C[(m*(i-1)+1):(i*m),] = RR$Co[,c(1,(m*p+2):(m*p+1+dim(X)[2]))];
				    Co[,,i]=RR$Co[,c(1,(m*p+2):(m*p+1+dim(X)[2]))]
      }
	resid[,(m*(i-1)+1):(i*m)] = RR$resid           
}
	Sigmao            = t(resid)%*%resid/(T-2*m*p)

	### put the estimates into the GVARData0 object 
	### c("Y","Uo","B","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
    	G = BoAoW2G(Bo,Ao,W,m,n,p)
      Gs = diag(n*m)
	for (L in 1:p) { Gs = Gs - G[,,L] }

	res$G	= G
	res$C = C
	res$Sigma = Sigmao
	res$r_npo = NA
	res$Ao    = Ao
	res$Bo    = Bo
      res$Co    = Co
 	res$mu    = solve(Gs)%*%C
      res$VAR_domestic = VAR_domestic
	return(res)
}



#' Data generating process of GVAR(n,m,p) 
#'
#' This function will generate data from a staionary global vector autoregressive process GVAR and return a GVAR object containing data and parameters used in the GVAR(p) process.
#'
#' @param m     : number of variables
#' @param n     : number of countries/units
#' @param p     : n x 3 matrix, each row contains the lag of the domestic variables, the lag of the foreign variables and the number of exogeneous variables. 
#' @param T     : number of observations
#' 
#' @param (m,n,p,T) are parameters which must be provided. 
#' @param W     : n x n weighting matrix. w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' @param r_npo : (m, p, n) array collecting the roots of the characteristic functions of Lag for each of the m domestic variable across n countries.
#' @param Ao    : (m, m, p, n) array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
#' @param Bo    : (m, m, p, n) array collecting the n country VAR(p) coefficients.  Bo are coefficents of stationary domestic VAR(p). 
#' @param Co    : (m, k+1, n) array collecting the coeeficients of the deterministic components of the n countries.
#' @param Uo    : an T x mn matrix of the temporally independent innovation processes
#' @param Sigmao: mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
#' @param (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly.
#' @param type  : deterministic component "const" and "none" are two options 
#' @param X     : (T x k) matrix of exogeneous variables.
#' @param mu    : if type = "const" mu has the same dimension as Co. is an muV is nm vector of the means of the time series in the system
#' @return      A GVAR object which is a list("Y","X","Uo","G","C","Sigmao","r_npo","Ao","Bo","Co","W","m","n","p","mu","check","type") containing the generated data, the parameters and the input exogeous variables. res = list("Y","Uo","G","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu") 
#' @param Y     : T x nm generated data via of the GVAR(m,n,p) 
#' @param X     : (T x k) matrix of exogeneous variables.
#' @param Uo    : T x mn array of the simulated innovations of the GVAR(m,n,p) 
#' @param G     : mn x mn x p  array of the GVAR(m,n,p) coefficients. G is contructed from Bo, Ao and W:  
#' @param         off-diagonal          G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L]  = Ao[,,L,i]*W[i,j]
#' @param         diagonal              G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]  = Bo[,,L,i] 
#' @param C     : the nm x 1 deterministic component of GVAR
#' @param Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients 
#' @param Bo    : m x m x p X n array collecting the n country VAR(p) coefficients. Bo are coefficents of stationary domestic VAR(p). 
#' @param check : maximum of the data for checking the stationarity
#' @examples 
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 2; p[,2]=1; 
#' res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
#' max(res_d$Y)
#' dim(res_d$Y)
#' res_e = GVARest(res = res_d)
#' summary_GVAR(res_e)  # the summary needs to be detailes commented
#'
#'
#'
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X3 = matrix(rnorm(200),200,1)
#' X4 = matrix(rnorm(200),200,1)
#' X  = cbind(X1,X2,X3,X4)
#' dim(X) = c(200,1,4)
#'
#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1;   p[,3]=1; p[2,2]=2;
#' p
#'
#' res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
#' summary_GVAR(res_e)  # the summary needs more detailed comments
#'
#' @export
GVARData = function(m,n,p,T,W=NA,r_npo=NA,Ao=NA,Bo=NA,Co=NA,Uo=NA,Sigmao=NA,type=NA,X=NA,mu=NA,d=d) {
### m     : number of variables in a country
### n     : number of countries
### p     : (n x 3) matrix, each raw contains the lag length of the domestic variables, the foreign variables and the number of the exogeneous variables
### T     : number of observations
###         (m,n,p,T) are parameters which must be provided.
### W     : n x n weiting matrix, w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
### r_npo : n x m x p array collecting the roots of the characteristic functions of the country VAR
### Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
### Bo    : m x m x p x n array collecting the n country VAR(p) coefficients.  Bo are coefficents of stationary domestic VAR(p). 
### Uo    : an T x mn matrix of the temporally independent innovation processes
### Sigmao: mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
###         (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly.
###
### type  : deterministic component "const" and "none" are two options 
### mu    : if type = "const" mu is an nm vector of the means of the time series in the system
### X     : T x k x n  array of exogenous varibales for each country. the number of exogenous variables can be different for each country, indicating by p[,3].
### output:
###  c("Y","Uo","G","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
###
### Y     : T x nm simulated data via of the GVAR(m,n,p) 
### Uo    : T x mn array of the simulated innovations of the GVAR(m,n,p) 
### G     : mn x mn x p  array of the GVAR(m,n,p) coefficients. B is contructed from Bo, Ao and W:  
###         off-diagonal          B[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L]  = Ao[,,L,i]*W[i,j]
###         diagonal              B[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]  = Bo[,,L,i] 
### Sigmao: mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
### r_npo : m x p x n matrix collecting the roots of the characteristic functions in L of the n dynamically independent domestic VAR(p)s.
### Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients 
### Bo    : m x m x p x n array collecting the n country VAR(p) coefficients. Bo are coefficents of stationary domestic VAR(p). 
### W     : n x n weiting matrix, w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
### m     : number of variables in a country
### n     : number of countries
### check : maximum of the data to check stationarity
###
### 
###
### Remarks: The n domestic VAR(p)s are A transformed from n dynamically independent ar processes, such that the stationary is guaranteed. 
###          But the coefficients of domestic VAR(p)s are restricted the A tranformation. 
###                        
  if (missing(Bo)) 	{Bo=NA}
  if (missing(Sigmao)) 	{Sigmao=NA}
  if (missing(Uo)) 	{Uo=NA}
  if (missing(W)) 	{W=NA}
  if (missing(Ao)) 	{Ao=NA}
  if (missing(type)) 	{type=NA}
  if (missing(Co)) 	{Co=NA}
  if (missing(mu)) 	{mu=NA}
  if (missing(d)) 	{d=NA}
  if (anyNA(d)) d = 1
  Pmax = max(p[,1:2])
  P    = max(p,d) 

  
  if (!anyNA(X)) k    = dim(X)[2]
  if (anyNA(Bo)) {
  	Bo = (1:(m*m*Pmax*n))*0
      dim(Bo) = c(m,m,Pmax,n)
      r_npo  = c(1:(m*Pmax*n))*0
      dim(r_npo) = c(m,Pmax,n)

      for (i in 1:n) {
      	VARD = VARData(m,p[i,1],T) 
      	Bo[,,1:p[i,1],i]   = VARD$B;     
            r_npo[,1:p[i,1],i] = VARD$r_np
      }      
  }

  if (anyNA(Sigmao)) {
  	Sigmao = matrix(0,m*n,m*n)
      VARD = VARData(m*n,p[1,1],T) 
      Sigmao = VARD$Sigma
  }

  if (anyNA(Uo)) {
      Uo= rnormSIGMA(T,Sigmao)
  }
  
  
  if (anyNA(W)) {
  	W = matrix(((1:(n*n))/(1:(n*n))*1/(n-1)),n,n)
      for (i in 1:n) W[i,i]= 0
  }

  if (anyNA(Ao)) {
  	Ao = (1:(m*m*Pmax*n))*0
      dim(Ao) = c(m,m,Pmax,n)
      for (i in 1:n) {
		for (L in 1:p[i,2]) Ao[,,L,i] = matrix((runif(m*m)-0.5)/10,m,m)
	}
  }


  if (anyNA(type)) {
  		type = "none"
  }

  if (type=="none") {
  		Co = matrix(0,m,n); dim(Co) = c(m,1,n) 
            mu = matrix(0,m,n); dim(mu) = c(m,1,n)
  }
  
  G = (1:(n*m*n*m*Pmax))*0
  dim(G) = c(n*m,n*m,Pmax)
  for (i in 1:n) {
      for (j in 1:n) { 
         for (L in 1:Pmax)    G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L] = Ao[,,L,i]*W[i,j]
      }
  }
  for (i in 1:n) {
	for (L in 1:Pmax)       G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L] = Bo[,,L,i]  
  }
  

  Ct = Uo*0
  if (type=="const") { 
  	if (anyNA(mu))   {  mu = matrix(rnorm(n*m),m,n); dim(mu) = c(m,1,n)}
	if (anyNA(Co))   {
            Co = mu
            muV = as.vector(mu)
            CoV = muV
 		for (L in 1:Pmax) CoV = CoV - G[,,L]%*%muV
      }   else   {
            H = diag(n*m)
		for (L in 1:Pmax) H = H - G[,,L]
            CoV = as.vector(Co)
            muV = solve(H)%*%CoV
            mu  = muV
            dim(mu) = c(m,1,n)
      }
      Ct = matrix(1,T,1)%*%t(CoV)   
  }

  if (type=="exog0" ) 	{
            if (anyNA(Co)) { 
			Co = matrix(rnorm(m*n*(k+1)),m*(k+1),n); 
			dim(Co) = c(m,k+1,n);      
			Co[,1,] = (1:m)*0
			for (i in 1:n) if (p[i,3]<k) Co[,(p[i,3]+2):(k+1),i] = Co[,(p[i,3]+2):(k+1),i]*0  
             }
            
            
            DMCo = dim(Co[,-1,])
            if (length(DMCo)<3) DMCo = c(dim(Co)[1],1,dim(Co)[3])
            if (!(DMCo[2]==dim(X)[2])|!(DMCo[1]==m)|!(DMCo[3]==n)) {print("dimension problem"); return("dimension")}
            CoV = matrix(0,dim(X)[2],m*n);
            for (i in 1:dim(X)[2])   CoV[i,] = as.vector(Co[,1+i,])
            Ct = matrix(0,dim(X)[1],m*n)
            for (i in 1:n) Ct[,((i-1)*m+1):((i-1)*m+m)] = as.matrix(X[,,i])%*%CoV[,((i-1)*m+1):((i-1)*m+m)];       #
		mu = NA 
  }   else 	{ 
 
  	if (type=="exog1" ) 	{
            if (anyNA(Co)) { 
			Co = matrix(rnorm(m*n*(k+1)),m*(k+1),n); 
			dim(Co) = c(m,k+1,n);
			for (i in 1:n) if (p[i,3]<k) Co[,(p[i,3]+2):(k+1),i] = Co[,(p[i,3]+2):(k+1),i]*0  
		}

            DMCo = dim(Co[,-1,])
            #if (!(DMCo[2]==dim(X)[2])|!(DMCo[1]==m)|!(DMCo[3]==n)) {print("dimension problem"); return("dimension")}
            CoV = matrix(0,dim(X)[2]+1,m*n);
            for (i in 1:(1+dim(X)[2]))   CoV[i,] = as.vector(Co[,i,])
            Ct = matrix(0,dim(X)[1],m*n)
            for (i in 1:n) Ct[,((i-1)*m+1):((i-1)*m+m)] = cbind(matrix(1,dim(X)[1],1),X[,,i])%*%CoV[,((i-1)*m+1):((i-1)*m+m)];       #
		mu = NA 
  	}   else 			{ 
		X  = NA  
  	}
  }
  Y = Uo + Ct  
  for ( t in (P+1):T )  {
      for (L in 1:Pmax)    Y[t,] = Y[t,] + Y[t-L,]%*%t(G[,,L])
  }

  check = max(abs(Y))
  result=list(Y,X,Uo,G,C,Sigmao,r_npo,Ao,Bo,Co,W,m,n,p,mu,check,type)
  names(result) = c("Y","X","Uo","G","C","Sigmao","r_npo","Ao","Bo","Co","W","m","n","p","mu","check","type")
  return(result)
}

#' Estimation of GVAR(m,n,p) 
#'
#' This function estimates the unknown parameters of a specified GVAR(m,n,p) model based on provided data.
#'
#' @param  res  : an GVAR object that is an output of GVARData including at least: m,n,p,type,Y and optionally X. 
#' @return res  : an GVAR object with estimated parameter values, AIC, BIC, AIC_g, BIC_g and LH, where AIC and BIC are the sum of the country equations AIC and BIC and AIG_g amd BIC_g are the GVAR information criteria respectively.
#' @examples 
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 2; p[,2]=1; 
#' res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
#' max(res_d$Y)
#' dim(res_d$Y)
#' res_e = GVARest(res = res_d)
#' summary_GVAR(res_e)  # the summary needs more detailed comments.
#'
#'
#' @export
GVARest = function(res) {
m	= res$m
n	= res$n
p	= res$p  # p is an n x 3 matrix of lag length of domestic, foreign and exogenoues varibles
Y	= res$Y
X     = res$X  # T x k x n
W 	= res$W
type  = res$type
Bo	= res$Bo
Ao	= res$Ao
Co    = res$Co


Pmax  = max(p[,1:2])
k     = max(p[,3])

AIC   = c(1:n)*0
BIC   = c(1:n)*0
LH    = c(1:n)*0

resid = Y*0
VAR_domestic = list()


T	= dim(Y)[1]
FY	= Y%*%t(W%x%diag(m))

Bo    = Bo = (1:(m*m*Pmax*n))*0
dim(Bo) = c(m,m,Pmax,n)
Ao 	= Bo*0


if (type == "none")  C   = matrix(0,m*n,1)*0 
if (type == "const") C   = matrix(0,m*n,1)*0 
if (type == "exog0") C   = matrix(0,m*n,dim(X)[2]+1)*0 
if (type == "exog1") C   = matrix(0,m*n,dim(X)[2]+1)*0 

#res_VAR = VARData(n=m,p=p,T=T) ### note the different use of (m n) in GVAR an n in VAR 

for (i in 1:n) {
      res_VAR = VARData(n=m,p=p[i,1],T=T) ### note the different use of (m n) in GVAR an n in VAR 
	Z         	=  matrix(0,T,p[i,2]*m)
      FYp	    	=  embed(FY[,(m*(i-1)+1):(i*m)],(Pmax+1))
	Z[(Pmax+1):T,] =  FYp[,(m+1):(p[i,2]*m+m)]
      if (type=="none")   Z = Z
	if (type=="const")  Z = Z
      if (type=="exog0" ) Z = cbind(Z,X[,1:p[i,3],i]) 	              	
      if (type=="exog1" ) Z = cbind(Z,X[,1:p[i,3],i])
      ki = dim(Z)[2] 	              	
	### put data into a VARData object in order to replace the data and order parameters
	res_VAR$Y 		=  Y[,(m*(i-1)+1):(i*m)]
	res_VAR$X		=  Z
	if ((type=="none" )|(type=="exog0")) {res_VAR$type 	=  "exog0"  }
	if ((type=="const")|(type=="exog1")) {res_VAR$type 	=  "exog1"  }

	RR			=  VARest(res_VAR)
      AIC[i]            =  RR$AIC
      BIC[i]            =  RR$BIC
      LH[i]             =  RR$LH

      VAR_domestic[[i]] =  RR$varp
	Bo[,,1:p[i,1],i]		=  RR$B
      if (type=="none")  Ai	=  RR$Co[,-1]
      if (type=="const") Ai	=  RR$Co[,-1]
      if (type=="exog0") Ai	=  RR$Co[,c(2:(m*p[i,2]+1))]
      if (type=="exog1") Ai	=  RR$Co[,c(2:(m*p[i,2]+1))]

	     #Ai		=  t(Ai)
       dim(Ai)		=  c(m,m,p[i,2])
	Ao[,,1:p[i,2],i]		=  Ai
	if (type=="none")  { C[(m*(i-1)+1):(i*m),] = RR$Co[,1]; Co[,,i]=RR$Co[,1] }
	if (type=="const") { C[(m*(i-1)+1):(i*m),] = RR$Co[,1]; Co[,,i]=RR$Co[,1] }

	if (type=="exog0" ) { C[(m*(i-1)+1):(i*m),1:(1+p[i,3])] = RR$Co[,c(1,(m*p[i,2]+2):(m*p[i,2]+1+p[i,3]))];
				    Co[,1:(p[i,3]+1),i]=RR$Co[,c(1,(m*p[i,2]+2):(m*p[i,2]+1+p[i,3]))]
      }
	if (type=="exog1" ) { C[(m*(i-1)+1):(i*m),] = RR$Co[,c(1,(m*p[i,2]+2):(m*p[i,2]+1+p[i,3]))];
				    Co[,              ,i] =  RR$Co[,c(1,(m*p[i,2]+2):(m*p[i,2]+1+p[i,3]))]
      }
	resid[,(m*(i-1)+1):(i*m)] = RR$resid           
}
	Sigmao            = t(resid)%*%resid/(T-m*(p[i,1]+p[i,2])-p[i,3])

	### put the estimates into the GVARData object 
	### c("Y","Uo","B","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
    	G = BoAoW2G(Bo,Ao,W,m,n,Pmax)
      Gs = diag(n*m)
	for (L in 1:Pmax) { Gs = Gs - G[,,L] }

      LH_g  = - T*m*n/2*log(2*pi) - T*m*n/2 + T/2*log(det(solve(Sigmao)))
      TwoN  = sum(AIC) + 2*sum(LH) - n * m*(m+1) + (n*m)*(n*m+1)   
      AIC_g = TwoN - 2*LH_g
      BIC_g = log(T)*TwoN/2 - 2*LH_g
 
	res$G	= G
	res$C = C
	res$Sigma = Sigmao
	res$r_npo = NA
	res$Ao    = Ao
	res$Bo    = Bo
      res$Co    = Co
 	res$mu    = solve(Gs)%*%C
      res$VAR_domestic = VAR_domestic
      res$AIC   = AIC
      res$BIC   = BIC
      res$LH    = LH
      res$LH_g  = LH_g
      res$AIC_g = AIC_g
      res$BIC_g = BIC_g
	return(res)
}




#' @export
BoAoW2G = function(Bo,Ao,W,m,n,p) {
  dim(Bo) = c(m,m,p,n)
  dim(Ao) = c(m,m,p,n)
  G = (1:(n*m*n*m*p))*0
  dim(G) = c(n*m,n*m,p)
  for (i in 1:n) {
      for (j in 1:n) { 
         for (L in 1:p)    G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L] = Ao[,,L,i]*W[i,j]
      }
  }
  for (i in 1:n) {
	for (L in 1:p)       G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]  = Bo[,,L,i]  
  }
  dim(G) = c(n*m,n*m,p)
return(G)  
}






#' @export
GW2BoAo = function(G,W) {
  n = dim(W)[1]
  m = dim(G)[1]/n 
 
  if (is.na(dim(G)[3])) dim(G) = c(m*n,m*n,1)
  p = dim(G)[3]
  Bo = (1:(m*m*p*n))*0
  dim(Bo) = c(m,m,p,n)
  Ao = Bo 
  for (i in 1:n) {
      for (j in 1:n) { 
         for (L in 1:p)   if (!j==i) Ao[,,L,i] = G[(1+(i-1)*m):(i*m),(1+(j-1)*m):(j*m),L]/W[i,j]
      }
  }
  for (i in 1:n) {
	for (L in 1:p)       Bo[,,L,i] = G[(1+(i-1)*m):(i*m),(1+(i-1)*m):(i*m),L]   
  }
  result = list(Bo,Ao)
  names(result) = c("Bo","Ao")
return(result)  
}



#' Summary of GVAR estimation results 
#'
#' This function sumerizes the estimation results of GVAR in its country VAR model
#'
#' @param  res  :a list of the output of GVARest
#' @examples 
#' Co = c(50,50,50,50,50,50,50,50);dim(Co) = c(2,1,4)
#' res_d = GVARData0(m=2,n=4,p=3,T=200,type="const",Co=Co)
#' res_e = GVARest0(res_d)
#' summary_GVAR(res_e)
#'
#' @export
summary_GVAR = function(res) {
  VAR_domestic = res$VAR_domestic
  n = length(VAR_domestic)
  for (i in 1:n) {
    print(summary(VAR_domestic[[i]]))
  }
}


#' Impulse Response Functions of GVAR 
#'
#' This function generates impulse response functions of an estimated GVAR 
#'
#' @param  res  : a list of the output of GVARest
#' @param  nstep: length of the impulse response functions
#' @param  comb : an mn vector specifying combined impulse such as global shocks, reginal shocks, or concerted actions.
#' @param  type : a list of the output of GVARest
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a matrix of (mn,mn,nstep) as the IRF colunms respesenting the impulse rows the responses. 
#' @examples 
#'
#' 
#' X1 = matrix(1,200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X3 = matrix(rnorm(200),200,1)
#' X4 = matrix(rnorm(200),200,1)
#' X  = cbind(X1,X2,X3,X4)
#' dim(X) = c(200,1,4)
#'
#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=2; p[,3]=1;
#' res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
#' summary_GVAR(res_e)  # the summary needs to be detailes commented
#' 
#' IRF = irf_GVAR(res_e,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_GVAR_CB1(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#' plott(IRF_CB,1,8) 
#'
#'
#' @export
irf_GVAR = function(res,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1")) {
	B 	= res$G
      neq 	= dim(B)[1]
	nvar	= dim(B)[2] 
	sigma = res$Sigma
      response <- array(0,dim=c(neq,nvar,nstep));
      response <- impulsdtrf_g(B,sigma,nstep,comb,irf=irf)
	return(response)
} 


#' Impulse Response Functions of GVAR 
#'
#' This function generates impulse response functions of an estimated GVAR with confidence bands 
#'
#' @param  res  : a list of the output of GVARest
#' @param  nstep: length of the impulse response functions
#' @param  comb : an mn vector specifying combined impulse such as global shocks, reginal shocks, or concerted actions.
#' @param  type : a list of the output of GVARest
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @param  runs : number of bootstrap runs to generate the confidence bands   
#' @param  conf :   
#' @return a matrix of (mn,mn,nstep,3) as the IRF colunms respesenting the impulse rows the responses. 
#' @examples 
#'
#' X1 = matrix(rnorm(200),200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X=cbind(X1,X2)
#' Co = c(1:(2*3*4))
#' dim(Co) = c(2,3,4)
#' Co[,1,] = Co[,1,]*100
#'
#' res_d = GVARData0(m=2,n=4,p=1,T=200,type="exog1",Co=Co,X=X)
#' max(abs(res_d$Y)) 
#' 
#' A = res_d$G[,,1]
#' eigen(A)$values
#'
#'
#' res_eee = GVARest0(res_d)
#' 
#' summary_GVAR(res_eee)
#' res_eee$Co
#' res_d$Co
#' res_eee$Ao 
#'
#' IRF = irf_GVAR(res_eee,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_GVAR_CB(res_eee,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#' plott(IRF_CB,1,8) 
#'
#'
#' @export
irf_GVAR_CB = function(res,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95)) {
	m = res$m
	n = res$n
      p = res$p
      T = dim(res$Y)[1]
      W = res$W
      Ao= res$Ao
      Bo= res$Bo
      Go= res$G
      Co= res$Co 
      type=res$type
      X   = res$X
      mu  = res$mu

	B 	= res$G
      neq 	= dim(B)[1]
	nvar	= dim(B)[2] 
	sigma = res$Sigmao
      response <- array(0,dim=c(neq,nvar,nstep,length(conf)+1))
     #response[,,,1] <- impulsdtrf_g(B,sigma,nstep,comb,irf=irf)
      response[,,,1] <- irf_GVAR(res,nstep,comb,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs))
      for (i in 1:runs) {
      	Uo_run = rnormSIGMA(T,sigma)
            res_run = GVARData0(m,n,p,T,W,r_npo=NA,Ao=NA,Bo=NA,Go,Co,Uo=Uo_run,Sigmao=NA,type,X,mu)
            res_e   = GVARest0(res_run)
            B_run   = res_e$G
            sigma_run = res_e$Sigmao
		#responseR[,,,i] <- impulsdtrf_g(B_run,sigma_run,nstep,comb,irf=irf)
            responseR[,,,i]  <- irf_GVAR(res_e,nstep,comb,irf)
	}

      for (tt in 1:(nstep) ) {
	for (i in 1:neq)           {
		for (j in 1:nvar)     {
			response[i,j,tt,-1] = quantile(responseR[i,j,tt,], conf)
			
	} 
	}
	}
	return(response)
} 


#' Impulse Response Functions of GVAR 
#'
#' This function generates impulse response functions of an estimated GVAR with confidence bands 
#'
#' @param  res  : a list of the output of GVARest
#' @param  nstep: length of the impulse response functions
#' @param  comb : an mn vector specifying combined impulse such as global shocks, reginal shocks, or concerted actions.
#' @param  type : a list of the output of GVARest
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol Cholezky decomposition, Chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse. 
#' @param  runs : number of bootstrap runs to generate the confidence bands   
#' @param  conf : a two componets vector of the tail probabilities of the confidence interval.  
#' @return a (mn,mn,nstep,3) array of the IRF with colunms respesenting the impulse rows the responses. 
#' @examples 
#'
#' 
#' X1 = matrix(rnorm(200),200,1)
#' X2 = matrix(rnorm(200),200,1)
#' X3 = matrix(rnorm(200),200,1)
#' X4 = matrix(rnorm(200),200,1)
#' X  = cbind(X1,X2,X3,X4)
#' dim(X) = c(200,1,4)
#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1;   p[,3]=1; p[2,2]=2;
#'
#' res_d = GVARData(m=2,n=4,p=p,T=200,type="exog1",X=X)
#' res_e = GVARest(res = res_d)
#' summary_GVAR(res_e)  
#'
#' IRF = irf_GVAR(res_e,nstep=10,comb=NA,irf="gen")
#' IRF_CB = irf_GVAR_CB1(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
#'
#' x11()
#' par(mfrow=c(3,3))
#' plott(IRF_CB,1,1) 
#' plott(IRF_CB,1,2) 
#' plott(IRF_CB,1,3) 
#' plott(IRF_CB,1,4) 
#' plott(IRF_CB,1,5) 
#' plott(IRF_CB,1,6) 
#' plott(IRF_CB,1,7) 
#' plott(IRF_CB,1,8) 
#' plott(IRF_CB,1,1) 
#'
#'
#' @export
irf_GVAR_CB1 =function(res,nstep,comb,irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95)) {
        m = res$m
        n = res$n
      p = res$p
      T = dim(res$Y)[1]
      W = res$W
      Ao= res$Ao
      Bo= res$Bo
      Go= res$G
      Co= res$Co 
      type=res$type
      X   = res$X
      mu  = res$mu

        B       = res$G
      neq       = dim(B)[1]
        nvar    = dim(B)[2] 
        sigma = res$Sigmao
      response <- array(0,dim=c(neq,nvar,nstep,length(conf)+1))
     #response[,,,1] <- impulsdtrf_g(B,sigma,nstep,comb,irf=irf)
      response[,,,1] <- irf_GVAR(res,nstep,comb,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs))
      for (i in 1:runs) {
        Uo_run = rnormSIGMA(T,sigma)
            res_run = GVARData(m,n,p,T,W,r_npo=NA,Ao,Bo,Co,Uo=Uo_run,Sigmao=NA,type,X,mu)
            res_e   = GVARest(res_run)
            B_run   = res_e$G
            sigma_run = res_e$Sigmao
                #responseR[,,,i] <- impulsdtrf_g(B_run,sigma_run,nstep,comb,irf=irf)
            responseR[,,,i]  <- irf_GVAR(res_e,nstep,comb,irf)
        }
        responseR[,,,1] = response[,,,1]
      for (tt in 1:(nstep) ) {
        for (i in 1:neq)           {
                for (j in 1:nvar)     {Ao
                        response[i,j,tt,-1] = quantile(responseR[i,j,tt,], conf)
                        
        } 
        }
        }
        return(response)
}



#' @export
plott = function(IRF_CB,i,j) {
  ylim = c(min(IRF_CB[i,j,,]),max(IRF_CB[i,j,,]))
  plot(IRF_CB[i,j,,1],type="l",ylim=ylim)
  lines(IRF_CB[i,j,,2],col="red")
  lines(IRF_CB[i,j,,3],col="red")
  lines(IRF_CB[i,j,,3]*0,col='black')
}




#' @export
SW2comb = function(SW,n,N,K)
{
	comb = matrix(0,n*N,N*n)
	SW = SW/sum(SW)
	for (i in 1:n) {
		comb[(i-1)*N+K,1] = SW[i]
	}
	return(comb)
}


############# GVAR_Selection ##################


#' Calculation of the model selection values for GVAR models  
#' 
#' 
#' @param  res  : a GVAR object obtained from GVARData or estimated from CIGVARest.
#' @param  L_V  : a two components vector containing the maxima of the domestic lag and the foreign lag, respectively. 
#' @param  I    : Index of the country under investigation. 

#' @return      : A matrix with different lag specifications and values of the model selection criteria.
#' @examples 
#' 
#' n = 4
#' p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1; p[2:3,2] = 2 
#' res_d = GVARData(m=2,n=4,p=p,T=4000,type="const")
#' 
#' 
#' I = 3
#' L_V = c(4,4)
#' res_d$p
#' GVARSelect = GVAR_Select(res=res_d,L_V=c(4,4),I=2)
#' GVARSelect[which.min(GVARSelect[,3]),]
#'
#' @export
GVAR_Select = function(res=res_e,L_V=L_v,I=i)  {
m = res$m
n = res$n
Y = as.matrix(res$Y)
T = dim(Y)[1]
type = res$type
XXX  = res$X

Wmat = res_d$W
Wnmat = kronecker(Wmat,diag(m))
FYI   = Y%*%Wnmat[,((I-1)*m+1):(I*m)]
Yi = as.matrix(Y[,((I-1)*m+1):(I*m)])
Criteria = matrix(0,L_V[1]*L_V[2],4)
idx = 0

for (l_d in 1: L_V[1] )   {
   for (l_f in 1:L_V[2] ) {
      idx = idx + 1
	X = embed(FYI,l_f+1)[,-(1:m)]
	XX = matrix(0,T,dim(X)[2])
	XX[(l_f+1):T,] = X 
      if (type=="none")  { type_holder = "exog0"; Co = (1:(m*(1+l_f*m))); dim(Co) = c(m,l_f*m+1); Co[,1]=0} 
      if (type=="const") { type_holder = "exog1"; Co = (1:(m*(1+l_f*m))); dim(Co) = c(m,l_f*m+1); Co[,1]=1} 


	if (l_d>=l_f) {
            if (type=="exog0") { type_holder = type; XX = cbind(XX,XXX[,,I]); Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=0} 
	      if (type=="exog1") { type_holder = type; XX = cbind(XX,XXX[,,I]); Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=1} 
		res_dd = VARData(n=m,p=l_d,T=T,Co=Co,type=type_holder,X=as.matrix(XX))
		res_dd$Y = Yi
	}   else  {
            XX = XX[(l_f-l_d+1):T,];
	   	if (type=="exog0") { type_holder = type; XX = cbind(XX,XXX[(l_f-l_d+1):T,,I]);Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=0 } 
	   	if (type=="exog1") { type_holder = type; XX = cbind(XX,XXX[(l_f-l_d+1):T,,I]);Co = matrix(1,m,1+dim(XX)[2]); Co[,1]=1 } 
		res_dd = VARData(n=m,p=l_d,T=(T-(l_f-l_d)),Co=Co,type=type_holder,X=XX)
		res_dd$Y = Yi[(l_f-l_d+1):T,]
	}
	res_e = VARest(res=res_dd)
	Criteria[idx,] = c (l_d,l_f,res_e$AIC,res_e$BIC)
    }
}
colnames(Criteria) = c("Domestic Lags","Foreign Lags","AIC","BIC")
return(Criteria)
}

#' @export
W2Wmat = function(W) {
 n = length(W)
 Wmat = matrix(0,n,n)
 for (i in 1:n ) {
    Wmat[i,-i] = W[-i]/sum(W[-i])
 }
 return(Wmat)
}



INVI = function(sigma,c,i) {
		n = dim(sigma)[1]
            sigmaout = matrix(0,n,n)
		cc = as.numeric(c>0)*(1:n)  
            sigmai = sigma    
            for (j in 1:n)         {
                   for (k in 1:n)  {
                       if ((j==i)|(k==i)) sigmai[j,k] = 0
            }
            }
            sigmai[i,i] = sigma[i,i]
		invsigmai = solve(sigmai[cc,cc])
            sigmaout[cc,cc] = invsigmai
		return(sigmaout)
}


