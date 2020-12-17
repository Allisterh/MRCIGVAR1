#' Data generating process of MRGVAR(m,n,p,S) 
#'
#' This function will generate data from an multi-regime staionary Global VAR(p) process and return an object of MRGVAR(m,n,p,S) process which contain the generated data and the used paramters.
#'
#' @param m     : number of variables in each a country/unit
#' @param n     : number of countries/units
#' @param p     : lag length
#' @param T     : number of observations. 
#' @param W     : (n x n) weighting matrix. w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' @param SESVI : n vector of indices of the switching variables across n countries 
#' @param TH    : (n x S-1) matrix of threshold values 
#' @param Go    : (mn,mn,p,S) array of the MRGVAR(m,n,p,S) coefficients. G is contructed from Bo, Ao and W.
#' @param Ao    : (m, m, p, n, S) array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
#' @param Bo    : (m, m, p, n, S) array collecting the n country VAR(p) coefficients. 
#' @param Co    : (m , k+1, n, S) array collecting the coeeficients of the deterministic components of the n countries.
#' @param Uo    : (T, mn, S) arrary of the temporally independent innovation processes
#' @param Sigmao: (mn, mn, S) array of the covariance matrix of the MRGVAR(m,n,p,S)
#' @param SV    : exogeneous switching variables
#' @param type	: deterministic component "const" and "none" are two options 
#' @param X	: (T x k) matrix of exogeneous variables.
#' @return      A MRGVAR object containing the generated data, the parameters and the input exogeous variables. res = list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d)
#' @field Y     : T x nm simulated data via of the MRGVAR(m,n,p,S). 
#' @field X     : (T x k) matrix of exogeneous variables.
#' @field Uo    : (T, mn,S) array of the simulated innovations of the MRGVAR(m,n,p,S) 
#' @field Go    : (mn, mn, p, S)  array of the MRGVAR(m,n,p,S) coefficients. G is contructed from Bo, Ao and W.
#' @field C     : (nm, (k+1),S) array containing the coefficients of the deterministic components.   
#' @field Sigmao: (mn, mn, S) array of the covariance matrix of the MRGVAR(m,n,p,S)
#' @field TH    : (n x S-1) matrix of threshold values
#' @field St    : simulated time path of states/regimes
#' @field SESVI : index if the switching variable
#' @field check : maximum of the data for checking the stationarity
#' @examples 
#' res_d <- MRGVARData2(m=2,n=2,p=1,T=1000,S=2,SESVI=c(1,3))
#' res_e = MRGVARest0(res)
#' res_d = MRGVARData2(m=2,n=2,p=2,T=1000,S=2,SESVI=c(1,3),type="const")
#' res_e = MRGVARest0(res)
#' res_d = MRGVARData2(m=2,n=3,p=3,T=100,S=2,SESVI=c(1,3,5))
#' res_e = MRGVARest0(res)
#' @export
MRGVARData2=function(m,n,p,T,S,W=NA,SESVI=NA,TH=NA,Go=NA,Ao=NA,Bo=NA,Sigmao=NA,Uo=NA,SV=NA,type=NA,Co=NA,X=NA,Yo=NA,d=NA) {
### m     : number of variables in a country
### n     : number of countries
### p     : lag length
### T     : number of observations
### S     : number of regimes
### SESVI : n vector of indeces of the switching variable in the endogeneous variables Y for the case of self-excited threshhold model.
###         (m,n,p,T,S,SESVI) are parameters which must be provided.
### TH    : n x S-1  matrix of threshold values of n countries.
### Go    : mn x mn x p x S GVAR(m,n,p) coefficients matrix of S different regimes 
### Sigmao: mn x mn x S array of the covariance matrix of the GVAR(m,n,p) in S different regimes
### Uo    : a T x mn x S  array of the temorally independent innovation processes
### X	    : a T x k matrix of exogenous veriables which are common to all countries 
###         (TH,Go,Sigmao,Uo,SV) if not provided, they will be generated randomly. 
### type  : deterministic component "const" and "none" are two options 
### Co    : if type = "const" mu is an m x n x S array of the intercepts of the time series in the different regimes
###
###
### output:
### c("Y","Uo","Go","Sigmao","TH","St","sv","SESVI","check")
###
### Y     : T x mn data matrix of the simulated data via of the MSGVAR(m,n,p,S) 
### Uo    : an T x mn x S   array of the temorally independent innovation processes of the MSGVAR(m,n,p,S)
### Go    : mn x mn x p x S array collecting the MSGVAR(m,n,p,S) coefficients in S different states 
### Sigmao: mn x mn x S array collecting the covariance matrices of the simulated MSGVAR(m,n,p,S) in S different states
### TH    : S-1 vector of thresholds
### St    : simulated time path of states/regimes
### SESVI : index if the switching variable
### check : maximum of the data to check stationarity
###
### Remarks: The states of each country at each time step is governed by the switching variables of each country and hence there is a large numbers
###          of possible combinations of gegimes. 
###          The coefficients matrix can be constructed from Go, given a regime combination. The stationarity of the GVAR(m,n,p)
###          at each time step is not guaranteed. But this is not relevant. (This is a open question.)
###          
###                   
  check = c(1:(S+1))*0

  if (missing(TH)) 	{TH = NA} 
  if (missing(W))  	{W = NA}
  if (missing(SV))      {SV = NA}
  if (missing(Sigmao)) 	{Sigmao=NA}
  if (missing(Uo))      {Uo = NA}
  if (missing(type))	{type = NA}
  if (missing(Co))      {Co = NA}
  if (missing(Go))      {Go = NA}
  if (missing(Bo))      {Bo = NA}
  if (missing(X))       {X = NA}
  if (missing(SV))      {SV = NA}
  if (missing(d))       { d = NA}
  if (missing(Yo))      {Yo = NA}
  if (anyNA(d))  { d = 1 }
  P  = max(p,d)
 

 if (anyNA(TH)) {
  	TH = matrix(runif(S-1)*n,S-1,n)
      for (i in 1:n)  {
         THh = runif(S-1)-0.5
         THh = THh[order(THh)]
         TH[,i] = THh
      }
  }

  if (anyNA(W)) {
  	W = matrix(((1:(n*n))/(1:(n*n))*1/(n-1)),n,n)
      for (i in 1:n) W[i,i]= 0
  }
  
  if (anyNA(Sigmao)) {
  	Sigmao = (1:(m*n*m*n*S))*0
      dim(Sigmao) = c(m*n,m*n,S)
      for (i in 1:S) {
      	GVARD = GVARData0(m,n,p,T) 
      	Sigmao[,,i] = GVARD$Sigmao
      }    

  }

  if (anyNA(Uo)) {
  	Uo = (1:(T*n*m*S))*0
      dim(Uo) = c(T,m*n,S)
      for (s in 1:S) Uo[,,s]= rnormSIGMA(T,Sigmao[,,s])
  }

  Ct = Uo*0;

  if (anyNA(type)) {
  		type = "none"
  }

  if (type=="none") {
	CDC = c(1:(m*n*S))*0; dim(CDC)=c(m,1,n,S)
	Co = CDC
	GDC = Co; dim(GDC) = c(n*m,1,S)
      Ct  = Uo*0
  }
  if (type=="const") {
	CDC = rnorm(m*n*S); dim(CDC)=c(m,1,n,S)
      if (anyNA(Co)) Co = CDC
	GDC = Co; dim(GDC) = c(n*m,1,S)
      CoV = as.vector(Co)
	for (s in 1:S) {
            CoV = as.vector(Co[,,,s])
		for (i in 1:n)  Ct[,,s] = matrix(1,T,1)%*%t(CoV)
      }   

  }

  if (type=="exog1")  {
	k 	   = dim(X)[2]+1
     	CDC      = rnorm(m*k*n*S)
	dim(CDC) = c(m,k,n,S) 

	if (anyNA(Co)) Co = CDC
	GDC	   = Co
      dim(GDC) = c(n*m,k,S);
      CoV = matrix(0,k,m*n);
     

      for (s in 1:S) {
            for (i in 1:k)  CoV[i,]  = as.vector(Co[,i,,s])
		Ct[,,s] = cbind(matrix(1,dim(X)[1],1),X)%*%CoV
      }      
  }   

  if (type=="exog0")  {
	k 	   = dim(X)[2]+1
     	CDC      = rnorm(m*k*n*S)
	dim(CDC) = c(m,k,n,S)
      CDC[,1,,] = 0 
    	if (anyNA(Co)) Co = CDC
	GDC	   = Co
      dim(GDC) = c(n*m,k,S);
      CoV = matrix(0,k,m*n);

      for (s in 1:S) {
  		for (i in 1:k)  CoV[i,]  = as.vector(Co[,i,,s])
		Ct[,,s] = X%*%CoV[-1,]
      }      
  }   
 

  if (anyNA(Go) & anyNA(Bo)) {
  	Go 	  = (1:(m*n*m*n*p*S))*0
      dim(Go) = c(m*n,m*n,p,S)
	Bo 	  = (1:(m*m*p*n*S))*0
	dim(Bo) = c(m,m,p,n,S)
	Ao      = Bo*0
      for (s in 1:S) {
            GVARD     = GVARData0(m,n,p,T,W=W)    #(m,n,p,T,W,r_npo,Ao,Go,Uo,Sigmao)  !!!!W must be the same for all regimes
      	Go[,,,s]  = GVARD$G
            check[s]  = max(abs(GVARD$Y))
		Ao[,,,,s] = GVARD$Ao 
		Bo[,,,,s] = GVARD$Bo
      }            
  }

  if (!anyNA(Go) & anyNA(Bo))        {
      Bo 	  = (1:(m*m*p*n*S))*0
	dim(Bo) = c(m,m,p,n,S)
	Ao      = Bo*0
      for (s in 1:S)  {
      	Bo[,,,,s] = GW2BoAo(Go[,,,s],W)$Bo
		Ao[,,,,s] = GW2BoAo(Go[,,,s],W)$Ao
      }
  }

  if (anyNA(Go) & !anyNA(Bo)) {
      Go 	  = (1:(m*n*m*n*p*S))*0
      dim(Go) = c(m*n,m*n,p,S)

  	for (s in 1:S )Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,p)
  }
	
  
  St = matrix(0,T,n)
  Y = Uo[,,1]
  if (!anyNA(Yo)) Y[1:P,] = Yo
  resid = Y*0
  Yo    = Y[1:P,]
  if (anyNA(SV)) {sv = Y[,SESVI]; SV=NA} else { sv = SV }
  ss = matrix(0,n,1)

  for ( tt in (P+1):T )  { 
      #if (anyNA(SV)) sv = Y[,SESVI]
      svt = sv[tt-d,];
      Bt=c(1:(m*n*m*n*p))*0
      dim(Bt) = c(n*m,n*m,p)
      for (i in 1:n) {
      	ss[i] = which.min(abs(TH[,i]-svt[i]))
      	if (svt[i]>TH[ss[i],i]) { ss[i] = ss[i]+1 }
            Bt[(1+(i-1)*m):(i*m),,]     = Go[(1+(i-1)*m):(i*m),,,ss[i]]  
	     	Y[tt,(1+(i-1)*m):(i*m)]     = Uo[tt,(1+(i-1)*m):(i*m),ss[i]] + Ct[tt,(1+(i-1)*m):(i*m),ss[i]]
		resid[tt,(1+(i-1)*m):(i*m)] = Uo[tt,(1+(i-1)*m):(i*m),ss[i]]	
	}
      for (L in 1:p)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bt[,,L])
      St[tt,] = ss 
  }
  check[S+1] = max(abs(Y))
  SigmaS = NA;   			# SigmaS is a (nmS x nmS) matrix used to construct state-dependent covariance matrices
  result=list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d)
  names(result) = c("Y","X","Uo","resid","Go","GDC","Co","Sigmao","TH","St","SV","SESVI","Ao","Bo","check","type","m","n","p","S","W","SigmaS","Yo","d")
  return(result)
}


#' Estimation of MRGVAR(m,n,p,S) 
#'
#' This function estimate the unknown parameters of a specified MRGVAR(m,n,p,S) model based on provided data.
#'
#' @param  res  :a list containing the components as the output of MRGVARData2 including at least: m,n,p,S,type,Y and optionally X. 
#' @return res  :a with the same structure as the input but filled with estimated paramters.
#' @examples 
#' Co = c(1,1,1,1,1,1,10,10,10,10,10,10)
#' dim(Co) = c(m,1,n,S)
#' res = MRGVARData2(m=2,n=3,p=1,T=200,S=2,SESVI=c(1,3,5),type="const",Co=Co)
#' res_est = MRGVARest0(res)
#' @export
MRGVARest0 = function(res) {
m	= res$m
n	= res$n
p	= res$p
Y	= res$Y
X     = res$X
W 	= res$W
type  = res$type
TH    = res$TH
SESVI = res$SESVI
Go	= res$Go
S     = res$S
Ao	= res$Ao
Bo    = res$Bo
Co    = res$Co
GDC   = res$GDC
Sigmao= res$Sigmao*0
SigmaS = (1:(n*m*S*n*m*S))*0
dim(SigmaS) = c(n*m*S,n*m*S)
VAR_domestic = list()

T	= dim(Y)[1]
FY	= Y%*%t(W%x%diag(m))
resid = (1:(T*m*n*S))*0; dim(resid) = c(T,m,n,S)
##PS  = (1:(T*n*S))*0;   dim(PS)    = c(T,n,S)

###   CC is the single country coefficients of determinsitic components i.e. the foreign variables and the deterministic variables
###   It is m x (m*p+k+1) where m*p is the number of foreign variables and k is the number of deterministic variables. with zeros for "none"  

if (type=="none" |type=="const")  k = 1;
if (type=="exog0"|type=="exog1")  k = dim(X)[2]+1;
CC  =  c(1:(m*(p*m+k)*S))*0
dim(CC) = c(m,m*p+k,S) 

#res_MSVAR = MSVARData(n=m,p=p,T=T,S=S,Co=CC,TM=TM,type="exog",X=Z) 
#res_MSVAR = MSVARData(n=m,p=p,T=T,S=S)

EMCVG       = c(1:n)*0
for (i in 1:n) {
      
	Z         	=  matrix(0,T,p*m)
      FYp	    	=  embed(FY[,(m*(i-1)+1):(i*m)],(p+1))
	Z[(p+1):T,] =  FYp[,(m+1):(p*m+m)]
	if (type=="const")  Z = Z
	if (type=="exog0")  Z = cbind(Z,X)
      if (type=="exog1")  Z = cbind(Z,X)
	if (type=="none"|type=="exog0")  type_holder = "exog0"
	if (type=="const"|type=="exog1") type_holder = "exog1"

	### put data into a MRVARData object in order to replace the data and order parameters
      res_MRVAR 		=  MRVARData4(n=m,p=p,T=T,S=S,Co=CC,SESVI=SESVI[i],type=type_holder,X=Z)

     	#res_MRVAR 		=  MRVARData4(n=m,p=p,T=T,S=S,Co=CC,SESVI=(SESVI[i]-SESVI[i]%/%m*m),type=type_holder,X=Z)
	res_MRVAR$Y 	=  Y[,(m*(i-1)+1):(i*m)]
	res_MRVAR$X		=  Z
      res_MRVAR$sv      =  Y[,SESVI[i]] 
      res_MRVAR$TH      =  TH[,i]
      RR			=  MRVARest1(res_MRVAR)
      VAR_domestic[[i]] =  RR
      #EMCVG[i] = RR$MinColSumsTM

	Bo[,,,i,]		=  RR$res$Bo
		Ai		=  RR$res$Co[,2:(m*p+1),]
	     #Ai		=  t(Ai)
       dim(Ai)		=  c(m,m,p,S)
	Ao[,,,i,]		=  Ai

      if (type=="none")    Co[,,i,]          		= RR$res$Co[,1,]
	if (type=="const")   Co[,,i,]			 	= RR$res$Co[,1,]
	if (type=="exog1")   Co[,,i,]			 	= RR$res$Co[,c(1,(m*p+2):(m*p+2+dim(X)[2]-1)),]
	if (type=="exog0")   Co[,,i,] 			= RR$res$Co[,c(1,(m*p+2):(m*p+2+dim(X)[2]-1)),]

      Sigmao[(m*(i-1)+1):(i*m),(m*(i-1)+1):(i*m),] 	= RR$res$Sigmao
	resid[,,i,]                           		= RR$res$resid
      GDC = Co; 
      dim(GDC) = c(m*n,k,S)
	##PS[,i,]                                    	= RR$PS
	#resid[,(m*(i-1)+1):(i*m)] = RR$res$resid
	           
}

      for (s in 1:S) {
		for (i in 1:n) {
			for (j in 1:n)  Sigmao[(m*(i-1)+1):(i*m),(m*(j-1)+1):(j*m),s] = t(resid[,,i,s])%*%(resid[,,j,s])/(dim(Z)[1]-dim(Z)[2])

		}
	}
 
      for (s in 1:S)  { 
          for (ss in 1:S)  {
               for (i in 1:n)  {
                   for (j in 1:n) SigmaS[ ((s-1)*n*m+(i-1)*m+1):((s-1)*n*m+i*m), ((ss-1)*n*m+(j-1)*m+1):((ss-1)*n*m+j*m)] = t(resid[,,i,s])%*%(resid[,,j,ss])/sum( (!resid[,1,i,s]==0)*(!resid[,1,j,ss]==0))         
		   }
	    }
	}
	
	### put the estimates into the GVARData0 object 
	### c("Y","Uo","B","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
    	for (s in 1:S) {
		BBo = Bo[,,,,s]
		AAo = Ao[,,,,s]
     		dim(BBo) = c(m,m,p,n)
		dim(AAo) = c(m,m,p,n)
		Go[,,,s] = BoAoW2G(BBo,AAo,W,m,n,p) 
	}
      
      
      dim(resid) = c(T,m*n,S)
	res$Go	= Go
	res$C       = CC
	res$Sigmao  = Sigmao
	res$r_npo   = NA
	res$Ao      = Ao
	res$Bo      = Bo
      res$Co      = Co
      res$GDC     = GDC
      res$SigmaS  = SigmaS	
      res$resid   = resid
      res$VAR_domestic = VAR_domestic
	return(res)
}

#' Data generating process of MRGVAR(m,n,p,S) 
#'
#' This function will generate data from an multi-regime staionary Global VAR(p) process and return an object of MRGVAR(m,n,p,S) process which contain the generated data and the used paramters.
#'
#' @param m     : number of variables in each a country/unit
#' @param n     : number of countries/units
#' @param p     : (n, 3, S) array, each raw contains the lag length of the domestic variables, the foreign variables and the number of the exogeneous variables
#' @param T     : number of observations. 
#' @param W     : (n x n) weighting matrix. w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' @param SESVI : n vector of indices of the switching variables across n countries 
#' @param TH    : (n x S-1) matrix of threshold values 
#' @param Go    : (mn,mn,p,S) array of the MRGVAR(m,n,p,S) coefficients. G is contructed from Bo, Ao and W.
#' @param Ao    : (m, m, p, n, S) array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
#' @param Bo    : (m, m, p, n, S) array collecting the n country VAR(p) coefficients. 
#' @param Co    : (m , k+1, n, S) array collecting the coeeficients of the deterministic components of the n countries.
#' @param Uo    : (T, mn, S) arrary of the temporally independent innovation processes
#' @param Sigmao: (mn, mn, S) array of the covariance matrix of the MRGVAR(m,n,p,S)
#' @param SV    : exogeneous switching variables
#' @param type	: deterministic component "const" and "none" are two options 
#' @param X	: (T x k) matrix of exogeneous variables.
#' @return      A MRGVAR object containing the generated data, the parameters and the input exogeous variables. res = list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d)
#' @field Y     : T x nm simulated data via of the MRGVAR(m,n,p,S). 
#' @field X     : (T x k) matrix of exogeneous variables.
#' @field Uo    : (T, mn,S) array of the simulated innovations of the MRGVAR(m,n,p,S) 
#' @field Go    : (mn, mn, p, S)  array of the MRGVAR(m,n,p,S) coefficients. G is contructed from Bo, Ao and W.
#' @field C     : (nm, (k+1),S) array containing the coefficients of the deterministic components.   
#' @field Sigmao: (mn, mn, S) array of the covariance matrix of the MRGVAR(m,n,p,S)
#' @field TH    : (n x S-1) matrix of threshold values
#' @field St    : simulated time path of states/regimes
#' @field SESVI : index if the switching variable
#' @field check : maximum of the data for checking the stationarity
#'
#' @examples 
#' ## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries 
#' p = rep(1,12); dim(p) = c(2,3,2)
#' p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=2000,S=2,SESVI=c(1,3))   
#' max(res_d$Y)
#' ### estimation of the MRGVAR model
#' res_e = MRGVARest(res=res_d) 
#' summary_MRGVAR(res_e)
#'
#'
#' @export
MRGVARData=function(m,n,p,T,S,W=NA,SESVI=NA,TH=NA,Go=NA,Ao=NA,Bo=NA,Sigmao=NA,Uo=NA,SV=NA,type=NA,Co=NA,X=NA,Yo=NA,d=NA) {
### m     : number of variables in a country
### n     : number of countries
### p     : lag length n x 3 x S array collecting lag length for each country with respect to domestic and foreigen variable for each state. The last column specifies the number of exogenous variabls for each state.
### T     : number of observations
### S     : number of regimes
### SESVI : n vector of indeces of the switching variable in the endogeneous variables Y for the case of self-excited threshhold model.
###         (m,n,p,T,S,SESVI) are parameters which must be provided.
### TH    : n x S-1  matrix of threshold values of n countries.
### Go    : mn x mn x p x S GVAR(m,n,p) coefficients matrix of S different regimes 
### Sigmao: mn x mn x S array of the covariance matrix of the GVAR(m,n,p) in S different regimes
### Uo    : a T x mn x S  array of the temorally independent innovation processes
### X	    : a T x k x n x S array of exogenous veriables which may be common/different to all countries and different across all states 
###         (TH,Go,Sigmao,Uo,SV) if not provided, they will be generated randomly. 
### type  : deterministic component "const" and "none" "exog0" "exog1" are foure options 
### Co    : if type = "const" mu is an m x n x S array of the intercepts of the time series in the different regimes
###
###
### output:
### c("Y","Uo","Go","Sigmao","TH","St","sv","SESVI","check")
###
### Y     : T x mn data matrix of the simulated data via of the MSGVAR(m,n,p,S) 
### Uo    : an T x mn x S   array of the temorally independent innovation processes of the MSGVAR(m,n,p,S)
### Go    : mn x mn x p x S array collecting the MSGVAR(m,n,p,S) coefficients in S different states 
### Sigmao: mn x mn x S array collecting the covariance matrices of the simulated MSGVAR(m,n,p,S) in S different states
### TH    : S-1 vector of thresholds
### St    : simulated time path of states/regimes
### SESVI : index if the switching variable
### check : maximum of the data to check stationarity
###
### Remarks: The states of each country at each time step is governed by the switching variables of each country and hence there is a large numbers
###          of possible combinations of gegimes. 
###          The coefficients matrix can be constructed from Go, given a regime combination. The stationarity of the GVAR(m,n,p)
###          at each time step is not guaranteed. But this is not relevant. (This is a open question.)
###          
###                   
  check = c(1:(S+1))*0

  if (missing(TH)) 	{TH = NA} 
  if (missing(W))  	{W = NA}
  if (missing(SV))      {SV = NA}
  if (missing(Sigmao)) 	{Sigmao=NA}
  if (missing(Uo))      {Uo = NA}
  if (missing(type))	{type = NA}
  if (missing(Co))      {Co = NA}
  if (missing(Go))      {Go = NA}
  if (missing(Bo))      {Bo = NA}
  if (missing(X))       {X = NA}
  if (missing(SV))      {SV = NA}
  if (missing(d))       { d = NA}
  if (missing(Yo))      {Yo = NA}
  if (anyNA(d))  { d = 1 }
  P  = max(p[,1:2,],d)
  Pmax = max(p[,1:2,])


 if (anyNA(TH)) {
  	TH = matrix(runif(S-1)*n,S-1,n)
      for (i in 1:n)  {
         THh = runif(S-1)-0.5
         THh = THh[order(THh)]
         TH[,i] = THh
      }
  }

  if (anyNA(W)) {
  	W = matrix(((1:(n*n))/(1:(n*n))*1/(n-1)),n,n)
      for (i in 1:n) W[i,i]= 0
  }
  
  if (anyNA(Sigmao)) {
  	Sigmao = (1:(m*n*m*n*S))*0
      dim(Sigmao) = c(m*n,m*n,S)
      for (s in 1:S) {
      	GVARD = GVARData0(m,n,max(p[,,s]),T) 
      	Sigmao[,,s] = GVARD$Sigmao
      }    

  }

  if (anyNA(Uo)) {
  	Uo = (1:(T*n*m*S))*0
      dim(Uo) = c(T,m*n,S)
      for (s in 1:S) Uo[,,s]= rnormSIGMA(T,Sigmao[,,s])
  }

  Ct = Uo*0;

  if (anyNA(type)) {
  		type = "none"
  }

  if (type=="none") {
	CDC = c(1:(m*n*S))*0; dim(CDC)=c(m,1,n,S)
	Co = CDC
	GDC = Co; dim(GDC) = c(n*m,1,S)
      Ct  = Uo*0
  }
  if (type=="const") {
	CDC = rnorm(m*n*S); dim(CDC)=c(m,1,n,S)
      if (anyNA(Co)) Co = CDC
	GDC = Co; dim(GDC) = c(n*m,1,S)
      CoV = as.vector(Co)
	for (s in 1:S) {
            CoV = as.vector(Co[,,,s])
		for (i in 1:n)  Ct[,,s] = matrix(1,T,1)%*%t(CoV)
      }   

  }
  
  if (type=="exog1")  {
	k 	   = dim(X)[2]+1  # wrong
     	CDC      = rnorm(m*k*n*S)
	dim(CDC) = c(m,k,n,S) 

	if (anyNA(Co)) Co = CDC
	GDC	   = Co
      dim(GDC) = c(n*m,k,S);
      CoV = matrix(0,k,m*n);
     
      for (s in 1:S) {
         for (i in 1:n) {
            if ( p[i,3,s]<k-1) Co[,(p[i,3,s]+1):k,i,s] = 0
		for (j in 1:k ) CoV[j,]  = as.vector(Co[,j,,s])
		Ct[,,s] = cbind(matrix(1,dim(X)[1],1),X[,,i,s])%*%CoV
	   }
      }      
  }   

  if (type=="exog0")  {
	k 	   = dim(X)[2]+1   # wrong
     	CDC      = rnorm(m*k*n*S)
	dim(CDC) = c(m,k,n,S)
      CDC[,1,,] = 0 
    	if (anyNA(Co)) Co = CDC
	GDC	   = Co
      dim(GDC) = c(n*m,k,S);
      CoV = matrix(0,k,m*n);

      for (s in 1:S) {
         for (i in 1:n) {
            if ( p[i,3,s]<k-1) Co[,(p[i,3,s]+2):k,i,s] = 0
		for (j in 1:k ) CoV[j,]  = as.vector(Co[,j,,s])
		Ct[,,s] = cbind(matrix(1,dim(X)[1],1),X[,,i,s])%*%CoV
	   }
      }        
  }   
 

  if (anyNA(Go) & anyNA(Bo)) {
  	Go 	  = (1:(m*n*m*n*Pmax*S))*0
      dim(Go) = c(m*n,m*n,Pmax,S)
	Bo 	  = (1:(m*m*Pmax*n*S))*0
	dim(Bo) = c(m,m,Pmax,n,S)
	Ao      = Bo*0
      for (s in 1:S) {
            GVARD     = GVARData(m,n,p[,,s],T,W=W)    #(m,n,p,T,W,r_npo,Ao,Go,Uo,Sigmao)  !!!!W must be the same for all regimes
      	Go[,,1:max(p[,1:2,s]),s]  = GVARD$G
            check[s]  = max(abs(GVARD$Y))
		Ao[,,1:max(p[,1:2,s]),,s] = GVARD$Ao 
		Bo[,,1:max(p[,1:2,s]),,s] = GVARD$Bo
      }            
  }

  if (!anyNA(Go) & anyNA(Bo))        {
      Bo 	  = (1:(m*m*Pmax*n*S))*0
	dim(Bo) = c(m,m,Pmax,n,S)
	Ao      = Bo*0
      for (s in 1:S)  {
      	Bo[,,,,s] = GW2BoAo(Go[,,,s],W)$Bo
		Ao[,,,,s] = GW2BoAo(Go[,,,s],W)$Ao
      }
  }

  if (anyNA(Go) & !anyNA(Bo)) {
      Go 	  = (1:(m*n*m*n*Pmax*S))*0
      dim(Go) = c(m*n,m*n,Pmax,S)

  	for (s in 1:S ) Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,Pmax)
  }
	
  
  St = matrix(0,T,n)
  Y = Uo[,,1]
  if (!anyNA(Yo)) Y[1:P,] = Yo
  resid = Y*0
  Yo    = Y[1:P,]
  if (anyNA(SV)) {sv = Y[,SESVI]; SV=NA} else { sv = SV }
  ss = matrix(0,n,1)

  for ( tt in (P+1):T )  { 
      if (anyNA(SV)) {sv = Y[,SESVI]; SV=NA} else { sv = SV }
      svt = sv[tt-d,];
      Bt=c(1:(m*n*m*n*Pmax))*0
      dim(Bt) = c(n*m,n*m,Pmax)
      for (i in 1:n) {
      	ss[i] = which.min(abs(TH[,i]-svt[i]))
      	if (svt[i]>TH[ss[i],i]) { ss[i] = ss[i]+1 }
            Bt[(1+(i-1)*m):(i*m),,]     = Go[(1+(i-1)*m):(i*m),,,ss[i]]  
	     	Y[tt,(1+(i-1)*m):(i*m)]     = Uo[tt,(1+(i-1)*m):(i*m),ss[i]] + Ct[tt,(1+(i-1)*m):(i*m),ss[i]]
		resid[tt,(1+(i-1)*m):(i*m)] = Uo[tt,(1+(i-1)*m):(i*m),ss[i]]	
	}
      for (L in 1:Pmax)    Y[tt,] = Y[tt,] + Y[tt-L,]%*%t(Bt[,,L])
      St[tt,] = ss 
  }
  check[S+1] = max(abs(Y))
  SigmaS = NA;   			# SigmaS is a (nmS x nmS) matrix used to construct state-dependent covariance matrices
  result=list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d)
  names(result) = c("Y","X","Uo","resid","Go","GDC","Co","Sigmao","TH","St","SV","SESVI","Ao","Bo","check","type","m","n","p","S","W","SigmaS","Yo","d")
  return(result)
}


#' Estimation of MRGVAR(m,n,p,S) 
#'
#' This function estimate the unknown parameters of a specified MRGVAR(m,n,p,S) model based on provided data.
#'
#' @param  res  :a MRGVAR object which is an output of MRGVARData.  
#' @return res  :a MRGVAR object with estimated paramters.
#' @examples 
#'
#' ### case of n = 2, m = 2, S = 2
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:2,,]/2; p[1,1,2]=2; p[2,2,2] = 2
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_dd <- MRGVARData(m=2,n=2,p=p,TH=TH,T=2000,S=2,SESVI=c(1,3))    ## m: number of variables, n: number of countries 
#' max(res_dd$Y)
#' res_dd$TH  
#' ### estimation of the MRGBAR model
#' res_ee = MRGVARest(res=res_dd) 
#' summary_MRGVAR(res_ee)
#'
#' @export
MRGVARest = function(res) {
m       = res$m
n       = res$n
p       = res$p
Y       = as.matrix(res$Y)
X     = res$X
W       = res$W
type  = res$type
TH    = res$TH
SESVI = res$SESVI
Go      = res$Go
S     = res$S
Ao      = res$Ao
Bo    = res$Bo
Co    = res$Co
GDC   = res$GDC
Sigmao= res$Sigmao*0
d     = res$d
P     = max(p,d) 
Pmax  = max(p[,1:2,])
SigmaS = (1:(n*m*S*n*m*S))*0
dim(SigmaS) = c(n*m*S,n*m*S)
VAR_domestic = list()

SESVIi = SESVI/SESVI+SESVI[1]-1

kmax  = max(p[,3,])

T       = dim(Y)[1]
FY      = Y%*%t(W%x%diag(m))
resid = (1:(T*m*n*S))*0; dim(resid) = c(T,m,n,S)
##PS  = (1:(T*n*S))*0;   dim(PS)    = c(T,n,S)

###   CC is the single country coefficients of determinsitic components i.e. the foreign variables and the deterministic variables
###   It is m x (m*p+k+1) where m*p is the number of foreign variables and k is the number of deterministic variables. with zeros for "none"  

if (type=="none" |type=="const")  k = 1;
if (type=="exog0"|type=="exog1")  k = dim(X)[2]+1;
CC  =  c(1:(m*(Pmax*m+k)*S))*0
dim(CC) = c(m,m*Pmax+k,S) 

#res_MSVAR = MSVARData(n=m,p=p,T=T,S=S,Co=CC,TM=TM,type="exog",X=Z) 
#res_MSVAR = MSVARData(n=m,p=p,T=T,S=S)

EMCVG       = c(1:n)*0
for (i in 1:n) {
      Pi          = max(p[i,1:2,])
        CCi             = c(1:(m*(Pi*m+k)*S))*0
        dim(CCi)        = c(m,m*Pi+k,S) 
        Z               = matrix(0,T,Pi*m)
      FYp               = embed(FY[,(m*(i-1)+1):(i*m)],(Pi+1))
        Z[(Pi+1):T,]= FYp[,(m+1):(Pi*m+m)]
      kz          = dim(Z)[2]
        if (type=="none")   {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz,S)); for (s in 1:S) ZZ[,1:kz,s]=Z[,,s]}
        if (type=="const")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz,S)); for (s in 1:S) ZZ[,1:kz,s]=Z[,,s]}
        if (type=="exog0")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz+kmax,S)); for (s in 1:S) ZZ[,1:(kz+p[i,3,s]),s] = cbind(Z[,,s],X[,1:p[i,3,s],i,s])}
      if (type=="exog1")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz+kmax,S)); for (s in 1:S) ZZ[,1:(kz+p[i,3,s]),s] = cbind(Z[,,s],X[,1:p[i,3,s],i,s])}
        if (type=="none"|type=="exog0")  type_holder = "exog0"
        if (type=="const"|type=="exog1") type_holder = "exog1"
      XX = (1:(T*dim(ZZ)[2]*S))*0; dim(XX) = c(T,dim(ZZ)[2],S)
      if (!anyNA(X)) for (s in 1:S) XX[,1:(m*p[i,2,s]+p[i,3,s]),s] = ZZ[,c(1:(m*p[i,2,s]),(m*Pi+1):(m*Pi+p[i,3,s])),s]
      if (anyNA(X))  for (s in 1:S) XX[,1:(m*p[i,2,s]),s] = ZZ[,c(1:(m*p[i,2,s])),s]

        ### put data into a MRVARData object in order to replace the data and order parameters
      pp = t(p[i,,]); pp[,2] = pp[,2]*m; ppp = pp[,1:2]; ppp[,2] = pp[,2]+pp[,3]
      res_MRVAR                 =  MRVARData(n=m,p=ppp,T=T,S=S,Co=NA,SESVI=SESVIi[i],type=type_holder,X=XX)

        #res_MRVAR              =  MRVARData4(n=m,p=p,T=T,S=S,Co=CC,SESVI=(SESVI[i]-SESVI[i]%/%m*m),type=type_holder,X=Z)
        res_MRVAR$Y     =  Y[,(m*(i-1)+1):(i*m)]
        res_MRVAR$X             =  XX
      res_MRVAR$SV      =  res_MRVAR$Y[,SESVIi[i]] 
      res_MRVAR$TH      =  TH[,i]
      RR                        =  MRVARest(res=res_MRVAR)
      VAR_domestic[[i]] =  RR
      #EMCVG[i] = RR$MinColSumsTM

        for (s in 1:S)  { 
                Bo[,,1:p[i,1,s],i,s]            =  RR$res$Bo[,,1:p[i,1,s],s]
                Ai              =  RR$res$Co[,2:(m*p[i,2,s]+1),s]
             #Ai                =  t(Ai)
            dim(Ai)             =  c(m,m,p[i,2,s])
                Ao[,,1:p[i,2,s],i,s]            =  Ai
      
      if (type=="none")    Co[,,i,s]                    = RR$res$Co[,1,s]
        if (type=="const")   Co[,,i,s]                  = RR$res$Co[,1,s]
        if (type=="exog1")   Co[,1:(1+p[i,3,s]),i,s]    = RR$res$Co[,c(1,(m*p[i,2,s]+2):(m*p[i,2,s]+p[i,3,s]+1)),s]
        if (type=="exog0")   Co[,1:(1+p[i,3,s]),i,s]    = RR$res$Co[,c(1,(m*p[i,2,s]+2):(m*p[i,2,s]+p[i,3,s]+1)),s]
        }
      Sigmao[(m*(i-1)+1):(i*m),(m*(i-1)+1):(i*m),]      = RR$res$Sigmao
        resid[,,i,]                                     = RR$res$resid
      GDC = Co; 
      dim(GDC) = c(m*n,k,S)
        ##PS[,i,]                                       = RR$PS
        #resid[,(m*(i-1)+1):(i*m)] = RR$res$resid
                   
}

      for (s in 1:S) {
                for (i in 1:n) {
                        for (j in 1:n)  Sigmao[(m*(i-1)+1):(i*m),(m*(j-1)+1):(j*m),s] = t(resid[,,i,s])%*%(resid[,,j,s])/(dim(Z)[1]-dim(Z)[2])

                }
        }
 
      for (s in 1:S)  { 
          for (ss in 1:S)  {
               for (i in 1:n)  {
                   for (j in 1:n) SigmaS[ ((s-1)*n*m+(i-1)*m+1):((s-1)*n*m+i*m), ((ss-1)*n*m+(j-1)*m+1):((ss-1)*n*m+j*m)] = t(resid[,,i,s])%*%(resid[,,j,ss])/sum( (!resid[,1,i,s]==0)*(!resid[,1,j,ss]==0))         
                   }
            }
        }
        
        ### put the estimates into the GVARData0 object 
        ### c("Y","Uo","B","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
        for (s in 1:S) {
                BBo = Bo[,,,,s]
                AAo = Ao[,,,,s]
                dim(BBo) = c(m,m,Pmax,n)
                dim(AAo) = c(m,m,Pmax,n)
                Go[,,,s] = BoAoW2G(BBo,AAo,W,m,n,Pmax) 
        }
      
      
      dim(resid) = c(T,m*n,S)
        res$Go  = Go
        res$C       = CC
        res$Sigmao  = Sigmao
        res$r_npo   = NA
        res$Ao      = Ao
        res$Bo      = Bo
        res$Co      = Co
        res$GDC     = GDC
        res$SigmaS  = SigmaS      
        res$resid   = resid
        res$VAR_domestic = VAR_domestic
        return(res)
}



#' Summary of MRGVAR estimation results 
#'
#' This function sumerizes the estimation results of GVAR in its country VAR model
#'
#' @param  res  :a list of the output of MRGVARest0
#' @examples 
#'
#' res_dd <- MRGVARData2(m=2,n=3,p=1,T=300,S=2,SESVI=c(1,1,1))
#' res_est = MRGVARest0(res_dd)
#'
#'  max(abs(res_dd$Y))
#'
#' summary_MRGVAR(res_est)
#' res_dd$Bo
#' res_est$Bo 
#' @export
summary_MRGVAR=function(res) {
  n = length(res$VAR_domestic) 
  S = res$S
  for ( i in 1:n ) {
    for ( s in 1:S ) {
        print(summary(res$VAR_domestic[[i]]$vars_obj[[s]])) 
    }
  }
}


#' @export
BoAoWs2Gs = function(Bo,Ao,W,m,n,p,state) {
  Bo_s=Bo[,,,,1]*0;dim(Bo_s) = c(m,m,p,n)
  Ao_s=Ao[,,,,1]*0;dim(Ao_s) = c(m,m,p,n)
  for (i in 1:n ) {
     Bo_s[,,,i] = Bo[,,,i,state[i]]
     Ao_s[,,,i] = Ao[,,,i,state[i]]
  }
  Go_s = BoAoW2G(Bo_s,Ao_s,W,m,n,p) 
  return(Go_s)
}


#' @export
SigmaNPD = function(res,StateT) {
      m = res$m
      p = res$p
      if (length(p)>1 ) pp = max(p[,1:2,]) else pp = p
      n = length(StateT);
      dim(StateT) = c(1,n);
      resid = res$resid[,,1]*0
      for (i in 1:n)  {
          resid[,(1+(i-1)*m):(m*i)]=res$resid[,(1+(i-1)*m):(m*i),StateT[i]]
      }
      
      sigmaT = diag(n*m)
	for (i in 1:n)        {
      	for (j in 1:n)  { 
                NN = sum((abs(resid[,((i-1)*m+1):(i*m)])>0)*(abs(resid[,((j-1)*m+1):(j*m)])>0))/m-pp*m
                if (NN>5)  sigmaT[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = t(resid[,((i-1)*m+1):(i*m)])%*%resid[,((j-1)*m+1):(j*m)]/NN
                   else    sigmaT[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = t(res$resid[,((i-1)*m+1):(i*m),1]+res$resid[,((i-1)*m+1):(i*m),2])%*%(res$resid[,((j-1)*m+1):(j*m),1]+res$resid[,((j-1)*m+1):(j*m),2])/dim(resid)[1]
      }
      }
      sigmanpd = as.matrix(Matrix::nearPD(sigmaT,conv.tol = 1e-10)[[1]])
      return(sigmanpd) 
}

#SigmaNPD(res_est,c(1,2,2))


#' Regime specific impulse response functions of MRGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRGVAR(n,p,S). 
#' Using the estimated G[,,,s] and Sigma[,,s] matrices of the MRGVAR, this function calculated the regime speicfic impulse response functions.
#' @param res a list of estimated MRGVAR as output of MRGVARest 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param state an n vector specifying the speciic state for each country. 
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol Cholezky decomposition, Chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse. 
#' @return a list constaining the impulse response functions and the accumulated impulse response function, and the boostrap parameters as well.  
#' @examples 
#'
#' ## case of n = 5, m =3, S = 2 #m: number of variables, n: number of countries
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' res_dd <- MRGVARData(m=3,n=5,p=p,T=500,S=2,SESVI=((1:5)*3-2))
#' max(abs(res_d$Y)) #to check stationarity 
#' res_e = MRGVARest(res=res_dd)
#' summary_MRGVAR(res_e)
#' IRFF2 = irf_MRGVAR(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#' 
#' IRF  = irf_MRGVAR_CB1(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf="gen1",runs=20,conf=c(0.05,0.95))
#' #### under the assumption of stationarity the IRF can be calculated.
#' x11()
#' par(mfrow=c(4,4))
#' plott(IRF[[1]],1,1)
#' plott(IRF[[1]],1,2)
#' plott(IRF[[1]],1,3)
#' plott(IRF[[1]],1,4)
#' plott(IRF[[1]],1,5)
#' plott(IRF[[1]],1,6)
#' plott(IRF[[1]],1,7)
#' plott(IRF[[1]],1,8)
#' plott(IRF[[1]],1,9)
#' plott(IRF[[1]],1,10)
#' plott(IRF[[1]],1,11)
#' plott(IRF[[1]],1,12)
#' plott(IRF[[1]],1,13)
#' plott(IRF[[1]],1,14)
#' plott(IRF[[1]],1,15)
#' plott(IRF[[1]],2,1)
#'
#' @export
#'
irf_MRGVAR = function(res=res_e,state=State,nstep,comb=comb_all,irf = c("gen", "chol", "chol1","gen1","genN1", "comb1","smat"),G=NA,smat=NA,sigmaNPDS=NA) {
      if (missing(G)) 		G = NA
      if (missing(sigmaNPDS))	sigmaNPDS=NA
      if (missing(smat))      smat = NA
      neq 	= dim(res$Go)[1]
	nvar	= dim(res$Go)[2]
      m 	= res$m
	n 	= res$n
      p 	= res$p
      if (length(p)>1) pp = max(p[,1:2,]) else pp = p
	Bo    = res$Bo
	Ao 	= res$Ao
      W	= res$W
      B     = BoAoWs2Gs(Bo,Ao,W,m,n,pp,state)
      #if (anyNA(sigmaNPDS))	sigma = SigmaNPD(res,state)  else  sigma = sigmaNPDS
      if (anyNA(sigmaNPDS))	sigma = SigmaNPDSelectR(res,state)   else  sigma = sigmaNPDS

      response <- array(0,dim=c(neq,nvar,nstep));
      response <- irf_B_sigma(B,sigma,nstep,comb,irf=irf,G=G,smat=smat)
	return(response)
} 

#' Regime specific impulse response functions of MRGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRGVAR(n,p,S). 
#' Using the estimated G[,,,s] and Sigma[,,s] matrices of the estimated MRGVAR, this function calculated the regime speicfic impulse response functions.
#' @param res a list of estimated MRGVAR as output of MRGVARest0 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param state an n vector specifying the speciic state for each country. 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a list of bootstrap result. The first component contains the impulse response functions with confidence bands. It is an (mn,mn,nstep,3) array where the IRF colunms respesent the impulse and rows represent the responses. 
#' @examples 
#'
#' 
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' res_dd <- MRGVARData(m=3,n=5,p=p,T=100,S=2,SESVI=rep(1,5))
#' max(abs(res_dd$Y))
#' ##Ao=res_dd$Ao
#' ##Bo=res_dd$Bo
#' res_d <- MRGVARData(m=3,n=5,p=p,T=500,S=2,SESVI=rep(1,5),Ao=Ao,Bo=Bo)
#' res_d$type
#' max(abs(res_d$Y)) #to check stationarity 
#' res_e = MRGVARest(res=res_d)
#' summary_MRGVAR(res_e)
#' IRFF2 = irf_MRGVAR(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#' IRF  = irf_MRGVAR_CB(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf="gen1",runs=20,conf=c(0.05,0.95))
#' ### under the assumption of stationarity the IRF can be calculated.
#' x11()
#' par(mfrow=c(2,1))
#' plott(IRF[[1]],1,1)
#' plott(IRF[[1]],1,2)
#'
#' @export
#'
irf_MRGVAR_CB = function(res,nstep,comb,state=c(2,1),irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95),NT=1) {
	m = res$m
	n = res$n
      p = res$p
      if (length(p)>1) pp = max(p[,1:2,]) else pp = p
      T = dim(res$Y)[1]*NT
      W = res$W
      S = res$S
      SV = res$SV
      SESVI = res$SESVI

      TH = res$TH

      Ao= res$Ao;
      Bo= res$Bo;
      Co= res$Co;
      Go= res$Go;
      GoColect = array(0,c(dim(Go),runs))
  	BoColect = array(0,c(dim(Bo),runs))
  	AoColect = array(0,c(dim(Ao),runs))
      UoColect = array(0,c(T,m*n,S,runs))
      YColect  = array(0,c(T,m*n,runs))

      type=res$type
      X   = res$X
      mu  = res$mu
      
	B     = BoAoWs2Gs(Bo,Ao,W,m,n,pp,state)

      neq 	= dim(B)[1]
	nvar	= dim(B)[2] 

	sigma = res$Sigmao
      response <- array(0,dim=c(neq,nvar,nstep,length(conf)+1))
      sigmaNPDS = SigmaNPD(res,state)

      response[,,,1] <- irf_MRGVAR(res,nstep,comb,state,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs))

      Uo_run = array(0,c(T,n*m,S))

      for (i in 1:runs) {
            for (s in 1:S)    Uo_run[,,s] = rnormSIGMA(T,res$Sigmao[,,s])
            if (length(p)>1)	{  
               		#res_run = MRGVARData(m,n,p,T,S,W,SESVI,TH,res$Go,Ao=NA,Bo=NA,Sigmao=NA,Uo=Uo_run,SV,type,Co=NA,X)
               		res_run  = MRGVARData(m,n,p,T,S,W,SESVI,TH,res$Go,Ao=NA,Bo=NA,Sigmao=NA,Uo=Uo_run,SV,type,res$Co,X)

			  	res_e    = MRGVARest(res_run)
		}

            if (length(p)==1) {
                        KK = 0
				repeat {
                           KK = KK +1  
				   res_run = MRGVARData2(m,n,p,T,S,W,SESVI,TH,res$Go,Ao=NA,Bo=NA,Sigmao=NA,Uo=Uo_run,SV,type,res$Co,X)
                           #if (((min(colSums(res_run$Y[,seq(2,65,3)]<1.6))> 20)&(max(colSums(res_run$Y[,seq(2,65,3)]<1.6))> 20))|(KK>50)) {
                           #	break
				   #}     
				}				
                        print(KK)
            		res_e   = MRGVARest0(res_run)
		}
	     #responseR[,,,i]  <- irf_MRGVAR(res_e,nstep,comb,state,irf,sigmaNPDS)
            responseR[,,,i]  <- irf_MRGVAR(res_e,nstep,comb,state,irf,SigmaNPD(res_e,state))
            GoColect[,,,,i]  <- res_e$Go
            BoColect[,,,,,i] <- res_e$Bo
            AoColect[,,,,,i] <- res_e$Ao
            UoColect[,,,i]   <- res_run$Uo 
 		YColect[,,i]     <- res_run$Y                   
	}
      
      for (tt in 1:(nstep) )      {
	for (i in 1:neq)            {
		for (j in 1:nvar)     {
			response[i,j,tt,-1] = quantile(responseR[i,j,tt,], conf)
			
	} 
	}
	}
	return(list(response,GoColect,BoColect,AoColect,UoColect,YColect))
} 



#' Generalized impulse response functions of MRGVAR(n,p,S) with regime migrations 
#'
#' This function calculates the generalized impulse response functions of an estimated MRGVAR(n,p,S) for given a shock vector.
#'                   GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK)) 
#'
#' @param  res   : a list containing the components as the output of MRGVARest0.
#' @param  shock : an mn vector containing the shocks as impulse.
#' @param  R     : the number runs to integrate out the random effects in order to obtain the means (see equation above).
#' @param  nstep : the length of the responses
#' @param  state : the state from which we start the shocks
#' @param  resid_method : resid_method = c("resid", "parametric"), It generate the random residuals from residuals bootstrap or parametric boostrap. 
#' @return an mn x mn x nstep+1 matrix of impulse response functions. The rows represent response the columns represent impulses. 
#' @examples 
#' res_dd <- MRGVARData2(m=2,n=3,p=1,T=300,S=2,SESVI=c(1,1,1))
#' res_est = MRGVARest0(res_dd)
#' max(abs(res_dd$Y)) to check stationarity 
#' summary_MRGVAR(res_est)
#' IRFF2 = irf_MRGVAR(res_est,nstep=10,comb=NA,state=c(2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#' IRF  = irf_MRGVAR_CB(res=res_est,nstep=10,comb=NA,state=c(1,1,1),irf="gen",runs=200,conf=c(0.05,0.95))
#' plott(IRF,1,1)
#'
#' GIRF = girf_MRGVAR_RM(res=res_est,shock=c(1,1,1,1,1,1),R=100,nstep=10,state=c(1,1,1),resid_method='parametric')
#' GIRF_CB = girf_RM_CB(res=res_est,shock=c(1,1,1,1,1,1),R=100,nstep=10,state=c(1,1,1),resid_method='parametric',conf=c(0.05,0.95),N=100)
#'
#' @export
girf_MRGVAR_RM = function(res,shock,R,nstep,state,resid_method) {
####  this function generate the impulse response funciton of MRVAR with regime migration
####
####                  GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK)) 
####   
####  
#### RES is the output of MRVAREST1 
      n 		= res$n
      m 		= res$m
      p 		= res$p
	d		= res$d
      S 		= res$S
      SESVI		= res$SESVI
      TH		= res$TH
      Go          = res$Go
      Bo		= res$Bo
      Ao          = res$Ao
      Co		= res$Co
      Sigmao	= res$Sigmao
      W           = res$W
      type		= res$type
      SV          = res$SV
      P           = max(d,p)
      X 		= res$X; if (!anyNA(X)) X = X[1:(P+nstep+1),];
      Yo		= res$Yo
      Uo          = res$Uo
      TT          = dim(res$Y)[1]
	d           = res$d
	IRF 		= list()
      YR          = list()
      YS          = list()
      residR   <-  list()
      residS   <-  residR
      IRFHLPR  	= (1:(m*n*(P+nstep+1)*R))*0;     dim(IRFHLPR) = c(m*n,P+nstep+1,R) 
      MeanR       = (1:(m*n*(P+nstep+1)))*0;       dim(MeanR)   = c(m*n,P+nstep+1)
	IRFHLPS  	= (1:(m*n*m*n*(P+nstep+1)*R))*0; dim(IRFHLPS) = c(m*n,m*n,P+nstep+1,R)
      MeanS       = (1:(m*n*m*n*(P+nstep+1)))*0;   dim(MeanS)   = c(m*n,m*n,P+nstep+1)

      #R 		= 400
      #shock       = (1:n)/(1:n)
      DIMresid    = dim(Uo)
      if (length(DIMresid) == 3)  residI = Uo[1:(P+1+nstep),,]*0 
      if (length(DIMresid) == 2)  residI = Uo[1:(P+1+nstep),]*0
      shock_mat = matrix(0,m*n,m*n)
      Sigma_s = SigmaNPD(res,state)

      for (i in 1:(m*n)) {   
          shock_mat[i,i]  = shock[i]; 
          shock_mat[-i,i] = rnormSIGMA_cond(T=1,sigma=Sigma_s,I=i, v=shock[i],switch=0)
      }
      for (i in 1:R) { 
      if ( resid_method=="resid" ) {
		if (length(DIMresid) == 2)  { 
			residI  = Uo[NoutofT(P+nstep+1,DIMresid[1]),]; 
            }
            if (length(DIMresid) == 3)  {
            	for (j in 1:DIMresid[3]) {
				residI[,,j]  = Uo[NoutofT(P+nstep+1,DIMresid[1]),,j];                        
                  }
            } 
      }
      if ( resid_method=="parametric" ) {
		if (length(DIMresid) == 2)  {
			residI   = rnormSIGMA(P+nstep+1,Sigmao) 
		}
            if (length(DIMresid) == 3)  {
			for (j in 1:DIMresid[3]) {
                        residI[,,j] = rnormSIGMA(P+nstep+1,Sigmao[,,j]) 
                  }
            }
      }
      residR[[i]] = residI
      residS[[i]] = residI 
	YR[[i]] = MRGVARData2(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residR[[i]],SV=SV,type=type,Co=Co,X=X,Yo=Yo*0,d=d)
	IRFHLPR[,,i] = t(YR[[i]]$Y)

	for (k in 1:DIMresid[2]) {
            	if (length(DIMresid) == 2)  residS[[i]][P+1,]  = t(shock_mat[,k])  
			if (length(DIMresid) == 3)  residS[[i]][P+1,,] = t(as.matrix((1:DIMresid[3])/(1:DIMresid[3]))%*%shock_mat[,k])
                  
            	YS[[i]] = MRGVARData2(m=m,n=n,p=p,T=(P+nstep+1),S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=residS[[i]],SV=SV,type=type,Co=Co,X=X,Yo=Yo*0,d=d)
                  IRFHLPS[,k,,i] = t(YS[[i]]$Y)

    	}
    }  # end of R loop	
    ############## integrating out the disturbances
      for (i in 1:(m*n)) {
          MeanR[i,] = rowMeans(IRFHLPR[i,,],dims=1)
	    for (j in 1:(m*n)) {
                 MeanS[i,j,]      = rowMeans(IRFHLPS[i,j,,],dims=1)
                 MeanS[i,j,]      = MeanS[i,j,] - MeanR[i,] 
      }
      }
   return(MeanS)
}


#' Generalized impulse response functions of MRGVAR(m,p,n,S,W,TH) with regime migrations and confidence intervals
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
#' @return an n x n x nstep+1 x 3 containing of impulse response functions with lower and upper confidence bonds. The rows represent response the columns represent impulses. 
#' @examples 
#' res_dd <- MRGVARData2(m=2,n=3,p=1,T=300,S=2,SESVI=c(1,1,1))
#' res_est = MRGVARest0(res_dd)
#' max(abs(res_dd$Y)) to check stationarity 
#' summary_MRGVAR(res_est)
#' IRFF2 = irf_MRGVAR(res_est,nstep=10,comb=NA,state=c(2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#' IRF  = irf_MRGVAR_CB(res=res_est,nstep=10,comb=NA,state=c(1,1,1),irf="gen",runs=200,conf=c(0.05,0.95))
#' plott(IRF,1,1)
#'
#' GIRF = girf_MRGVAR_RM(res=res_est,shock=c(1,1,1,1,1,1),R=100,nstep=10,state=c(1,1,1),resid_method='parametric')
#' GIRF_CB = girf_RM_CB(res=res_est,shock=c(1,1,1,1,1,1),R=100,nstep=10,state=c(1,1,1),resid_method='parametric',conf=c(0.05,0.95),N=100)
#'
#' @export
girf_RM_CB = function(res,shock,R,nstep,state,resid_method="parametric",conf,N) {
##### this is to generate bootstrap GIRF for MRVAR with regime migrations 
#####
      n 		= res$n
      m 		= res$m
      p 		= res$p
      S 		= res$S
      SESVI		= res$SESVI
      TH		= res$TH
      Go          = res$Go
      Bo		= res$Bo
      Ao          = res$Ao
      Co		= res$Co
      Sigmao	= res$Sigmao
      W           = res$W
      type		= res$type
      SV          = res$SV
      X 		= res$X; 
      Yo		= res$Yo
	d           = res$d
      Uo          = res$Uo
      T           = dim(res$Y)[1]
      P           = max(d,p)
      GIRF 		= (1:(m*n*n*m*(P+nstep+1)*N))*0; dim(GIRF) = c(m*n,m*n,P+nstep+1,N)

      GIRFBd	= (1:(m*n*n*m*(P+nstep+1)*(length(conf)+1)))*0; dim(GIRFBd) = c(m*n,m*n,P+nstep+1,length(conf)+1)

	GIRFBd[,,,1]= girf_MRGVAR_RM(res,shock,R,nstep,state,resid_method)
      
      for (i in 1:N) {
	      res_run = MRGVARData2(m=m,n=n,p=p,T=T,S=S,W=W,SESVI=SESVI,TH=TH,Go=Go,Ao=Ao,Bo=Bo,Sigmao=Sigmao,Uo=NA,SV=SV,type=type,Co=Co,X=X,Yo=Yo,d=d)
      	res_e   = MRGVARest0(res_run)
      	RF3     = girf_MRGVAR_RM(res=res_e,shock,R,nstep,state,resid_method)
		GIRF[,,,i]  = RF3
      }
	

      for (tt in 1:(P+nstep+1) ) {
	for (i in 1:(n*m))           {
		for (j in 1:(n*m))     {
			GIRFBd[i,j,tt,-1] = quantile(GIRF[i,j,tt,], conf)
			
	} 
	}
	}
	return(GIRFBd)	
}



#' Calculation of information criteria for a MRGVAR model  
#' 
#' @param res  : a MRGVAR object obtained from MRGVARData or estimated from MRGVARest.
#' @param I    : index of the country under investigation.
#' @param L_V  : a four components vector containing the maxima of the domestic lag and the foreign lag for each regime, respectively. 
#' @param TH_V : a vector containing possible threshold values. 
#' @return      : A matrix with different lag specifications and values of the model selection criteria.
#' @examples 
#' 
#' ## case of n = 5, m =3, S = 2 #m: number of variables, n: number of countries
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p[,2,] = 1; p[1,1,1] = 1; p[3,1,1] = 1; p[2,1,2] = 1
#'
#' TH = c(1:5)*0; dim(TH) = c(1,5)
#' res_dd <- MRGVARData(m=3,n=5,p=p,TH =TH,T=100,S=2,SESVI=((1:5)*3-2))
#' max(abs(res_dd$Y))		# to make sure it is not explosive
#' ##Ao=res_dd$Ao
#' ##Bo=res_dd$Bo
#' #res_d <- MRGVARData(m=3,n=5,p=p,TH= TH, T=2000,S=2,SESVI=((1:5)*3-2),Ao=Ao,Bo=Bo)
#' res_d$type
#' max(abs(res_d$Y)) #to check stationarity 
#' res_e = MRGVARest(res=res_d)
#' #summary_MRGVAR(res_e)
#' 
#' res_dd$TH
#' ### four numbers for the maxima lag length in counttry I: regime 1: (domesti foreign regime 2: domestic and foreign) 
#' L_v  = c(3,3,3,3)
#' ### a vector eontaining possible threhsold values 
#' TH_v = c(-0.1, -0.05, 0,0.05, 0.1  )
#' 
#' 
#' CC = MRGVAR_Select(res=res_dd,I=1,L_V=L_v,TH_V=TH_v)
#' CCC = CC[[1]]
#'
#' CCC[which.min(CCC[,9]),]
#' CCC[which.min(CCC[,18]),]
#'  
#' 
#' CC = MRGVAR_Select(res=res_dd,I=2,L_V=L_v,TH_V=TH_v)
#' CCC = CC[[1]]
#' 
#' CCC[which.min(CCC[,9]),]
#' CCC[which.min(CCC[,18]),]
#' 
#' 
#' CC = MRGVAR_Select(res=res_dd,I=3,L_V=L_v,TH_V=TH_v)
#' CCC = CC[[1]]
#' 
#' CCC[which.min(CCC[,9]),]
#' CCC[which.min(CCC[,18]),]
#'
#' @export
#'
MRGVAR_Select = function(res,I,L_V,TH_V) {
##    res   object of MRGVARData2
##    I     index of the equation under investigation, I =1 we investigate the first equation, I = 10 we investigate the 10th equation.  
##    L_V   maximum of the lags of regimes (4 4 4 4) Regime 1 domestic lags 4   foreign lags 4   regime 2 lags 4  4 
##    TH_V  vector of threshold values that will be estimated for each chosen lag combinations  
m= res$m
n= res$n
p= res$p
Y= as.matrix(res$Y)
X     = res$X
W = res$W
type  = res$type
TH    = res$TH
SESVI = res$SESVI
Go= res$Go
S     = res$S
Ao= res$Ao
Bo    = res$Bo
Co    = res$Co
GDC   = res$GDC
Sigmao= res$Sigmao*0
d     = res$d
P     = max(p,d) 
Pmax  = max(p[,1:2,])
SigmaS = (1:(n*m*S*n*m*S))*0
dim(SigmaS) = c(n*m*S,n*m*S)
VAR_domestic = list()

foreignLagSame = 1


k1 = p[I,3,1]
k2 = p[I,3,2]
kmax  = max(p[,3,])

T= dim(Y)[1]
FY= Y%*%t(W%x%diag(m))
resid = (1:(T*m*n*S))*0; dim(resid) = c(T,m,n,S)
##PS  = (1:(T*n*S))*0;   dim(PS)    = c(T,n,S)
###   CC is the single country coefficients of determinsitic components i.e. the foreign variables and the deterministic variables
###   It is m x (m*p+k+1) where m*p is the number of foreign variables and k is the number of deterministic variables. with zeros for "none"  

P     = max(L_V,d)
P_candidate = array(0,c(L_V[1],L_V[2],L_V[3],L_V[4]))
Criteria   = matrix(0,prod(L_V),4+length(TH_V))
Criteria_b = matrix(0,prod(L_V),4+length(TH_V))
Criteria_c = matrix(0,prod(L_V),4+length(TH_V))
Criteria_d = matrix(0,prod(L_V),4+length(TH_V))

if (type=="none" |type=="const")  k = 1;
if (type=="exog0"|type=="exog1")  k = dim(X)[2]+1;
CC  =  c(1:(m*(Pmax*m+k)*S))*0
dim(CC) = c(m,m*Pmax+k,S) 

#res_MSVAR = MSVARData(n=m,p=p,T=T,S=S,Co=CC,TM=TM,type="exog",X=Z) 
#res_MSVAR = MSVARData(n=m,p=p,T=T,S=S)

EMCVG       = c(1:n)*0

Pppp = matrix(0,3,2) 

for (h in 1:length(TH_V)) {
    idx = 0
for (l_1 in 1:L_V[1])             {
for (l_2 in 1:L_V[2])       {
for (l_3 in 1:L_V[3])  {
for (l_4 in 1:L_V[4]){
      if ( foreignLagSame == 1) l_4=l_2 
      Pppp[,1] = c(l_1,l_2,k1);
      Pppp[,2] = c(l_3,l_4,k2);

      idx = idx+1;
for (i in I:I) {
      p[i,,]      = Pppp
      Pi          = max(p[i,1:2,])
CCi  = c(1:(m*(Pi*m+k)*S))*0
dim(CCi) = c(m,m*Pi+k,S) 
Z         = matrix(0,T,Pi*m)
      FYp    = embed(FY[,(m*(i-1)+1):(i*m)],(Pi+1))
Z[(Pi+1):T,]= FYp[,(m+1):(Pi*m+m)]
      kz          = dim(Z)[2]
if (type=="none")   {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz,S)); for (s in 1:S) ZZ[,1:kz,s]=Z[,,s]}
if (type=="const")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz,S)); for (s in 1:S) ZZ[,1:kz,s]=Z[,,s]}
if (type=="exog0")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz+kmax,S)); for (s in 1:S) ZZ[,1:(kz+p[i,3,s]),s] = cbind(Z[,,s],X[,1:p[i,3,s],i,s])}
      if (type=="exog1")  {Z = rep(Z,S); dim(Z) = c(T,kz,S); ZZ = array(0,c(T,kz+kmax,S)); for (s in 1:S) ZZ[,1:(kz+p[i,3,s]),s] = cbind(Z[,,s],X[,1:p[i,3,s],i,s])}
if (type=="none"|type=="exog0")  type_holder = "exog0"
if (type=="const"|type=="exog1") type_holder = "exog1"
      XX = (1:(T*dim(ZZ)[2]*S))*0; dim(XX) = c(T,dim(ZZ)[2],S)
      if (!anyNA(X)) for (s in 1:S) XX[,1:(m*p[i,2,s]+p[i,3,s]),s] = ZZ[,c(1:(m*p[i,2,s]),(m*Pi+1):(m*Pi+p[i,3,s])),s]
      if (anyNA(X))  for (s in 1:S) XX[,1:(m*p[i,2,s]),s] = ZZ[,c(1:(m*p[i,2,s])),s]

### put data into a MRVARData object in order to replace the data and order parameters
      pp = t(p[i,,]); pp[,2] = pp[,2]*m; ppp = pp[,1:2]; ppp[,2] = pp[,2]+pp[,3]
      if (SESVI[i]/m > SESVI[i]%/%m)   SESVIi = SESVI[i]- (SESVI[i]%/%m)*m  else SESVIi = m

     if (max(ppp[,1])>=max(ppp[,2])%/%m) {
      res_MRVAR =  MRVARData(n=m,p=ppp,T=T,S=S,Co=NA,SESVI=SESVIi,type=type_holder,X=XX)
      res_MRVAR$Y       =  Y[,(m*(i-1)+1):(i*m)]
      res_MRVAR$X       =  XX
      res_MRVAR$SV      =  res_MRVAR$Y[,SESVIi] 
      res_MRVAR$TH      =  TH_V[h]
      
      }  else  {
      
      startT = max(ppp[,2])%/%m - max(ppp[,1])+1 
      res_MRVAR =  MRVARData(n=m,p=ppp,T=(T-startT+1),S=S,Co=NA,SESVI=SESVIi,type=type_holder,X=XX[startT:T,,])
      res_MRVAR$Y       =  Y[startT:T,(m*(i-1)+1):(i*m)]
      res_MRVAR$X       =  XX[startT:T,,]
      res_MRVAR$SV      =  res_MRVAR$Y[,SESVIi] 
      res_MRVAR$TH      =  TH_V[h]
      }
      RR=  MRVARest(res=res_MRVAR)
      Criteria[idx,c(c(1:4),4+h)]     = c(l_1,l_2,l_3,l_4,RR$LH_AIC)
      Criteria_b[idx,c(c(1:4),4+h)]   = c(l_1,l_2,l_3,l_4,RR$LH_BIC)
      Criteria_c[idx,c(c(1:4),4+h)]   = c(l_1,l_2,l_3,l_4,RR$LH_P)
      Criteria_d[idx,c(c(1:4),4+h)]   = c(l_1,l_2,l_3,l_4,RR$LH_N)
    
      Criterion = cbind(Criteria,Criteria_b,Criteria_c,Criteria_d)
     
}
}
}
}
}
}
return(list(Criterion,RR))
}



#' Regime specific impulse response functions of MRGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRGVAR(n,p,S). 
#' Using the estimated G[,,,s] and Sigma[,,s] matrices of the estimated MRGVAR, this function calculated the regime speicfic impulse response functions. A feature of this function is that 
#' it will gnerate sufficient boostrap data to calculate the boostrap confidence interval.  
#' 
#' @param res a list of estimated MRGVAR as output of MRGVARest0 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param state an n vector specifying the speciic state for each country. 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a list of bootstrap result. The first component contains the impulse response functions with confidence bands. It is an (mn,mn,nstep,3) array where the IRF colunms respesent the impulse and rows represent the responses. 
#' @examples 
#' ## case of n = 5, m =3, S = 2 #m: number of variables, n: number of countries
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' res_dd <- MRGVARData(m=3,n=5,p=p,T=500,S=2,SESVI=((1:5)*3-2))
#' max(abs(res_d$Y)) #to check stationarity 
#' res_e = MRGVARest(res=res_dd)
#' summary_MRGVAR(res_e)
#' IRFF2 = irf_MRGVAR(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#' 
#' IRF  = irf_MRGVAR_CB1(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf="gen1",runs=20,conf=c(0.05,0.95))
#' #### under the assumption of stationarity the IRF can be calculated.
#' x11()
#' par(mfrow=c(4,4))
#' plott(IRF[[1]],1,1)
#' plott(IRF[[1]],1,2)
#' plott(IRF[[1]],1,3)
#' plott(IRF[[1]],1,4)
#' plott(IRF[[1]],1,5)
#' plott(IRF[[1]],1,6)
#' plott(IRF[[1]],1,7)
#' plott(IRF[[1]],1,8)
#' plott(IRF[[1]],1,9)
#' plott(IRF[[1]],1,10)
#' plott(IRF[[1]],1,11)
#' plott(IRF[[1]],1,12)
#' plott(IRF[[1]],1,13)
#' plott(IRF[[1]],1,14)
#' plott(IRF[[1]],1,15)
#' plott(IRF[[1]],2,1)
#'
#'
#' @export
irf_MRGVAR_CB1 = function (res, nstep, comb, state = c(2, 1), irf = c("gen", "chol", "chol1", "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95), NT = 1) 
{
    m = res$m
    n = res$n
    p = res$p
    if (length(p) > 1) 
        pp = max(p[, 1:2, ])
    else pp = p
    T = dim(res$Y)[1] * NT
    W = res$W
    S = res$S
    SV = res$SV
    SESVI = res$SESVI
    TH = res$TH
    Ao = res$Ao
    Bo = res$Bo
    Co = res$Co
    Go = res$Go
    GoColect = array(0, c(dim(Go), runs))
    BoColect = array(0, c(dim(Bo), runs))
    AoColect = array(0, c(dim(Ao), runs))
    UoColect = array(0, c(T, m * n, S, runs))
    YColect = array(0, c(T, m * n, runs))
    type = res$type
    X = res$X
    mu = res$mu
    B = BoAoWs2Gs(Bo, Ao, W, m, n, pp, state)
    neq = dim(B)[1]
    nvar = dim(B)[2]
    sigma = res$Sigmao
    response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
    accresponse <- response
    sigmaNPDS = SigmaNPD(res, state)
    response[, , , 1] <- irf_MRGVAR(res=res, nstep=nstep, comb=comb, state=state, irf=irf)
    accresponse[, , , 1] <- response[, , , 1]
    for (tt in 2:nstep) {  accresponse[, ,tt, 1] = accresponse[, ,tt-1, 1] + response[, ,tt, 1] }

    responseR <- array(0, dim = c(neq, nvar, nstep, runs))
    accresponseR <- responseR
    Uo_run = array(0, c(T, n * m, S))
    for (i in 1:runs) {
        for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, res$Sigmao[,, s])
        if (length(p) > 1) {
            res_run = MRGVARDataR(res)
            res_e   = MRGVARest(res_run)
        }
        if (length(p) == 1) {
            KK = 0
            repeat {
                KK = KK + 1
                res_run = MRGVARData2(m, n, p, T, S, W, SESVI, TH, res$Go, Ao = NA, Bo = NA, Sigmao = NA, Uo = Uo_run, SV, type, res$Co, X)
            }
            print(KK)
            res_e = MRGVARest0(res_run)
        }
        
        
        responseR[, , , i] <- irf_MRGVAR(res=res_e, nstep=nstep, comb=comb, state=state, irf=irf, sigmaNPDS=sigmaNPDS)
        accresponseR[, , , i] <-  responseR[, , , i]
        for (tt in 2:nstep) {  accresponseR[, ,tt, i] = accresponseR[, ,tt-1, i] + responseR[, ,tt, i] }
        GoColect[, , , , i] <- res_e$Go
        BoColect[, , , , , i] <- res_e$Bo
        AoColect[, , , , , i] <- res_e$Ao
        #UoColect[, , , i] <- res_run$Uo
        #YColect[, , i] <- res_run$Y
    }
    responseR[, , , 1]    = response[, , , 1]
    accresponseR[, , , 1] = accresponse[, , , 1]

    for (tt in 1:(nstep)) {
        for (i in 1:neq) {
            for (j in 1:nvar) {
                response[i, j, tt, -1]    = quantile(responseR[i, j, tt, ], conf)
 		    accresponse[i, j, tt, -1] = quantile(accresponseR[i, j, tt, ], conf)
            }
        }
    }
    return(list(response, accresponse, GoColect, BoColect, AoColect, UoColect, YColect,responseR))
}

#' Data generating process of MRGVARDataR(res) 
#'
#' This function will generate data from a MRGVAR object. It will generate enough data for estimation purpose. 
#'
#' @param res     : an output of MRGVARest
#' @return	: a MRGVAR object. 
#'
#' @export
MRGVARDataR=function(res) {
### res_e : an estimated MRGVAR model that is an output of MRGVARest 
### T     : number of observations
### Remarks: MRGVARDataR is used in bootstrapin to generate suffiicent observtions such that in the regime of rare occurance also contains
###          sufficient observations for estimation purposed
###
###  MINH:   minimum number of observations of a nation in the rare occurance regime.  

    m = res$m
    n = res$n
    p = res$p
    if (length(p) > 1) 
        pp = max(p[, 1:2, ])
    else pp = p
    T = dim(res$Y)[1]
    W = res$W
    S = res$S
    SV = res$SV
    SESVI = res$SESVI
    TH = res$TH
    Ao = res$Ao
    Bo = res$Bo
    Co = res$Co
    Go = res$Go
    type = res$type
    X = res$X
    mu = res$mu
    MINH = T
    repeat {
	 Uo_run = array(0, c(T, n * m, S))
       for (s in 1:S) Uo_run[, , s] = rnormSIGMA(T, res$Sigmao[,,s])       #### zhe ge you wenti
    	 res_run = MRGVARData(m, n, p, T, S, W, SESVI, TH, res$Go, Ao = NA, Bo = NA, Sigmao = NA, Uo = Uo_run, SV, type, Co, X)
       for (i in 1:n) { MINH = min(MINH,sum(res_run$Y[,SESVI[i]]>TH[1,i]),sum(res_run$Y[,SESVI[i]]<TH[1,i]))}
       if (MINH > (max(p)*(m*2)+2*m)) break
      else {T=T+T; print(c(MINH,T)); MINH=T}   
    }
  return(res_run)
}



#' @export
ACCIRFconfR = function(IRF) {
      ACCirf = IRF
        dm = dim(IRF)
      for (t in 2:dm[3])          {
                 ACCirf[,,t,] =    ACCirf[,,t-1,] + IRF[,,t,]  
      }
      return(ACCirf)      
}


#' @export
SigmaNPDSelectR = function(res,StateT) {
      n = length(StateT);
      dim(StateT) = c(1,n);
      Ranking = 1:n;
      SigmaS = res$SigmaS
      N = dim(SigmaS)[1]/(2*n)
      sigmaT = diag(n*N)
        for (i in 1:n)        {
        for (j in 1:n)  { 
            sigmaT[((i-1)*N+1):(i*N),((j-1)*N+1):(j*N)] = SigmaS[((StateT[1,i]-1)*n*N+(i-1)*N+1):((StateT[1,i]-1)*n*N+i*N),((StateT[1,j]-1)*n*N+(j-1)*N+1):((StateT[1,j]-1)*n*N+j*N)]   
      }
      }
      if (anyNA(sigmaT)) { sigmanpd = SigmaRKSelect(Model,StateT,Ranking)[[1]] }
      else {
          sigmanpd = as.matrix(Matrix::nearPD(sigmaT,conv.tol = 1e-10)[[1]])
      }
      return(sigmanpd) 
}





#' @export
plott = function (IRF_CB, i, j) 
{
    ylim = c(min(IRF_CB[i, j, , ]), max(IRF_CB[i, j, , ]))
    ylim = c(-10,10)
    plot(IRF_CB[i, j, , 1], type = "l", ylim = ylim)
    lines(IRF_CB[i, j, , 2], col = "red")
    lines(IRF_CB[i, j, , 3], col = "red")
    lines(IRF_CB[i, j, , 3] * 0, col = "black")
}
#' @export
plot35 = function(ACCIRF_band=IRFC[[2]],i,r,ylimm1=c(-1,1),ylimm2=c(-1,1))  {
par(mfrow=c(3,5))
 plot(ACCIRF_band[i[1],r[1],,1],type='l',ylim=ylimm1)
lines(ACCIRF_band[i[1],r[1],,2], col='red')
lines(ACCIRF_band[i[1],r[1],,3], col='blue')
 plot(ACCIRF_band[i[2],r[2],,1], ,type='l',ylim=ylimm1)
lines(ACCIRF_band[i[2],r[2],,2], col='red')
lines(ACCIRF_band[i[2],r[2],,3], col='blue')
 plot(ACCIRF_band[i[3],r[3],,1], type='l',ylim=ylimm1)
lines(ACCIRF_band[i[3],r[3],,2], col='red')
lines(ACCIRF_band[i[3],r[3],,3], col='blue')
 plot(ACCIRF_band[i[4],r[4],,1],type='l',ylim=ylimm1)
lines(ACCIRF_band[i[4],r[4],,2], col='red')
lines(ACCIRF_band[i[4],r[4],,3], col='blue')
 plot(ACCIRF_band[i[5],r[5],,1],type='l',ylim=ylimm1)
lines(ACCIRF_band[i[5],r[5],,2], col='red')
lines(ACCIRF_band[i[5],r[5],,3], col='blue')
 plot(ACCIRF_band[i[6],r[6],,1], ,type='l',ylim=ylimm1)
lines(ACCIRF_band[i[6],r[6],,2], col='red')
lines(ACCIRF_band[i[6],r[6],,3], col='blue')
 plot(ACCIRF_band[i[7],r[7],,1],type='l',ylim=ylimm2)
lines(ACCIRF_band[i[7],r[7],,2], col='red')
lines(ACCIRF_band[i[7],r[7],,3], col='blue')
 plot(ACCIRF_band[i[8],r[8],,1], ,type='l',ylim=ylimm2)
lines(ACCIRF_band[i[8],r[8],,2], col='red')
lines(ACCIRF_band[i[8],r[8],,3], col='blue')
 plot(ACCIRF_band[i[9],r[9],,1], type='l',ylim=ylimm2)
lines(ACCIRF_band[i[9],r[9],,2], col='red')
lines(ACCIRF_band[i[9],r[9],,3], col='blue')
 plot(ACCIRF_band[i[10],r[10],,1],type='l',ylim=ylimm2)
lines(ACCIRF_band[i[10],r[10],,2], col='red')
lines(ACCIRF_band[i[10],r[10],,3], col='blue')
 plot(ACCIRF_band[i[11],r[11],,1],type='l',ylim=ylimm2)
lines(ACCIRF_band[i[11],r[11],,2], col='red')
lines(ACCIRF_band[i[11],r[11],,3], col='blue')
 plot(ACCIRF_band[i[12],r[12],,1], ,type='l',ylim=ylimm2)
lines(ACCIRF_band[i[12],r[12],,2], col='red')
lines(ACCIRF_band[i[12],r[12],,3], col='blue')
 plot(ACCIRF_band[i[13],r[13],,1],type='l',ylim=ylimm2)
lines(ACCIRF_band[i[13],r[13],,2], col='red')
lines(ACCIRF_band[i[13],r[13],,3], col='blue')
 plot(ACCIRF_band[i[14],r[14],,1],type='l',ylim=ylimm2)
lines(ACCIRF_band[i[14],r[14],,2], col='red')
lines(ACCIRF_band[i[14],r[14],,3], col='blue')
 plot(ACCIRF_band[i[15],r[15],,1], ,type='l',ylim=ylimm2)
lines(ACCIRF_band[i[15],r[15],,2], col='red')
lines(ACCIRF_band[i[15],r[15],,3], col='blue')
}

#' @export
plotttt35 = function(ACCIRF,ACIRF,i,ylimm1,k=0)  {
par(mfrow=c(3,5))
 plot(ACCIRF[1+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "USA")
 lines(ACIRF[1+k,i,],col='red')
 plot(ACCIRF[3+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "JPN")
 lines(ACIRF[3+k,i,],col='red')
 plot(ACCIRF[5+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "DEU")
 lines(ACIRF[5+k,i,],col='red')
 plot(ACCIRF[7+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "FRA")
 lines(ACIRF[7+k,i,],col='red')
 plot(ACCIRF[9+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "GBR")
 lines(ACIRF[9+k,i,],col='red')
 plot(ACCIRF[11+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "ITA")
 lines(ACIRF[11+k,i,],col='red')
 plot(ACCIRF[13+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "ESP")
 lines(ACIRF[13+k,i,],col='red')
 plot(ACCIRF[15+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "CAN")
 lines(ACIRF[15+k,i,],col='red')
 plot(ACCIRF[17+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "DNK")
 lines(ACIRF[17+k,i,],col='red')
 plot(ACCIRF[19+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "FIN")
 lines(ACIRF[19+k,i,],col='red')
 plot(ACCIRF[21+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "SWE")
 lines(ACIRF[21+k,i,],col='red')
 plot(ACCIRF[23+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "NOR")
 lines(ACIRF[23+k,i,],col='red')
 plot(ACCIRF[25+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "AUT")
 lines(ACIRF[25+k,i,],col='red')
 plot(ACCIRF[27+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "BEL")
 lines(ACIRF[27+k,i,],col='red')
 plot(ACCIRF[29+k,i,],xlab="",ylab="",  type='l',ylim=ylimm1);title(main = "NLD")
 lines(ACIRF[29+k,i,],col='red')
 
}






