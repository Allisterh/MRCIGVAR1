#' Data generating process of MRCIGVAR(m,n,p,S) 
#'
#' This function will generate data from an multi-regime cointegrated Global VAR(p) process and return a list containing data and parameters used in the MRCIGVAR(m,n,p,S) process.
#'
#' @param m     : number of variables in each a country/unit
#' @param n     : number of countries/units
#' @param S     : number of regimes
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
#' @param r     : n vector containing the number of unit root process in each country
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
#'
#' @examples 
#' ## case of n = 2, m = 2, S = 2
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:2,,]/2; p[1,1,2]=2; p[2,2,2] = 2
#' TH = c(1:2)*0; dim(TH) = c(1,2)
#' res_dd <- MRGVARData3(m=2,n=2,p=p,TH=TH,T=2000,S=2,SESVI=c(1,3))    ## m: number of variables, n: number of countries 
#' #cbind( as.numeric(res_dd$Y[,1]>TH[1,1]), as.numeric(res_dd$Y[,3]>TH[1,2]))[3:8,] - (res_dd$St[4:9,]-1)
#' max(res_dd$Y)
#' res_dd$TH  
#' ### estimation of the MRGBAR model
#' res_ee = MRGVARest1(res=res_dd) 
#' summary_MRGVAR(res_ee)
#'
#'
#' @export
MRCIGVARData0=function(m,n,p,T,S,W=NA,SESVI=NA,TH=NA,Go=NA,Ao=NA,Bo=NA,Sigmao=NA,Uo=NA,SV=NA,type=NA,Co=NA,X=NA,Yo=NA,d=NA,r=NA) {
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
  if (missing(r))       {r = NA}

  if (anyNA(d))  { d = 1 }
  if (anyNA(r))  { r = rep(1,n) }

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
      	GVARD = GVARData(m,n,p[,,s],T) 
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
      CIVARB =  (1:(2*m*2*m*Pmax*S*n))*0; dim(CIVARB)=c(2*m,2*m,Pmax,n,S)
  	Go 	  = (1:(m*n*m*n*Pmax*S))*0
      dim(Go) = c(m*n,m*n,Pmax,S)
	Bo 	  = (1:(m*m*Pmax*n*S))*0
	dim(Bo) = c(m,m,Pmax,n,S)
	Ao      = Bo*0

            #CIGVARD     = GVARData1(m,n,p[,,s],T,W=W)    #(m,n,p,T,W,r_npo,Ao,Go,Uo,Sigmao)  !!!!W must be the same for all regimes
      	#Go[,,1:max(p[,1:2,s]),s]  = GVARD$G
            #check[s]  = max(abs(GVARD$Y))
		#Ao[,,1:max(p[,1:2,s]),,s] = GVARD$Ao 
		#Bo[,,1:max(p[,1:2,s]),,s] = GVARD$Bo
            for (i in 1:n) {
                  pi = matrix(0,2,2)
                  pi[,1] = c(max(p[,,1]),max(p[,,2]));
                  sigma  = cbind(diag(2*m),diag(2*m));dim(sigma) = c(2*m,2*m,2) 
			CIVARD = MRCIVARData(n=2*m,p=pi,T =100,S=S,SESVI=SESVI[i]%%m,TH=TH[i],Sigmao=sigma,type="none",r=r[i])
			#plot(ts(CIVARD$Y))
                  CIVARB[,,,i,] = CIVARD$Bo
			#MRCIVARData(n=4,p=pi,T =100,S=2,SESVI=1,TH=0,Sigmao=sigma,type="none",r=1)
     			Bo[,,1:max(pi[,1]),i,] = CIVARD$Bo[1:m,1:m,,] 
		      Ao[,,1:max(pi[,1]),i,] = CIVARD$Bo[1:m,(m+1):(2*m),,]
            }
           

  	for (s in 1:S ) Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,Pmax)
 	## check  
	#B2CIB(Go[,,,1])[[1]]-B2CIB(Go[,,,2])[[1]]
      #eigen(B2CIB(Go[,,,1])[[1]][,,1])
      #abs(eigen(B2CIB(Go[,,,2])[[1]][,,1])$values)
          
  } else CIVARB <-NA



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
  if (anyNA(SV)) {sv = Y[,SESVI]*0; SV=NA} else { sv = SV }
  ss = matrix(0,n,1)

  for ( tt in (P+1):T )  { 
      if (anyNA(SV)) {sv[2:T,] = Y[2:T,SESVI]-Y[1:(T-1),SESVI] ; SV=NA} else { sv = SV }
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
  result=list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d,r,CIVARB)
  names(result) = c("Y","X","Uo","resid","Go","GDC","Co","Sigmao","TH","St","SV","SESVI","Ao","Bo","check","type","m","n","p","S","W","SigmaS","Yo","d","r","CIVARB")
  return(result)
}




#' Data generating process of MRCIGVAR(m,n,p,S) 
#'
#' This function will generate data from an multi-regime cointegrated Global VAR(p) process and return a list containing data and parameters used in the MRCIGVAR(m,n,p,S) process.
#'
#' @param m     : number of variables in each a country/unit
#' @param n     : number of countries/units
#' @param S     : number of regimes
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
#' @param r     : an n-vector containing the number of unit root process in each country
#' @param r_np        : (m, Pmax, n, S) array containing the roots of the characteristic functions in the lag operator. 
#' @param Ncommtrend  : number of common stochastic trends in the system.
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
#'
#' @examples 
#' ## case of n = 2, m = 2, S = 2
#' #### for MRCIGVAR
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 2; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#' r  = rep(1,n)
#' i = 1
#' S = 2
#' 
#'
#' ## case of n = 3, m = 2, S = 2
#' 
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=200,S=2, SESVI=c(1,3,5),r=rep(1,3))    ## m: number of variables, n: number of countries 
#' max(res_d$Y)
#' plot(ts(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' max(abs(res_d$Y))
#' #Bo = res_d$Bo
#' #Ao = res_d$Ao
#' #res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=10000,S=2, SESVI=c(1,3,5),Bo=Bo,Ao=Ao,r=rep(1,3))    ## m: number of variables, n: number of countries 
#' max(res_d$Y)
#' res_e  = MRCIGVARest(res=res_d) 
#'
#'
#' @export
MRCIGVARData = function(m,n,p,T,S,W=NA,SESVI=NA,TH=NA,Go=NA,Ao=NA,Bo=NA,Sigmao=NA,Uo=NA,SV=NA,type=NA,Co=NA,X=NA,Yo=NA,d=NA,r=NA,r_np=NA,Ncommtrend=NA) {
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
### X       : a T x k x n x S array of exogenous veriables which may be common/different to all countries and different across all states 
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

  if (missing(TH))      	{TH = NA} 
  if (missing(W))       	{W = NA}
  if (missing(SV))      	{SV = NA}
  if (missing(Sigmao))  	{Sigmao=NA}
  if (missing(Uo))      	{Uo = NA}
  if (missing(type))    	{type = NA}
  if (missing(Co))      	{Co = NA}
  if (missing(Go))      	{Go = NA}
  if (missing(Bo))      	{Bo = NA}
  if (missing(Ao))      	{Ao = NA}
  if (missing(X))       	{X = NA}
  if (missing(d))       	{ d = NA}
  if (missing(Yo))      	{Yo = NA}
  if (missing(r))       	{r = NA}
  if (missing(r_np))    	{r_np = NA}
  if (missing(Ncommtrend))	{Ncommtrend = NA}


  if (anyNA(d))  { d = 1 }
  if (anyNA(r))  { r = rep(1,n) }

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
        GVARD = GVARData(m,n,p[,,s],T) 
        Sigmao[,,s] = GVARD$Sigmao
      }    

  }


  if (anyNA(r_np)) {
	r_np = (1:(m*Pmax*n*S))/(1:(m*Pmax*n*S))*(1.1+abs(rnorm(m*Pmax*n*S)))
	dim(r_np) = c(m,Pmax,n,S)
	r_np[1,,,] = 1
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
        k          = dim(X)[2]+1  # wrong
        CDC      = rnorm(m*k*n*S)
        dim(CDC) = c(m,k,n,S) 

        if (anyNA(Co)) Co = CDC
        GDC        = Co
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
        k          = dim(X)[2]+1   # wrong
        CDC      = rnorm(m*k*n*S)
        dim(CDC) = c(m,k,n,S)
      CDC[,1,,] = 0 
        if (anyNA(Co)) Co = CDC
        GDC        = Co
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
      CIVARB =  (1:(2*m*2*m*Pmax*S*n))*0; dim(CIVARB)=c(2*m,2*m,Pmax,n,S)
        Go        = (1:(m*n*m*n*Pmax*S))*0
        dim(Go) = c(m*n,m*n,Pmax,S)
        Bo        = (1:(m*m*Pmax*n*S))*0
        dim(Bo) = c(m,m,Pmax,n,S)

        Bo = VARBS_commtrend(m,p,T,r,Ncommtrend,n,S,r_npo=r_np)[[1]]
        if (anyNA(Ao)) Ao = Bo2Ao(m,p,T,r, Ncommtrend,n,S,r_npo=r_np)
        BoAo_h = BoAo2TBoAo(Bo,Ao,W)
        Bo = BoAo_h[[1]]
        Ao = BoAo_h[[2]]
        
        ##### country-wise transformation step 
        for (i in 1:n) {
		Bi =  matrix(rnorm(m*m),m,m)

              
	  }



        for (s in 1:S ) Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,Pmax)
        ## check  
        #B2CIB(Go[,,,1])[[1]]-B2CIB(Go[,,,2])[[1]]
      #eigen(B2CIB(Go[,,,1])[[1]][,,1])
      #abs(eigen(B2CIB(Go[,,,2])[[1]][,,1])$values)
          
  } else CIVARB <-NA



  if (!anyNA(Go) & anyNA(Bo))        {
      Bo          = (1:(m*m*Pmax*n*S))*0
        dim(Bo) = c(m,m,Pmax,n,S)
        Ao      = Bo*0
      for (s in 1:S)  {
        Bo[,,,,s] = GW2BoAo(Go[,,,s],W)$Bo
                Ao[,,,,s] = GW2BoAo(Go[,,,s],W)$Ao
      }
  }

  if (anyNA(Go) & !anyNA(Bo)) {
      Go          = (1:(m*n*m*n*Pmax*S))*0
      dim(Go) = c(m*n,m*n,Pmax,S)

        for (s in 1:S ) Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,Pmax)
  }
        
  
  St = matrix(0,T,n)
  Y = Uo[,,1]
  if (!anyNA(Yo)) Y[1:P,] = Yo
  resid = Y*0
  Yo    = Y[1:P,]
  if (anyNA(SV)) {sv = Y[,SESVI]*0; SV=NA} else { sv = SV }
  ss = matrix(0,n,1)

  for ( tt in (P+1):T )  { 
      if (anyNA(SV)) {sv[2:T,] = Y[2:T,SESVI]-Y[1:(T-1),SESVI] ; SV=NA} else { sv = SV }
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
  SigmaS = NA;                          # SigmaS is a (nmS x nmS) matrix used to construct state-dependent covariance matrices
  result=list(Y,X,Uo,resid,Go,GDC,Co,Sigmao,TH,St,SV,SESVI,Ao,Bo,check,type,m,n,p,S,W,SigmaS,Yo,d,r,CIVARB)
  names(result) = c("Y","X","Uo","resid","Go","GDC","Co","Sigmao","TH","St","SV","SESVI","Ao","Bo","check","type","m","n","p","S","W","SigmaS","Yo","d","r","CIVARB")
  return(result)
}

#' @export
VARBS_commtrend = function(m,p,T,r, Ncommtrend,n,S,r_npo) {  
  alpha   = list()
  beta    = list()
  VAR_Ione = list()
  if (missing(Ncommtrend))  {Ncommtrend=NA}
  Pmax = max(p[,1:2,])
  Pmin = min(p[,1:2,][p[,1:2,]>0])
  Bo = (1:(m * m * Pmax * n * S)) * 0
  dim(Bo) = c(m, m, Pmax, n, S)
  
  
  for (i in 1:n) {
      Pmini             = min(p[i,1:2,][p[i,1:2,]>0])
	r_npi 		= c(1:(r[i]*Pmini))/c(1:(r[i]*Pmini))*1.5
      dim(r_npi) 		= c(r[i],Pmini)
      r_npi[1:r[i],1]  	= 1  
	VAR_Ione[[i]]  	= VARData(n=r[i],p=Pmini,T=100,r_np=r_npi)
  }



  if (!anyNA(Ncommtrend)) {
     	r_npc  = c(1:(Ncommtrend*Pmin))/c(1:(Ncommtrend*Pmin))*1.5
      dim(r_npc) = c(Ncommtrend,Pmin)
      r_npc[1:Ncommtrend,1]  = 1
      VAR_IoneC = VARData(n=Ncommtrend,p=Pmin,T=100,r_np=r_npc)
   	for (i in 1:n) {
		VAR_Ione[[i]]$B[1:Ncommtrend,1:Ncommtrend,]=0		
            VAR_Ione[[i]]$B[1:Ncommtrend,1:Ncommtrend,1:Pmin] = VAR_IoneC$B			
	}
  }
  	     
      for (i in 1:n) {
 		NIzero = m - r[i]
         	#B =  matrix(rnorm(m*m),m,m)                 ### jizhu 0 nadiao
      	#B[(NIzero+1):m,] = 0 
      	#B = B + diag(m)
            B = diag(m)
            #B = eigen(t(B)%*%B)$vectors
            alphai = list()
            betai = list()
            for (s in 1:S) {
      		r_np = c(1:(NIzero * p[i,1,s]))*0      
            	dim(r_np) = c(NIzero,p[i,1,s])
            	r_np <- r_npo[(r[i]+1):m,1:p[i,1,s],i,s]
            	dim(r_np) = c(NIzero,p[i,1,s])
                  #A = matrix(rnorm(NIzero*NIzero),NIzero,NIzero); A = eigen(t(A)%*%A)[[2]]
	            VAR_IZero = VARData(NIzero, p[i, 1,s],T=100,r_np=r_np,)
            
            	for (j in 1: p[i,1,s]) {
               		Dblock    = matrix(0,m,m)
               		Dblock[1:NIzero,1:NIzero]     = VAR_IZero$B[,,j]
               		if (dim(VAR_Ione[[i]]$B)[3]>=j) Dblock[(NIzero+1):m,(NIzero+1):m] = VAR_Ione[[i]]$B[,,j]
               		Bo[, , j, i,s] = B%*%Dblock%*%solve(B)
      		}
            	#alpha[[i]] = B2CIB(VARD$B)[[2]]
            	#beta[[i]]  = B2CIB(VARD$B)[[3]]
			B2CIBab  	= B2CIB(Bo[, , 1:p[i, 1,s], i,s])
            	alphai[[s]] = B2CIBab[[2]]
            	betai[[s]]  = B2CIBab[[3]]
      	}
            alpha[[i]] = alphai
		beta[[i]]  = betai
  	}



  return(list(Bo,alpha,beta))
}


############################################
#' @export
Bo2Ao = function(m,p,T,r, Ncommtrend,n,S,r_npo) {
    d = 1 
    Boab  = VARBS_commtrend(m,p,T,r,Ncommtrend,n,S,r_npo)
    Pmax = max(p[,1:2,])
    beta  = Boab[[3]]
    alpha = Boab[[2]]
    Ao = (1:(m * m * Pmax * n * S)) * 0
    dim(Ao) = c(m, m, Pmax, n, S)
    for (i in 1:n) {
        for (s in 1:S) { 
        	if (p[i,2,s] < 2) Ao = Ao  
        	if (p[i,2,s] >= 2) {
                VARD = VARData(m, p=(p[i, 2,s]-1),T)
            	BB = VARD$B
            	Ao[,,1,i,s] = BB[,,1]+alpha[[i]][[s]]%*%t(beta[[i]][[s]])*d
            	Ao[,,p[i,2,s],i,s] = -BB[,,p[i,2,s]-1]
            	if ((p[i, 2,s]-1)>=2) { for (L in 2:(p[i, 2,s]-1)) Ao[, ,L, i,s] = BB[,,L]-BB[,,L-1]}
        	}
        }
    }
    
    return(Ao)
}

#' @export
BoAo2TBoAo = function(Bo,Ao,W) {
      TBo = Bo*0
      TAo = Ao*0
	n = dim(W)[1]
 	m = dim(Bo)[1]
      L = dim(Bo)[3]
      BD = matrix(0,m*n,m*n)
      B =  matrix(rnorm(m*m),m,m) 
      for (i in 1:n) {
		BD[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)] = B		
	}
      for (j in 1:L) {
	   for ( s in 1: S ) {

      	BoDj =  matrix(0,m*n,m*n)
		AoDj =  matrix(0,m*n,m*n)
		
		for (i in 1:n)  {
    			BoDj[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]=Bo[,,j,i,s]
			AoDj[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]=Ao[,,j,i,s]
      	}
		B_h = BD%*%BoDj%*%solve(BD)
            A_h = BD%*%AoDj%*%kronecker(W,diag(m))%*%solve(BD)%*%solve(kronecker(W,diag(m)))
		for (i in 1:n)  {
    			TBo[,,j,i,s] = B_h[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]
			TAo[,,j,i,s] = A_h[((i-1)*m+1):((i-1)*m+m),((i-1)*m+1):((i-1)*m+m)]
      	}
         }				
	}
	return(list(TBo,TAo))
}








#' Estimation of MRCIGVAR(m,n,p,S,r) 
#'
#' This function will estimate the parameters of a multi regime cointegrated VAR model. The adjustment speeds are assumed to be identical in different regimes.
#'
#' @param res   : an MRCIVAR object of the output of MRCIGVARData
#' @return      : an MRCIVAR object with estimated parameters and statistics.
#' @examples 
#' ## case of n = 2, m = 2, S = 2
#' #### for MRCIGVAR
#' m = 2
#' n = 3
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:n,,]; p[,1,] = 2; p[,2,] = 2
#'
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI=rep(1,3,5)
#' r  = rep(1,n)
#' i = 1
#' S = 2
#' 
#'
#' ## case of n = 3, m = 2, S = 2
#' 
#' res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=200,S=2, SESVI=c(1,3,5),r=rep(1,3))    ## m: number of variables, n: number of countries 
#' max(abs(res_d$Y))
#' res_e  = MRCIGVARest(res=res_d) 
#' res_e$CRK_tst
#'
#'
#' @export
MRCIGVARest=function(res) {
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
 

    m = res$m
    n = res$n
    p = res$p
    Y = as.matrix(res$Y)
    X = res$X
    W = res$W
    type = res$type
    TH = res$TH
    SESVI = res$SESVI
    Go = res$Go
    S = res$S
    Ao = res$Ao
    Bo = res$Bo
    Co = res$Co
    GDC = res$GDC
    Sigmao = res$Sigmao * 0
    d = res$d
    P  = max(p[,1:2,],d)
    Pmax = max(p[,1:2,])
    St = res$St
    r  = res$r
    Uo = res$Uo
    SigmaS = (1:(n*m*S*n*m*S))*0
    dim(SigmaS) = c(n*m*S,n*m*S)
    VAR_domestic = list()
    CRK_tst = list()

    T = nrow(Y)
    check = c(1:(S+1))*0

  Ct = Uo*0;

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
 
      WMATn = W%x%diag(m)
      FY      = Y%*%WMATn
      
  	Go 	  = (1:(m*n*m*n*Pmax*S))*0
      dim(Go) = c(m*n,m*n,Pmax,S)
	Bo 	  = (1:(m*m*Pmax*n*S))*0
	dim(Bo) = c(m,m,Pmax,n,S)
	Ao      = Bo*0
	if (type=="none")   Model="I"
      if (type=="const")  Model="III"
            for (i in 1:n) {
                  pi = t(p[i,,])
                  y = Y[,((i-1)*m+1):(i*m)]
                  x = FY[,((i-1)*m+1):(i*m)]
                  crk = m-r[i]
                  Sti = 2-St[,i]
			#CIVARD = MRCIVARData(n=2*m,p=pi,T =100,S=S,SESVI=SESVI[i],TH=TH[i],Sigmao=sigma,type="none",r=r[i])
			tst <- MRCVECMest2(y,x,s=Sti,model=Model,type="trace",P=pi,crk=crk,q = 0.95)
			VAR_domestic[[i]] = tst[[2]]
 			CRK_tst[[i]] = tst[[1]]
			
			nrowresid = nrow(VAR_domestic[[i]]$residuals)
                  Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),1] = (2-St[(T-nrowresid+1):T,i])*VAR_domestic[[i]]$residuals
                  Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),2] = (St[(T-nrowresid+1):T,i]-1)*VAR_domestic[[i]]$residuals
                  
			AB1 = VECM2VAR(param = tst[[2]][[1]][1:(crk+m*(pi[1,1]-1+pi[1,2]-1)),],beta=tst$beta,p=c(crk,(pi[1,1]-1),(pi[1,2]-1)))
			AB2 = VECM2VAR(param = tst[[2]][[1]][c(1:crk,(crk+m*(pi[1,1]-1+pi[1,2]-1))+1:(m*(pi[2,1]+pi[2,2]-2))),],beta=tst$beta,p=c(crk,(pi[2,1]-1),(pi[2,2]-1)))

			Bo[,,1:pi[1,1],i,1] = AB1[[1]]
		      Ao[,,1:pi[1,2],i,1] = AB1[[2]]
			Bo[,,1:pi[2,1],i,2] = AB2[[1]]
		      Ao[,,1:pi[2,2],i,2] = AB2[[2]]
			if (type=="const") {
				Co[,,i,1] = tst[[2]][[1]][crk+(sum(pi[,-3])-4)*m+1,] ;
				Co[,,i,2] = tst[[2]][[1]][crk+(sum(pi[,-3])-4)*m+2,]
			}

            }
           

  	for (s in 1:S ) Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,Pmax)
 	
      for (s in 1:S) Sigmao[,,s] = t(Uo[,,s])%*%Uo[,,s]/sum(St[,i]==s)

      for (s in 1:S) {
		for (i in 1:n) {
			for (j in 1:n)  Sigmao[(m*(i-1)+1):(i*m),(m*(j-1)+1):(j*m),s] = t(Uo[,(m*(i-1)+1):(i*m),s])%*%Uo[,(m*(j-1)+1):(j*m),s]/sum((St[,i]==s)*(St[,j]==s))

		}
	}
 
      for (s in 1:S)  { 
          for (ss in 1:S)  {
               for (i in 1:n)  {
                   for (j in 1:n) SigmaS[ ((s-1)*n*m+(i-1)*m+1):((s-1)*n*m+i*m), ((ss-1)*n*m+(j-1)*m+1):((ss-1)*n*m+j*m)] =t(Uo[,(m*(i-1)+1):(i*m),s])%*%Uo[,(m*(j-1)+1):(j*m),ss]/sum((St[,i]==s)*(St[,j]==ss))

     
		   }
	    }
	}




    res$Go = Go
    res$Sigmao = Sigmao
    res$r_npo = NA
    res$Ao = Ao
    res$Bo = Bo
    res$Co = Co
    res$GDC = GDC
    res$SigmaS = SigmaS
    res$resid = Uo
    res$VAR_domestic = VAR_domestic
    res$CRK_tst      = CRK_tst
    return(res)
}





#' Estimation of MRCIGVAR(res) 
#'
#' This function will estimate the parameters of a multi regime cointegrated VAR model. The adjustment speeds to the cointegration relations are assumed to be different in different regimes.
#'
#' @param res   : an MRCIVAR object of the output of MRCIGVARData
#' @return      : an MRCIVAR object with estimated parameters and statistics.
#' @examples 
#' ## case of n = 2, m = 2, S = 2
#' m = 3
#' S = 2 
#' n  = 10
#' p = (1:(3*n*S))*0
#' dim(p) = c(n,3,S)
#' p[,1:2,] = 2; 
#' TH = c(1:n)*0; dim(TH) = c(1,n)
#' SESVI= (1:n)*m-m+1
#' r  = rep(1,n)
#' 
#' res_d <- MRCIGVARData(m=3,n=10,p=p,TH=TH,T=200,S=2, SESVI=(1:10)*3-2,r=rep(1,10),Ncommtrend=1)    ## m: number of variables, n: number of countries 
#' max(abs(res_d$Y))
#' STAT(res_d$Go[,,,2])
#' STAT(res_d$Go[,,,1])
#' 
#' res_e  = MRCIGVARest(res=res_d) 
#' res_em = MRCIGVARestm(res=res_d) 
#' summary(res_e$VAR_domestic[[1]])
#' summary(res_em$VAR_domestic[[1]])
#'
#' @export
MRCIGVARestm=function(res) {
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
 

    m = res$m
    n = res$n
    p = res$p
    Y = as.matrix(res$Y)
    X = res$X
    W = res$W
    type = res$type
    TH = res$TH
    SESVI = res$SESVI
    Go = res$Go
    S = res$S
    Ao = res$Ao
    Bo = res$Bo
    Co = res$Co
    GDC = res$GDC
    Sigmao = res$Sigmao * 0
    d = res$d
    P  = max(p[,1:2,],d)
    Pmax = max(p[,1:2,])
    St = res$St
    r  = res$r
    Uo = res$Uo
    SigmaS = (1:(n*m*S*n*m*S))*0
    dim(SigmaS) = c(n*m*S,n*m*S)
    VAR_domestic = list()

    CRK_tst = list()
    T = nrow(Y)
    check = c(1:(S+1))*0
    LH_P   = (1:n)*0
    LH_AIC = (1:n)*0
    LH_BIC = (1:n)*0
    ORLHP  = (1:n)*0
    ORAIC  = (1:n)*0
    ORBIC  = (1:n)*0

  Ct = Uo*0;

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
 
      WMATn = W%x%diag(m)
      FY      = Y%*%WMATn
      
  	Go 	  = (1:(m*n*m*n*Pmax*S))*0
      dim(Go) = c(m*n,m*n,Pmax,S)
	Bo 	  = (1:(m*m*Pmax*n*S))*0
	dim(Bo) = c(m,m,Pmax,n,S)
	Ao      = Bo*0
	if (type=="none")   Model="I"
      if (type=="const")  Model="III"
            for (i in 1:n) {
                  pj = t(p[i,,])
                  pjOR = pj*0
                  pjOR[1,1]=max(pj[,1])
  			pjOR[1,2]=max(pj[,2])
            
                  y = Y[,((i-1)*m+1):(i*m)]
                  x = FY[,((i-1)*m+1):(i*m)]
                  crk = m-r[i]
                  Sti = 2-St[,i]
                  Pmaxi = max(p[i,1,])
			#CIVARD = MRCIVARData(n=2*m,p=pj,T =100,S=S,SESVI=SESVI[i],TH=TH[i],Sigmao=sigma,type="none",r=r[i])

                  ORtst  <- MRCVECMest2(y, x,model=Model,type = "trace",P =pjOR, crk = crk , q = 0.95, Dxflag = 0)


                  nrowresidOR = nrow(ORtst[[2]]$residuals)
                  UOR = ORtst[[2]]$residuals
                  SigmaOR = t(UOR)%*%UOR/(T-P*m+crk-ncol(x)+(Model=="III"))
                  ORLHP[i] = -(T*m/2)*log(2*pi) -(T*m/2) + T/2*log(det(solve(SigmaOR)))
                  ORAIC[i] =  2*m*(m*pjOR[1,1]-m+crk+(m+1)/2) - 2*ORLHP[i]
                  ORBIC[i] =  log(T)*m*(m*pjOR[1,1]-m+crk+(m+1)/2)  - 2*ORLHP[i]

    			tst <- MRCVECMestm(y,x,s=Sti,model=Model,type="trace",P=pj,crk=crk,q = 0.95)
			VAR_domestic[[i]] = tst$VECM1
                  CRK_tst[[i]] = tst[[1]]

			nrowresid = nrow(VAR_domestic[[i]]$residuals)
                  Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),1] = (2-St[(T-nrowresid+1):T,i])*VAR_domestic[[i]]$residuals
                  Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),2] = (St[(T-nrowresid+1):T,i]-1)*VAR_domestic[[i]]$residuals
                  

                  #Sigmao[(m*(i-1)+1):(i*m),(m*(i-1)+1):(i*m),1] = t(Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),1])%*%Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),1]/(sum(2-St)-(m*sum(p[i,1:2,1])+m-crk)) 
                  #Sigmao[(m*(i-1)+1):(i*m),(m*(i-1)+1):(i*m),2] = t(Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),2])%*%Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),2]/(sum(St-1)-(m*sum(p[i,1:2,2])+m-crk)) 
			sigma1 = t(Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),1])%*%Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),1]/(sum(2-St)-(m*sum(p[i,1:2,1])+m-crk)) 
                  sigma2 = t(Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),2])%*%Uo[(T-nrowresid+1):T,((i-1)*m+1):(i*m),2]/(sum(St-1)-(m*sum(p[i,1:2,2])+m-crk)) 
             
    			LH_P[i]   = -(T*m/2)*log(2*pi) -(T*m/2) + (sum(2-St))/2*log(det(solve(sigma1)))+(sum(St-1))/2*log(det(solve(sigma2)))
                  LH_AIC[i] =  2*(m*(sum(p[i,1:2,1])-m+crk) +m*(m+1)/2) + 2*(m*(sum(p[i,1:2,2])-m+crk)+m*(m+1)/2) - 2*LH_P[i]
                  LH_BIC[i] =  log(sum(2-St))*(m*(sum(p[i,1:2,1])-m+crk) +m*(m+1)/2)+ log(sum(St-1))*(m*(sum(p[i,1:2,2])-m+crk)+m*(m+1)/2)  - 2*LH_P[i]


			AB1 = VECM2VAR(param = tst[[2]][[1]][1:(crk+m*(pj[1,1]-1+pj[1,2]-1)),],beta=tst$beta,p=c(crk,(pj[1,1]-1),(pj[1,2]-1)))
			AB2 = VECM2VAR(param = tst[[2]][[1]][c(1:crk,(crk+m*(pj[1,1]-1+pj[1,2]-1))+1:(m*(pj[2,1]+pj[2,2]-2))),],beta=tst$beta,p=c(crk,(pj[2,1]-1),(pj[2,2]-1)))

			Bo[,,1:pj[1,1],i,1] = AB1[[1]]
		      Ao[,,1:pj[1,2],i,1] = AB1[[2]]
			Bo[,,1:pj[2,1],i,2] = AB2[[1]]
		      Ao[,,1:pj[2,2],i,2] = AB2[[2]]
			if (type=="const") {
				Co[,,i,1] = tst[[2]][[1]][crk+(sum(pj[,-3])-4)*m+1,] ;
				Co[,,i,2] = tst[[2]][[1]][crk+(sum(pj[,-3])-4)*m+2,]
			}

            }
           

  	for (s in 1:S ) Go[,,,s] = BoAoW2G(Bo[,,,,s],Ao[,,,,s],W,m,n,Pmax)
 	
      for (s in 1:S) Sigmao[,,s] = t(Uo[,,s])%*%Uo[,,s]/sum(St[,i]==s)

      for (s in 1:S) {
		for (i in 1:n) {
			for (j in 1:n)  Sigmao[(m*(i-1)+1):(i*m),(m*(j-1)+1):(j*m),s] = t(Uo[,(m*(i-1)+1):(i*m),s])%*%Uo[,(m*(j-1)+1):(j*m),s]/sum((St[,i]==s)*(St[,j]==s))

		}
	}
 
      for (s in 1:S)  { 
          for (ss in 1:S)  {
               for (i in 1:n)  {
                   for (j in 1:n) SigmaS[ ((s-1)*n*m+(i-1)*m+1):((s-1)*n*m+i*m), ((ss-1)*n*m+(j-1)*m+1):((ss-1)*n*m+j*m)] =t(Uo[,(m*(i-1)+1):(i*m),s])%*%Uo[,(m*(j-1)+1):(j*m),ss]/sum((St[,i]==s)*(St[,j]==ss))

     
		   }
	    }
	}
    
   


    res$Go = Go
    res$Sigmao = Sigmao
    res$r_npo = NA
    res$Ao = Ao
    res$Bo = Bo
    res$Co = Co
    res$GDC = GDC
    res$SigmaS = SigmaS
    res$resid = Uo
    res$VAR_domestic = VAR_domestic
    res$CRK_tst      = CRK_tst
    res$ORAIC = sum(ORAIC)
    res$ORBIC = sum(ORBIC)
    res$ORLHP = sum(ORLHP)
    res$LH_BIC = sum(LH_BIC)
    res$LH_AIC = sum(LH_AIC)
    res$LH_P  = sum(LH_P)
    return(res)
}









#' Calculation of the model selection values of for MRCIGVAR models  
#' 
#' @param  res  : an object generated from CIGVARData or estimated from CIGVARest.
#' @param  L_V  : a 2 components vector containing the maxima of the domestic lag and the foreign lag respectively. 
#' @param  TH_V  : a vector containing the maxima of the domestic lag and the foreign lag respectively. 

#' @return      : A matrix with different lag specifications and the model selection criteria.
#' @examples 
#' 
#' m = 2
#' n = 4
#' p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
#' p = p[1:4,,]; p[,1:2,2] = 3
#' 
#' TH = c(1:4)*0; dim(TH) = c(1,4)
#' SESVI=c(1,1)
#' r  = c(1,1,1,1)
#' i = 1
#' S = 2
#' 
#' res_d <- MRCIGVARData(m=2,n=4,p=p,TH=TH,T=4000,S=2, SESVI=rep(1,4),r=rep(1,4))    ## m: number of variables, n: number of countries 
#' 
#' res_e = MRCIGVARest(res=res_d) 
#' res_em = MRCIGVARestm(res=res_d) 
#'
#' L_v = c(2,3)
#' TH_v = c(0,0.1)
#' res_d$p
#' 
#' MRCIGVARSelect = MRCIGVAR_Select(res=res_d,L_V=L_v,TH_V=TH_v)
#' dim(MRCIGVARSelect)
#' MRCIGVARSelect[which.min(MRCIGVARSelect[,28]),]
#'
#' @export
MRCIGVAR_Select = function(res,L_V=L_v,TH_V=TH_v)  {
Criteria = matrix(0,n*(L_V[1]-1)*(L_V[2]-1)*2*length(TH_V),length(as.vector(res$p))+6)
idx = 0
res_dd   = res
idx = 0
n   = res$n

for (i in 1:n )           	{
for (l_d in 2: L_V[1] )   	{
for (l_f in 2: L_V[2] )   	{
for (s  in 1:2       )    	{
for (l_th in 1: length(TH_V) ){ 
      idx = idx + 1
      res_dd$p[i,1,s] = l_d 
      res_dd$p[i,2,s] = l_f       
	res_dd$TH       = TH_V[l_th]

    	res_s = MRCIGVARestm(res_dd)
      Criteria[idx,] = c(c(as.vector(res_dd$p),TH_V[l_th],s),sum(res_s$LH_AIC),sum(res_s$LH_BIC),sum(res_s$ORAIC),sum(res_s$ORBIC))
      #colnames(Criteria) = c("l_d","l_f","AIC","BIC","ORAIC","ORBIC")
    }
 }
}
}
}
return(Criteria)
}


#' Regime specific impulse response functions of MRCIGVAR(m,n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRCIGVAR(n,p,S). Using G[,,,s] and Sigma[,,s] matrices of the estimated MRCIGVAR, this function can produce impulse response functions for any possible combinations of states. 
#' @param res a list of estimated MRGVAR as output of MRGVARest 
#' @param nstep the length of impulse response function 
#' @param comb a vector specify the concerted action in policy-simulation impulse response function 
#' @param state an n vector specifying the speciic state for each country. 
#' @param irf types of the impulse response function c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol cholezky decomposition, chol1 cholezky decomposition with unit impulse, comb1 concerted action with unit impulse.
#' @return a matrix of (mn,mn,nstep) as the IRF colunms respesenting the impulse and rows the responses. 
#' @examples 

#'
#' @export
#'
irf_MRCIGVAR_CB = function(res,nstep,comb,state=c(2,1),irf=c("gen","chol","chol1","gen1","comb1"),runs=200,conf=c(0.05,0.95),NT=1) {
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
      Sigmao = res$Sigmao
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
      sigmaNPD1 = SigmaNPD(res,rep(1,n))
 	sigmaNPD2 = SigmaNPD(res,rep(2,n))
      Sigmao_sim = Sigmao*0
      Sigmao_sim[,,1] = sigmaNPD1
      Sigmao_sim[,,2] = sigmaNPD2
      response[,,,1] <- irf_MRCIGVAR(res,nstep,comb,state,irf)
      responseR <- array(0,dim=c(neq,nvar,nstep,runs))

      Uo_run = array(0,c(T,n*m,S))

      for (i in 1:runs) {
            for (s in 1:S)    Uo_run[,,s] = rnormSIGMA(T,sigmaNPDS)
            if (length(p)>1)	{  
               		#res_run = MRGVARData3(m,n,p,T,S,W,SESVI,TH,res$Go,Ao=NA,Bo=NA,Sigmao=NA,Uo=Uo_run,SV,type,Co=NA,X)
               		res_run  = MRCIGVARData(m,n,p,T,S,W,SESVI,TH,res$Go,Sigmao=Sigmao_sim,SV,type="none")

			  	res_e    = MRCIGVARest(res_run)
		}

            if (length(p)==1) {
				res_run = MRGVARData2(m,n,p,T,S,W,SESVI,TH,res$Go,Ao=NA,Bo=NA,Sigmao=NA,Uo=Uo_run,SV,type,res$Co,X)
                        #(min(colSums(res_run$Y[,seq(2,68,3)]<1.6))> 20) $ (max(colSums(res_run$Y[,seq(2,68,3)]<1.6)))> 20)  
            		res_e   = MRCIGVARest(res_run)
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


#' Regime specific impulse response functions of MRCIGVAR(n,p,S)
#'
#' This function calculates the regime specific impulse response functions of an estimated MRCIGVAR(n,p,S). 
#' Using the estimated G[,,,s] and Sigma[,,s] matrices of the MRGVAR, this function calculated the regime speicfic impulse response functions.
#' @param res a list of estimated MRCIGVAR as output of MRCIGVARest or MRCIGVARestm 
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
#' IRFF2 = irf_MRCIGVAR(res=res_e,nstep=10,comb=NA,state=c(2,2,2,2,2),irf='gen')
#' plot(IRFF2[1,1,],type="l")
#' 
#'
#' @export
#'
irf_MRCIGVAR = function(res,nstep,comb,state=state,irf=c("gen","chol","chol1","gen1","comb1"),sigmaNPDS=NA) {
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
      if (anyNA(sigmaNPDS))	sigma = SigmaNPD(res,state)  else  sigma = sigmaNPDS
      response <- array(0,dim=c(neq,nvar,nstep));
      response <- impulsdtrf_g(B,sigma,nstep,comb,irf=irf)
	return(response)
} 




