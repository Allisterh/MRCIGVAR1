#' Data generating process of CIVAR(p) 
#'
#' This function will generate data from a cointegrated vector autoregressive process CIVAR(p) and return a list containing data and parameters used in the CIVAR(p) process.
#'
#' @param n     : number of variables
#' @param p     : lag length
#' @param T     : number of observations
#' @param r_np  : n x p matrix of roots outside and on the unit circle for an n dimensional independent ar(p)-processes (the roots of the characteristic functions in the lag operator.). An element one in a row implies a unit root process. If the matrix is not provided, it will be generated randomly with one unit root in the first row.
#' @param A     : an n x n full rank matrix of transformation to generate correlated a VAR(p) from the n independent AR(p)  
#' @param B     : (n,n,p) array of the CIVAR(p) process. If B is not given, it will be calculated out of r_np and A.
#' @param Co    : (n,k+1) vector of intercept of the CIVAR(p) process. for type="none" Co = O*(1:n), for const Co is an n vector, exog0 Co is a (n,k) matrix for exog1 Co is an (n,1+k) matrix.  
#' @param U     : residuals, if it is not NA it will be used as input to generate the VAR(p)
#' @param Sigma : The covariance matrix of the CIVAR(p)
#' @param type  : deterministic component "none", "const" "exog0" and "exog1" are four options
#' @param X     : (T x k) matrix of exogeneous variables.
#' @param mu    : an n vector of the expected mean of the VAR(p) process
#' @param Yo    : (p x n) matrix of initial values of the process
#' @param crk   : the numer of cointegration relations. It equals n-r, where r is the number of unit roots.
#' @return      : An object of CIVAR containing the generated data, the parameters and the input of exogeous variables. res = list(n,p,type,r_np,Phi,A,B,Co,Sigma,Y,X,resid,U,Y1,Yo,check)
#' @examples 
#' res_d = CIVARData(n=2,p=2,T=100,type="const")  
#' res_d = CIVARData(n=2,p=2,T=10,Co=c(1:2)*0,type="none") 
#' res_d = CIVARData(n=2,p=2,T=10,Co=c(1:2)*NA,  type="const") 
#'
#' p = 3
#' n = 4
#' r_np = matrix(c(1,2,1.5,1.5,2.5,2.5,2,-1.5,2,-4,1.9,-2.1),4,3)
#' 
#' res_d  = CIVARData(n=4,p=3,T=200,r_np=r_np,Co=matrix(0,n,1) ,type="none",crk=3)
#' res_e = CIVARest(res=res_d)
#' sum(abs(res_e$B-res_d$B))
#' sum(abs(res_e$Co-res_d$Co))
#'
#'
#' res_d  = CIVARData(n=4,p=3,T=200)
#' res_e = CIVARest(res=res_d)
#' sum(abs(res_e$B-res_d$B))
#' sum(abs(res_e$Co-res_d$Co))#' plot(ts(res_d$Y))
#' @export
CIVARData=function (n, p, T, r_np, A, B, Co, U, Sigma, type, X, mu, Yo,crk) 
{
    if (missing(r_np)) {
        r_np = NA
    }
    if (missing(A)) {
        A = NA
    }
    if (missing(B)) {
        B = NA
    }
    if (missing(Co)) {
        Co = NA
    }
    if (missing(U)) {
        U = NA
    }
    if (missing(Sigma)) {
        Sigma = NA
    }
    if (missing(type)) {
        type = NA
    }
    if (missing(X)) {
        X = NA
    }
    if (missing(mu)) {
        mu = NA
    }
    if (missing(Yo)) {
        Yo = NA
    }
    if (missing(crk)) {
        crk = NA
    }

    if (anyNA(r_np)) {
        r_np = matrix(0, n, p)
        for (i in 1:n) r_np[i, ] <- 0.51/(runif(p) - 0.5)
        if (anyNA(crk)) {r_np[1,1] = 1; crk = n-1 }
          else {
            if (crk<n) for ( j in 1:(n-crk) ) r_np[j,1] = 1
        }
    }
    d = dim(r_np)
    if (!(n == d[1] & p == d[2])) {
        print("dimension problem")
        return("dimension")
    }
    Phi = NA
    if (anyNA(Phi)) {
        Phi = matrix(0, n, p)
        for (i in 1:n) Phi[i, ] = Roots2coef(p, r_np[i, ])
    }
    if (anyNA(A)) {
        A = matrix(runif(n * n), n, n)
        A = A %*% t(A)
        A = eigen(A)[[2]]
    }
    if (anyNA(B)) {
        B = matrix(0, n, n * p)
        dim(B) = c(n, n, p)
        if (n == 1) {
            B = Phi
            dim(B) = c(n, n, p)
        }
        if (n > 1) {
            for (i in 1:p) B[, , i] = A %*% diag(Phi[, i]) %*% 
                t(A)
        }
    }
    if (anyNA(Sigma)) {
        Sigmao = matrix(rnorm(n * n), n, n)
        Sigmao = Sigmao %*% t(Sigmao)
        Sigma = A %*% Sigmao %*% t(A)
    }
    else {
        Sigmao = solve(A) %*% Sigma %*% solve(t(A))
    }
    if (anyNA(U)) {
        Uh = rnormSIGMA(T, Sigmao)
        U = Uh %*% t(A)
    }
    else {
        Uh = U %*% solve(t(A))
    }
    if (anyNA(Yo)) {
        Yo = Uh
    }
    else {
        Uh[1:p, ] = Yo
        Yo = Uh
    }
    for (i in 1:n) {
        for (t in (p + 1):T) {
            for (L in 1:p) Yo[t, i] = Yo[t, i] + Yo[t - L, i] * 
                Phi[i, L]
        }
    }
    Y1 = Yo %*% t(A)
    Ct = U * 0
    if (anyNA(type)) {
        type = "const"
    }
    type_soll = Type(Co, X)
    if (!type_soll == type) 
        return(list("type mismatch", type_soll, type))
    if (!anyNA(Co) & (!nrow(as.matrix(Co)) == n)) 
        return("dimension mismatch")

    if ((type == "none") & (anyNA(Co))) {
        Co = matrix(0,n,1)
    }

    if (type == "const") {
        if (anyNA(mu)) {
           
            mu = matrix(rnorm(n), n, 1)
        }
        if (anyNA(Co)) {
            Co = mu
            for (L in 1:p) Co = Co - B[, , L] %*% mu
        }
        else {
            H = diag(n)
            
            for (L in 1:p) H = H - B[, , L]
            if (min(abs(eigen(H)$values))>0.0000000000001) mu = solve(H) %*% Co  else   mu = NA
            
               
        }
        Ct = matrix(1, T, 1) %*% t(Co)
    }
    if (type == "exog0" & anyNA(Co)) {
        k = ncol(as.matrix(X))
        Co = matrix(rnorm(n * (k + 1)), n, k + 1)
        Co[, 1] = 0
    }
    if (type == "exog1" & anyNA(Co)) {
        k = ncol(as.matrix(X))
        Co = matrix(rnorm(n * (k + 1)), n, k + 1)
    }
    if (!anyNA(X)) {
        if (type == "exog0") 
            Ct = X %*% t(Co[, -1])
        if (type == "exog1") 
            Ct = as.matrix(cbind((1:T)/(1:T), X)) %*% t(Co)
    }
    else {
        X = NA
    }
    if (n > 0) {
        Y = U + Ct
        for (tt in (p + 1):T) {
            for (L in 1:p) Y[tt, ] = Y[tt, ] + Y[tt - L, ] %*% 
                t(B[, , L])
        }
    }
    check = max(abs(Y1))
    resid = Uh
    result = list(n, p, type, r_np, Phi, A, B, Co, Sigma, Y, 
        X, resid, U, Y1, Yo, check, crk)
    names(result) = c("n", "p", "type", "r_np", "Phi", "A", "B", 
        "Co", "Sigma", "Y", "X", "resid", "U", "Y1", "Yo", "check", "crk")
    return(result)
}


#' Estimation of CIVAR(p) 
#'
#' This function estimates the unknown parameters of a specified CIVAR(p) model based on provided data.
#'
#' @param  res  :an object of CIVAR(p) containing the components which are the output of CIVARData including as least: n, p, Y, crk and optionally X and type. 
#' @return res  :an object oc CIVAR(p) containing estimated parameter values, AIC, BIC, LH and the estimated VECM in regression format.
#' @examples 
#' p = 3
#' n = 4
#' r_np = matrix(c(1,2,1.5,1.5,2.5,2.5,2,-1.5,2,-4,1.9,-2.1),4,3)
#'
#' res_d  = CIVARData(n=4,p=3,T=200,r_np=r_np,Co =(1:n)/(1:n)*0,type="none",crk=3)
#' res_e = CIVARest(res=res_d)
#'
#' res_d  = CIVARData(n=4,p=3,T=200)
#' B = res_d$B
#' plot(ts(res_d$Y))
#' res_d$Co
#' res_d$type
#' res_d$crk
#' res_e = CIVARest(res=res_d)
#'
#'
#'
#' @export
CIVARest = function (res) 
{
    n = res$n
    p = res$p
    crk = res$crk
    Y = as.matrix(res$Y)
    type = res$type
    T = dim(Y)[1]
    crk = res$crk
    if (type == "none") 
        Model = "I"
    if (type == "const") 
        Model = "III"
    P = matrix(0, 2, 3)
    P[, 1] = p
    tst <- MRCVECMest2(y = Y, x = 0, model = Model, type = "eigen", 
        P = P, crk = crk, q = 0.95, Dxflag = 0)
    param = tst[[2]][[1]]
    B = VECM2VAR(param = param, beta = tst$beta, p = c(crk, p - 
        1, 0))[[1]]
    C = VECM2VAR(param = param, beta = tst$beta, p = c(crk, p - 
        1, 0))[[3]]
    LREG = tst[[2]]
    sigma = t(LREG$residuals) %*% (LREG$residuals)/(T)
    resid = Y * 0
    resid[(p + 1):T, ] = LREG$residuals
    if (type == "const") 
        Co = C
    if (type == "none") 
        Co = (1:n)/(1:n) * 0
    res$B <- B
    res$Co <- Co
    res$Sigma <- sigma
    res$resid <- resid
    LH = -(T * n/2) * log(2 * pi) - (T * n/2) + ((T)/2) * log(det(solve(sigma)))
    res$LH = LH
    AIC = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1) - 2 * LH
    BIC = log(T) * (n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)/2) - 2 * LH
    res$AIC = AIC
    res$BIC = BIC
    res$LH_N = 2 * n * (dim(tst[[2]][[1]])[1]) + n * (n + 1)
    if (is.null(colnames(Y))) 
        colnames(Y) = sprintf("Y%s", 1:n)
    res$varp = LREG
    res$tst  = tst
    result = res
    return(res)
}



#' Estimation of a multi regime conditional vector error correction process. 
#'
#' This function estimates the unknown parameters of a multi regime conditional VECM based on provided data.
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
#'
#' @export
MRCVECMest2 = function (y,x,s,model = c("I","II","III","IV","V"),type = c("eigen", "trace"), constant = TRUE, ret = c("statistic", 
    "test"), ctable = c("A3", "A1", "A5"), crk=2, P = matrix(2,2,2), q = 0.95,Dxflag=0) 
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
          ZyI <- Zy[,-Aa]
      }   else ZyI = Zy

      if (P[1,2]<p) {
          Ba = P[1,2]*NN1+(1:((p-P[1,2])*NN1))    
          ZxI = Zx[,-Ba]
      }   else ZxI = Zx

      if (P[2,1]<p) {
          Aa = P[2,1]*N1+(1:((p-P[2,1])*N1))    
          ZyII = Zy[,-Aa]
      }   else ZyII = Zy
      if (P[2,2]<p) {
          Bb = P[2,2]*NN1+(1:((p-P[2,2])*NN1))    
          ZxII = Zx[,-Bb]
      }   else ZxII = Zx

      Z = cbind(ZyI,ZxI); 
      Z_2 = cbind(ZyII,ZxII)
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
    beta   <- V[,1:crk]

    if (crk>0) {
       CI = Z1%*%beta
       VECM<-lm(Y0~0+CI+Z2) 
    }
    if (crk == 0 ) {
       CI = 0
       VECM<-lm(Y0~0+Z2) 
    } 
    E = -T1 * log(1 - lambda)
    E = E[1:N]
    ##### Estimation of VECM
  
    #####
    ##### transform to level model
    #####
    ##### 
    resultsvecm <- summary(VECM)
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
        result$estimation = VECM
        result$lambda     = E
        result$z          = z
        result$Z2         = Z2
        result$beta       = beta
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
        result$estimation = VECM
        result$lambda     = E
        result$z          = z
        result$Z2         = Z2
        result$beta       = beta
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

 

#' ########
#' @export
VECM2VAR = function(param,beta,p=c(1,2,2,2,2),s=NA) {
 m = dim(param)[2]
 VECB = t(param)
 if (anyNA(s)) {
 	if (!p[3]==0) {
 		B = (1:(m*m*(p[2]+1)))*0; dim(B) = c(m,m,(p[2]+1))
 		A = (1:(m*m*(p[3]+1)))*0; dim(A) = c(m,m,(p[3]+1))
 		AB = VECB[,1:p[1]]%*%t(beta)
 		B[,,1] = AB[,1:m]
 		A[,,1] = AB[,(m+1):(2*m)]
 		for ( i in 2:(p[2]+1)) B[,,i] = VECB[, p[1]+((i-2)*m+1):((i-2)*m+m)]
 		for ( i in 2:(p[3]+1)) A[,,i] = VECB[, (p[1]+p[2]*m)+((i-2)*m+1):((i-2)*m+m)]
 		#CIB3B(B)
 		#CIA2A(A) 
      	B = CIB3B(B)
      	A = CIA2A(A) 	
 	} else {
		B = (1:(m*m*(p[2]+1)))*0; dim(B) = c(m,m,(p[2]+1))
 		AB = VECB[,1:p[1]]%*%t(beta)
 		B[,,1] = AB[,1:m]
 		for ( i in 2:(p[2]+1)) B[,,i] = VECB[, p[1]+((i-2)*m+1):((i-2)*m+m)]
 		B = CIB3B(B)
      	A = NA
  	}
  	if ( dim(param)[1]> p[1]+(p[2]+p[3])*m ) C = as.matrix(t(param)[,p[1]+(p[2]+p[3])*m+1]) else C = NA
 }    else   {

          if (length(p)==5) {
            P = max(p[-1])
 		B = (1:(m*m*(P+1)*2))*0; dim(B) = c(m,m,(P+1),2)
 		A = B
 		AB = VECB[,1:p[1]]%*%t(beta)
 		B[,,1,1] = AB[,1:m]
 		A[,,1,1] = AB[,(m+1):(2*m)]
 		B[,,1,2] = AB[,1:m]
 		A[,,1,2] = AB[,(m+1):(2*m)]
 		for ( i in 2:(p[2]+1)) B[,,i,1] = VECB[, p[1]+((i-2)*m+1):((i-2)*m+m)]
 		for ( i in 2:(p[3]+1)) A[,,i,1] = VECB[, (p[1]+p[2]*m)+((i-2)*m+1):((i-2)*m+m)]
 		for ( i in 2:(p[4]+1)) B[,,i,2] = VECB[, (p[1]+(p[2]+p[3])*m)+((i-2)*m+1):((i-2)*m+m)]
 		for ( i in 2:(p[5]+1)) A[,,i,2] = VECB[, (p[1]+(p[2]+p[3]+p[4])*m)+((i-2)*m+1):((i-2)*m+m)] 

      	B[,,1:(p[2]+1),1] = CIB3B(B[,,1:(p[2]+1),1])
            B[,,1:(p[4]+1),2] = CIB3B(B[,,1:(p[4]+1),2])
      	A[,,1:(p[3]+1),1] = CIA2A(A[,,1:(p[3]+1),1])
            A[,,1:(p[5]+1),2] = CIA2A(A[,,1:(p[5]+1),2])
         	if ( dim(param)[1]> p[1]+(sum(p[-1]))*m)  {
                  C = (1:(m*2))*0; dim(C) = c(m,1,2)  
                  C[,1,1] = as.matrix(t(param)[,p[1]+sum(p[-1])*m+1]) 
			C[,1,2] = as.matrix(t(param)[,p[1]+sum(p[-1])*m+2])
			}  else {
                  C = NA
		}
          }
          if (length(p)==3)  {
             P = max(p[-1])
 		 B = (1:(m*m*(P+1)*2))*0; dim(B) = c(m,m,(P+1),2)
		 AB = VECB[,1:p[1]]%*%t(beta)
 		 B[,,1,1] = AB[,1:m]
   		 B[,,1,2] = AB[,1:m]
 		 for ( i in 2:(p[2]+1)) B[,,i,1] = VECB[, p[1]+((i-2)*m+1):((i-2)*m+m)]
 		 for ( i in 2:(p[3]+1)) B[,,i,2] = VECB[, (p[1]+p[2]*m)+((i-2)*m+1):((i-2)*m+m)]
      	 B[,,1:(p[2]+1),1] = CIB3B(B[,,1:(p[2]+1),1])
             B[,,1:(p[3]+1),2] = CIB3B(B[,,1:(p[3]+1),2])
		 A = NA
		 if ( dim(param)[1]> p[1]+(sum(p[-1]))*m)  {
                  C = (1:(m*2))*0; dim(C) = c(m,1,2)  
                  C[,1,1] = as.matrix(t(param)[,p[1]+sum(p[-1])*m+1]) 
			C[,1,2] = as.matrix(t(param)[,p[1]+sum(p[-1])*m+2])
			}  else {
                  C = NA
		}


          }
 }
 
 return(list(B,A,C))
}


#' ########
#' @export
B2CIB = function(B) {
    CIB = B*0
    L = dim(B)[3]
    n = dim(B)[1]
    for (i in L:1) for (j in L:i) CIB[,,i] =  CIB[,,i]-B[,,j]
    CIB[,,1] = -diag(n)- CIB[,,1]
    P = eigen(CIB[,,1])$vectors
    lambda = eigen(CIB[,,1])$values
    k = 0
    for (i in 1:(n-1)) {
      sm = sum(abs(CIB[,,1] - P[,i:n]%*%diag(lambda[i:n])%*%solve(P)[i:n,])) 
      if (sm < 0.000001) k = i
    }
    sm = sum(abs(CIB[,,1] - as.matrix(P[,n:n])%*%lambda[n:n]%*%t(as.matrix(solve(P)[n:n,])))) 
    if (sm < 0.000001) k = n
    if (k<n) {
       alpha =  P[,k:n]%*%diag(lambda[k:n])
       beta  =  t(solve(P)[k:n,]) 
       sm    =  sum(abs(CIB[,,1] - alpha%*%t(beta)))
    }
    if (k==n) {
       alpha = as.matrix(P[,n:n]*lambda[n:n])
       beta  = as.matrix(solve(P)[n:n,]) 
       sm    = sum(abs(CIB[,,1] - alpha%*%t(beta)))

    }
    return(list(CIB,alpha,beta,sm))
}


#' ########
#' @export
CIB3B = function(CIB) { 
	BB = CIB*0
	L = dim(CIB)[3]
      n = dim(BB)[1]
	if (L==1) BB = CIB + diag(n)
	if (L==2) {BB[,,1] = CIB[,,1]+CIB[,,2]+diag(n); BB[,,2] = - CIB[,,2]}
	if (L>2)  {
     		BB[,,1] = CIB[,,1]+CIB[,,2]+diag(n)
     		for (i in 2:(L-1))  BB[,,i] = CIB[,,i+1]-CIB[,,i]
     		BB[,,L] = - CIB[,,L]
	}
return(BB)
}

#' ########
#' @export
CIA2A = function(CIA) { 
	AA = CIA*0
	L = dim(CIA)[3]
      n = dim(AA)[1]
	if (L==1) AA = CIA 
	if (L==2) {AA[,,1] = CIA[,,1]+CIA[,,2]; AA[,,2] = - CIA[,,2]}
	if (L>2)  {
     		AA[,,1] = CIA[,,1]+CIA[,,2]
     		for (i in 2:(L-1))  AA[,,i] = CIA[,,i+1]-CIA[,,i]
     		AA[,,L] = - CIA[,,L]
	}
return(AA)
}

#' ########
#' @export
CIB2B = function(tst) {
    P  = tst$P
    r  = ncol(as.matrix(tst$beta))
    m  = ncol(tst[[2]]$coefficients)
    C  = matrix(0,m,1)
    C1 = C
    C2 = C
    mx = tst$NN1
    if (r==1) PI = as.matrix(tst[[2]]$coefficients[1:r,])%*%t(tst$beta)
    if (r>1 ) PI = t(tst[[2]]$coefficients[1:r,])%*%t(tst$beta)
    if (tst$model=="II"|tst$model=="IV") C1 = as.matrix(PI[,m+mx+1])
    PI = PI[,1:(m+mx)]
    BB = t(tst[[2]]$coefficients[(r+1):(r+(P[1]-1)*m),])
    BA = t(tst[[2]]$coefficients[(r+(P[1]-1)*m+1):(r+(P[1]-1)*m+(P[2]-1)*mx) ,])
    #BBAA  = t(tst[[2]]$coefficients[(r+1:(r+(P[1]-1)*m+(P[2]-1)*mx),])

    if (tst$model=="III"|tst$model=="IV") C2 = as.matrix(tst[[2]]$coefficients[nrow(tst[[2]]$coefficients),])
    C  = C1 + C2
    B  = (1:(m*m*P[1]))*0; dim(B) = c(m,m,P[1])
    A  = (1:(m*mx*P[2]))*0; dim(A) = c(m,mx,P[2])

    if ( P[1]==1 )  B[,,1] = PI[,1:m]+diag(m)
    if ( P[1]==2 ) {B[,,1] = PI[,1:m] + BB[,1:m] +diag(m); B[,,2] = - BB[,1:m]}	
    
    if ( P[1]>2 ) {
    	B[,,1] = PI[,1:m] + BB[,1:m] + diag(m)
      for (i in 2:(P[1]-1)) {
		B[,,i]    = - BB[,((i-2)*m+1):((i-2)*m+m)] + BB[,((i-1)*m+1):((i-1)*m+m)];	
      }
      B[,,P[1]] = - BB[,((P[1]-2)*m+1):((P[1]-2)*m+m)]
    }
    	

    if ( P[2]==1 )  A[,,1] = PI[,(m+1):(m+mx)]
    if ( P[2]==2 ) {A[,,1] = PI[,(m+1):(m+mx)] + BA[,1:mx]; A[,,2] = - BA[,1:mx]}	
    
    if ( P[2]>2 ) {
    	A[,,1] = PI[,(m+1):(m+mx)] + BA[,1:mx]
      for (i in 2:(P[2]-1)) {
		A[,,i]    = - BA[,((i-2)*mx+1):((i-2)*mx+mx)] + BA[,((i-1)*mx+1):((i-1)*mx+mx)];	
      }
      A[,,P[2]] = - BA[,((P[2]-2)*mx+1):((P[2]-2)*mx+mx)]
    }
 	
    return(list(B,A,C))
}




#' @export 
STAT = function(G) {
   L = dim(G)[3]
   m = dim(G)[1]
   PP=matrix(0,m*L,m*L) 
   PP[1:m,] = G[,,]
   if (L>1 ) {
     for ( i in 1:(L-1) ) PP[(i*m+1):(i*m+m),((i-1)*m+1):((i-1)*m+m)] = diag(m)
     eigen(PP)$values
   }  else {
     eigen(G[,,1])$values
   }
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




#' Generating impulse and response functions of an estimated CIVAR. 
#'
#' This function generates the impulse response functions of an estimated CIVAR model
#'
#' @param res	: output of CIVARest. 
#' @param nstep	: the length of the impulse response functions.
#' @param comb  : a weighting matrix specifying the weights used in the impulse response functions of a global VAR. Its default value is NA for CIVAR(p).  
#' @param irf   : types of the generated impusle response functions. 
#' @param nrun  : number of runs used in the the calculation of the boostrap confidence interval.
#' @param conf  : a two component vector containing the tail probabilities of the bootstrap confidence interval. 
#' @return       an array of dimesion (n, n, nstep, 3). 
#' @examples 
#'
#' res_d2 = CIVARData(n=4,p=2,T=84,Co=matrix(c(1,1,1,1),4,1)*0,type="none",crk=1)
#' res_e2 = CIVARest(res=res_d2)
#' res_e2$BIC
#' res_e2$AIC
#' res_e2$tst[[1]]
#' 
#' IRF = irf_B_sigma(B=res_e$B,sigma=res_e$Sigma,nstep=20,irf="gen1")
#' 
#' par(mfrow=c(4,4))
#' plot(IRF[1,1,],type="l")
#' plot(IRF[2,1,],type="l")
#' plot(IRF[3,1,],type="l")
#' plot(IRF[4,1,],type="l")
#' plot(IRF[1,2,],type="l")
#' plot(IRF[2,2,],type="l")
#' plot(IRF[3,2,],type="l")
#' plot(IRF[4,2,],type="l")
#' plot(IRF[1,3,],type="l")
#' plot(IRF[2,3,],type="l")
#' plot(IRF[3,3,],type="l")
#' plot(IRF[4,3,],type="l")
#' plot(IRF[1,4,],type="l")
#' plot(IRF[2,4,],type="l")
#' plot(IRF[3,4,],type="l")
#' plot(IRF[4,4,],type="l")
#'
#' 
#' IRF_CB = irf_CIVAR_CB(res=res_e2, nstep=20, comb=NA, irf = "gen1", runs = 20, conf = c(0.05, 0.95)) 
#' 
#' 
#' par(mfrow=c(4,4))
#' plott(IRF_CB,1,1) 
#' plott(IRF_CB,2,1) 
#' plott(IRF_CB,3,1) 
#' plott(IRF_CB,4,1) 
#' plott(IRF_CB,1,2) 
#' plott(IRF_CB,2,2) 
#' plott(IRF_CB,3,2) 
#' plott(IRF_CB,4,2) 
#' plott(IRF_CB,1,3) 
#' plott(IRF_CB,2,3) 
#' plott(IRF_CB,3,3) 
#' plott(IRF_CB,4,3) 
#' plott(IRF_CB,1,4) 
#' plott(IRF_CB,2,4) 
#' plott(IRF_CB,3,4) 
#' plott(IRF_CB,4,4) 
#'
#' @export
irf_CIVAR_CB = function (res, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95)) 
{
    n = res$n
    p = res$p
    T = dim(res$Y)[1]
    A = res$A
    B = res$B
    Co = res$Co
    type  = res$type
    X     = res$X
    crk   = res$crk
    neq   = dim(B)[1]
    nvar  = dim(B)[2]
    sigma = res$Sigma
    response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
    response[, , , 1] <- irf_B_sigma(B=B,sigma=sigma,nstep,comb,irf)
    responseR <- array(0, dim = c(neq, nvar, nstep, runs))
    for (i in 1:runs) {
        Uo_run    = rnormSIGMA(T, sigma)
        res_run   = CIVARData(n, p, T, r_np = NA, A=NA, B, Co, U = Uo_run, Sigma = NA, type,X=NA,mu=NA,Yo=NA,crk)
        res_e     = CIVARest(res_run)
        B_run     = res_e$B
        sigma_run = res_e$Sigma
        responseR[, , , i] <- irf_B_sigma(B=res_e$B, sigma=res_e$Sigma, nstep, comb, irf)
    }
    for (tt in 1:(nstep)) {
        for (i in 1:neq) {
            for (j in 1:nvar) {
               
                response[i, j, tt, -1] = quantile(responseR[i,j, tt, ], conf)
            }
        }
    }
    return(response)
}