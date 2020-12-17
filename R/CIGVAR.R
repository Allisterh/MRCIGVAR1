#' Data generating process of cointegrated global VAR CIGVAR(m,n,p,T,...) 
#'
#' This function will generate data from a cointegrated GVAR(p) process and return a list containing data and parameters used in the CIGVAR(p) process.
#'
#' @param m     : number of variables in each a country/unit
#' @param n     : number of countries/units
#' @param p     : (n x 3) matrix, each raw contains the lag length of the domestic variables, the foreign variables and the number of the exogeneous variables
#' @param T     : number of observations. 
#' @param               (m,n,p,T) are parameters which must be provided. 
#' @param W     : (n x n) weighting matrix. w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' @param r_npo : (m, p, n) array collecting the roots of the characteristic functions in the lag operator of the country VAR. The number of ones in r_npo[,,i] is the number of unit roots i-th country/unit. 
#' @param Ao    : (m, m, p, n) array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
#' @param Bo    : (m, m, p, n) array collecting the n country VAR(p) coefficients. 
#' @param Co    : (m , k+1, n) array collecting the coeeficients of the deterministic components of the n countries.
#' @param Uo    : an (T x mn) matrix of the temporally independent innovation processes
#' @param Sigmao	    : (mn x mn) matrix of the covariance matrix of the CIGVAR(m,n,p)
#' @param 		      (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly. The default assumption is one unit root in one country. Hence m-1 cointegration relations in each country.
#' @param type	: deterministic component "const" and "none" are two options 
#' @param X	: (T x k) matrix of exogeneous variables.
#' @param mu    : if type = "const" mu has the same dimension as Co. is an muV is nm vector of the means of the time series in the system
#' @param d	: d = 0 implies foreign variables are not in the cointegration space. d = 1 allows foreign variables enter the cointegration relation.
#' @param crk   : n vector containing the cointegration rank in each country/unit. crk is to specified for estimation of parameters. 
#' @return      A CIGVAR object containing the generated data, the parameters and the input exogeous variables. res = list("Y","Uo","G","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
#' @field Y     : T x nm simulated data via of the GVAR(m,n,p). 
#' @field X     : (T x k) matrix of exogeneous variables.
#' @field Uo    : (T x mn) matrix of the simulated innovations of the GVAR(m,n,p) 
#' @field G     : (mn, mn, p)  array of the GVAR(m,n,p) coefficients. G is contructed from Bo, Ao and W.
#' @field C     : (nm x (k+1)) matrix containing the coefficients of the deterministic components.   
#' @field Sigmao: (mn x mn) matrix of the covariance matrix of the CIGVAR(m,n,p)
#' @field r_npo : (m, p, n) array collecting the roots of the characteristic functions in L of the n dynamically independent domestic VAR(p)s.
#' @examples 
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3; 
#' res_d = CIGVARData(m=4,n=5,p=p,T=1000,type="const")
#' max(abs(res_d$Y))
#' plot(ts(res_d$Y[,1:10]))
#' STAT(res_d$G)
#' res_e = CIGVARest(res_d)
#' res_e$CRKtst
#' summary_GVAR(res_e)
#' res_d$crk = rep(3,5)
#' res_e = CIGVARest(res_d)
#' summary_GVAR(res_e)
#' STAT(res_d$G)
#' plot(ts(res_d$Y[,1:10]))
#' @export
CIGVARData = function (m, n, p, T, W = NA, r_npo = NA, Ao = NA, Bo = NA, Co = NA, 
    Uo = NA, Sigmao = NA, type = NA, X = NA, mu = NA, d = 0,crk=NA) 
{
    if (missing(Bo)) {
        Bo = NA
    }
    if (missing(Sigmao)) {
        Sigmao = NA
    }
    if (missing(Uo)) {
        Uo = NA
    }
    if (missing(W)) {
        W = NA
    }
    if (missing(Ao)) {
        Ao = NA
    }
    if (missing(type)) {
        type = NA
    }
    if (missing(Co)) {
        Co = NA
    }
    if (missing(mu)) {
        mu = NA
    }
    if (missing(d)) {
        d = NA
    }
    if (missing(X)) {
       X=NA
    }

    if (missing(r_npo)) {
        #return("missing r_npo")
        maxPd = max(p[,1])
        r_npo = c(1:(m*maxPd*n))/c(1:(m*maxPd*n))*1.3
        dim(r_npo) = c(m,maxPd,n)
        r_npo[1,1,] = 1  
    }

    if (missing(crk)) {
       crk = (1:n)/(1:n)
    }


    alpha = list()
    beta  = list()

    if (anyNA(d)) 
        d = 1
    Pmax = max(p[, 1:2])
    P = max(p)
    if (!anyNA(X)) 
        k = dim(X)[2]
    if (anyNA(Bo)) {
        Bo = (1:(m * m * Pmax * n)) * 0
        dim(Bo) = c(m, m, Pmax, n)
        for (i in 1:n) {
                r_np = c(1:(m * p[i,1])) * 0
                dim(r_np) = c(m,p[i,1])
                r_np = r_npo[,1:p[i,1],i]
            VARD = VARData(m, p[i, 1], T,r_np=r_np)
            Bo[, , 1:p[i, 1], i] = VARD$B
            #r_npo[, 1:p[i, 1], i] = VARD$r_np
            alpha[[i]] = B2CIB(VARD$B)[[2]]
            beta[[i]]  = B2CIB(VARD$B)[[3]]
        }
    }

    
    if (anyNA(Sigmao)) {
        Sigmao = matrix(0, m * n, m * n)
        VARD = VARData(m * n, p[1, 1], T)
        Sigmao = VARD$Sigma
    }
    if (anyNA(Uo)) {
        Uo = rnormSIGMA(T, Sigmao)
    }
    if (anyNA(W)) {
        W = matrix(((1:(n * n))/(1:(n * n)) * 1/(n - 1)), n, 
            n)
        for (i in 1:n) W[i, i] = 0
    }
    if (anyNA(Ao)) {
        Ao = (1:(m * m * Pmax * n)) * 0
        dim(Ao) = c(m, m, Pmax, n)
        for (i in 1:n) {
        if (p[i,2] < 2) Ao = Ao  
        if (p[i,2] >= 2) {
                VARD = VARData(m, p=(p[i, 2]-1),T)
            BB = VARD$B/1
            Ao[,,1,i] = BB[,,1]+alpha[[i]]%*%t(beta[[i]])*d
            Ao[,,p[i,2],i] = -BB[,,p[i,2]-1]
            if ((p[i, 2]-1)>=2) { for (L in 2:(p[i, 2]-1)) Ao[, ,L, i] = BB[,,L]-BB[,,L-1]}
        }
        }
    }
    if (anyNA(type)) {
        type = "none"
    }
    if (type == "none") {
        Co = matrix(0, m, n)
        dim(Co) = c(m, 1, n)
        mu = matrix(0, m, n)
        dim(mu) = c(m, 1, n)
    }
    G = (1:(n * m * n * m * Pmax)) * 0
    dim(G) = c(n * m, n * m, Pmax)
    for (i in 1:n) {
        for (j in 1:n) {
            for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + 
                (j - 1) * m):(j * m), L] = Ao[, , L, i] * W[i, 
                j]
        }
    }
    for (i in 1:n) {
        for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + (i - 
            1) * m):(i * m), L] = Bo[, , L, i]
    }
    Ct = Uo * 0
    if (type == "const") {
        if (anyNA(mu)) {
            mu = matrix(rnorm(n * m), m, n)
            dim(mu) = c(m, 1, n)
        }
        if (anyNA(Co)) {
            Co = mu
            muV = as.vector(mu)
            CoV = muV
            for (L in 1:Pmax) CoV = CoV - G[, , L] %*% muV
            Co = CoV; dim(Co) = c(m, 1, n)
        }
        else {
            H = diag(n * m)
            for (L in 1:Pmax) H = H - G[, , L]
            CoV = as.vector(Co)
            if ( min(abs(eigen(H)$values)) > 0.001 )  {
                muV = solve(H) %*% CoV
                mu = muV
                dim(mu) = c(m, 1, n)
              }   else   {
                mu = NA

            }

              
        }
        Ct = matrix(1, T, 1) %*% t(CoV)
    }
    if (type == "exog0") {
        if (anyNA(Co)) {
            Co = matrix(rnorm(m * n * (k + 1)), m * (k + 1), 
                n)
            dim(Co) = c(m, k + 1, n)
            Co[, 1, ] = (1:m) * 0
            for (i in 1:n) if (p[i, 3] < k) 
                Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i, 3] + 
                  2):(k + 1), i] * 0
        }
        DMCo = dim(Co[, -1, ])
        if (length(DMCo) < 3) 
            DMCo = c(dim(Co)[1], 1, dim(Co)[3])
        if (!(DMCo[2] == dim(X)[2]) | !(DMCo[1] == m) | !(DMCo[3] == 
            n)) {
            print("dimension problem")
            return("dimension")
        }
        CoV = matrix(0, dim(X)[2], m * n)
        for (i in 1:dim(X)[2]) CoV[i, ] = as.vector(Co[, 1 + 
            i, ])
        Ct = matrix(0, dim(X)[1], m * n)
        for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * m + 
            m)] = as.matrix(X[, , i]) %*% CoV[, ((i - 1) * m + 
            1):((i - 1) * m + m)]
        mu = NA
    }
    else {
        if (type == "exog1") {
            if (anyNA(Co)) {
                Co = matrix(rnorm(m * n * (k + 1)), m * (k + 
                  1), n)
                dim(Co) = c(m, k + 1, n)
                for (i in 1:n) if (p[i, 3] < k) 
                  Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i, 
                    3] + 2):(k + 1), i] * 0
            }
            DMCo = dim(Co[, -1, ])
            if (!(DMCo[2] == dim(X)[2]) | !(DMCo[1] == m) | !(DMCo[3] == 
                n)) {
                print("dimension problem")
                return("dimension")
            }
            CoV = matrix(0, dim(X)[2] + 1, m * n)
            for (i in 1:(1 + dim(X)[2])) CoV[i, ] = as.vector(Co[, 
                i, ])
            Ct = matrix(0, dim(X)[1], m * n)
            for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * 
                m + m)] = cbind(matrix(1, dim(X)[1], 1), X[, 
                , i]) %*% CoV[, ((i - 1) * m + 1):((i - 1) * 
                m + m)]
            mu = NA
        }
        else {
            X = NA
        }
    }
    Y = Uo + Ct
    for (t in (P + 1):T) {
        for (L in 1:Pmax) Y[t, ] = Y[t, ] + Y[t - L, ] %*% t(G[, 
            , L])
    }
    check = max(abs(Y))
    result = list(Y, X, Uo, G, C, Sigmao, r_npo, Ao, Bo, Co, 
        W, m, n, p, mu, check, type,crk)
    names(result) = c("Y", "X", "Uo", "G", "C", "Sigmao", "r_npo", 
        "Ao", "Bo", "Co", "W", "m", "n", "p", "mu", "check", 
        "type","crk")
    return(result)
}

#' Estimation of CIGVAR(m,n,p) 
#'
#' This function estimate the unknown parameters of a CIGVAR(m,n,p) object based on provided data.
#' It runs a VECMX estimation country/unit by country/unit under given crk and pieces the results together to obtain a CIGVAR. 
#' @param  res  :aCIGVAR object that can be a output of CIGVARData or CIGVARData1 including at least values of m,n,p,type,Y,crk and optionally X. 
#' @return res  :a list as the input, but filled with estimated parameter values, AIC, BIC and LH
#' @examples 
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3; 
#' r_npo = (1:(3 * 3 * 5))*0; dim(r_npo) = c(3,3,5)
#' r_npo[,,1] = matrix(c(1,1,3,2,3,3,3,3,3),3,3)
#' r_npo[,,2] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,3] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,4] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,5] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' 
#' res_d = CIGVARData1(m=3,n=5,p=p,T=500,r_npo=r_npo,type="const",d=0,Ncommtrend=1)
#' max(res_d$Y)
#' STAT(res_d$G)
#' plot(ts(res_d$Y[,1:10]))
#' res_e = CIGVARest(res_d)
#' summary_GVAR(res_e)
#' res_e$CRKtst
#' res_d$crk = c(1,2,2,2,2)
#' res_e = CIGVARest(res_d)
#' summary_GVAR(res_e)
#' #### unit roots on the GVAR
#' STAT(res_e$G)
#' @export
CIGVARest=function (res) 
{
    m = res$m
    n = res$n
    p = res$p
    Y = res$Y
    X = res$X
    W = res$W
    crk = res$crk
    type = res$type
    Bo = res$Bo
    Ao = res$Ao
    Co = res$Co
    Pmax = max(p[, 1:2])
    k = max(p[, 3])
    AIC = c(1:n) * 0
    BIC = c(1:n) * 0
    LH = c(1:n) * 0
    resid = Y * 0
    Tresid = resid[,1:m]
    VAR_domestic = list()
    CRKtst       = list()

    T = dim(Y)[1]
    FY = Y %*% t(W %x% diag(m))
    Bo = Bo = (1:(m * m * Pmax * n)) * 0
    dim(Bo) = c(m, m, Pmax, n)
    Ao = Bo * 0
    if (type == "none") 
        C = matrix(0, m * n, 1) * 0
    if (type == "const") 
        C = matrix(0, m * n, 1) * 0
    if (type == "exog0") 
        C = matrix(0, m * n, dim(X)[2] + 1) * 0
    if (type == "exog1") 
        C = matrix(0, m * n, dim(X)[2] + 1) * 0
    for (i in 1:n) {
        y = Y[,((i-1)*m+1):((i-1)*m+m)]
        x = FY[,((i-1)*m+1):((i-1)*m+m)]
        if (type=="none")  Model ="I"
        if (type=="const") Model ="III"
        P = matrix(0,2,3); P[1,1:2]=p[i,1:2];P[2,1:2]=p[i,1:2]
        #tst<-cvecmgvar(y,x,model = Model, type = "eigen",P = p[i,1:2], r = crk[i], q = 0.95)
        tst<-MRCVECMest2(y,x,model = Model, type = "eigen",P = P , crk = crk[i], q = 0.95,Dxflag=0)

        BAC = VECM2VAR(param=tst[[2]][[1]],beta=tst$beta,p=c(crk[i],P[1,1]-1,P[1,2]-1))
        #BAC = CIB2B(tst)  
        
        Tresid[(T-dim(tst[[2]]$residuals)[1]+1):T,]=tst[[2]]$residuals
        Sigma_one = t(Tresid)%*%Tresid/(T-dim(tst[[2]][[1]])[1])
    
        LH_P     = -(T*m/2)*log(2*pi) -(T*m/2) + (T/2)*log(det(solve(Sigma_one)))
        LH_AIC   =           2*(n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2) - 2*LH_P
        LH_BIC   =      log(T)*(n*(dim(tst[[2]][[1]])[1])+n*(n+1)/2) - 2*LH_P
        LH_N     =     2*n*(dim(tst[[2]][[1]])[1])+n*(n+1)
     
        AIC[i]   = LH_AIC
        BIC[i]   = LH_BIC
        LH[i]    = LH_P
        VAR_domestic[[i]] = tst[[2]]
        CRKtst[[i]]       = tst[[1]]

        Bo[, , 1:p[i, 1], i] = BAC[[1]]
        Ao[, , 1:p[i, 2], i] = BAC[[2]]
        Co[, , i]            = BAC[[3]]
        resid[1:nrow(tst[[2]]$residuals), (m * (i - 1) + 1):(i * m)] = tst[[2]]$residuals
    }
    Sigmao = t(resid) %*% resid/(T - m * (p[i, 1] + p[i, 2]) -  p[i, 3])
    G = BoAoW2G(Bo, Ao, W, m, n, Pmax)
    Gs = diag(n * m)
    for (L in 1:Pmax) {
        Gs = Gs - G[, , L]
    }
    LH_g = -T * m * n/2 * log(2 * pi) - T * m * n/2 + T/2 * log(det(solve(Sigmao)))
    TwoN =  2*(sum(p[,1:2])*m*m+sum(p[,3])*m + (type=="const")*n*m) 
    AIC_g = TwoN - 2 * LH_g
    BIC_g = log(T) * TwoN/2 - 2 * LH_g
    res$G = G
    res$C = C
    res$Sigmao = Sigmao
    res$r_npo = NA
    res$Ao = Ao
    res$Bo = Bo
    res$Co = Co
    #res$mu = solve(Gs) %*% C
    res$VAR_domestic = VAR_domestic
    res$AIC = AIC
    res$BIC = BIC
    res$LH = LH
    res$LH_g = LH_g
    res$AIC_g = AIC_g
    res$BIC_g = BIC_g
    res$CRKtst = CRKtst
    return(res)
}


#' Data generating process of a cointegrated global VAR CIGVAR(n,m,p,T,...) 
#'
#' This function will generate CIGVAR process with common and idiosycratic stochastic trends.
#'
#' @param m     : number of variables in each a country/unit
#' @param n     : number of countries/units
#' @param p     : (n x 3) matrix, each raw contains the lag length of the domestic variables, the foreign variables and the number of the exogeneous variables
#' @param T     : number of observations. 
#'                (m,n,p,T) are parameters which must be provided. 
#' @param W     : n x n weighting matrix. w_ij is the weight of foreign country j in the foreign variables of ith country diag(W)=0
#' @param r_npo : m x p x n array collecting the roots of the characteristic functions of the country VAR. The number of ones in r_npo[,,i] is the number of unit roots i-th country/unit. 
#' @param Ao    : m x m x p x n array collecting the off-diagonal block of coefficents which represent the inter-country lag coefficients (coefficients of foreign variables)
#' @param  Bo   : m x m x p x n array collecting the n country VAR(p) coefficients.  Bo are coefficents of stationary domestic VAR(p). 
#' @param  Co   : m x (k+1) x n array collecting the coeeficients of the deterministic components of the n countries.
#' @param  Uo   : an T x mn matrix of the temporally independent innovation processes
#' @param  Sigmao	    : mn x mn matrix of the covariance matrix of the GVAR(m,n,p)
#'  		      (W,r_npo,Ao,Bo,Uo,Sigmao) if not provided, they will be generated randomly.
#' @param  type	: deterministic component "const" and "none" are two options 
#' @param  X	: (T x k) matrix of exogeneous variables.
#' @param  mu   : if type = "const" mu has the same dimension as Co. is an muV is nm vector of the means of the time series in the system
#' @param  d	: d = 0 implies foreign variables are not in the cointegration space. d = 1 allows foreign variables enter the cointegration relation.
#' @param  crk  : n vector containing the cointegration rank in each country/unit. crk is used only in estimation.   
#' @param  Ncommtrend : number of common stochastic trends across all countries. 
#' @return      An CIGVAR object containing the generated data, the parameters and the input exogeous variables. res = list("Y","Uo","G","C","Sigmao","r_npo","Ao","Bo","W","m","n","p","check","type","mu")
#' @field Y     : T x nm simulated data via of the GVAR(m,n,p). 
#' @field X     : (T x k) matrix of exogeneous variables.
#' @field Uo    : T x mn array of the simulated innovations of the GVAR(m,n,p) 
#' @field G     : mn x mn x p  array of the GVAR(m,n,p) coefficients. G is contructed from Bo, Ao and W.
#' @field C     : nm x (k+1) matrix containing the coefficients of the deterministic components.   
#' @field Sigmao: mn x mn matrix of the covariance matrix of the CIGVAR(m,n,p)
#' @field r_npo : m x p x n matrix collecting the roots of the characteristic functions in L of the n dynamically independent domestic VAR(p)s.
#' @examples 
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3; 
#' r_npo = (1:(3 * 3 * 5))*0; dim(r_npo) = c(3,3,5)
#' r_npo[,,1] = matrix(c(1,1,3,2,3,3,3,3,3),3,3)
#' r_npo[,,2] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,3] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,4] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' r_npo[,,5] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
#' 
#' res_d = CIGVARData1(m=3,n=5,p=p,T=500,r_npo=r_npo,type="const",d=0,Ncommtrend=1)
#' max(res_d$Y)
#' STAT(res_d$G)
#' plot(ts(res_d$Y[,1:10]))
#' res_e = CIGVARest(res_d)
#' summary_GVAR(res_e)
#' res_e$CRKtst
#' res_d$crk = c(1,2,2,2,2)
#' res_e = CIGVARest(res_d)
#' summary_GVAR(res_e)
#' #### unit roots on the GVAR
#' STAT(res_e$G)
#' @export
CIGVARData1 = function (m, n, p, T, W = NA, r_npo = NA, Ao = NA, Bo = NA, Co = NA, 
    Uo = NA, Sigmao = NA, type = NA, X = NA, mu = NA, d = 0,crk=NA,Ncommtrend=1) 
{
    if (missing(Bo)) {
        Bo = NA
    }
    if (missing(Sigmao)) {
        Sigmao = NA
    }
    if (missing(Uo)) {
        Uo = NA
    }
    if (missing(W)) {
        W = NA
    }
    if (missing(Ao)) {
        Ao = NA
    }
    if (missing(type)) {
        type = NA
    }
    if (missing(Co)) {
        Co = NA
    }
    if (missing(mu)) {
        mu = NA
    }
    if (missing(d)) {
        d = NA
    }
    if (missing(X)) {
       X=NA
    }

    if (missing(r_npo)) {
        #return("missing r_npo")
        maxPd = max(p[,1])
        r_npo = c(1:(m*maxPd*n))/c(1:(m*maxPd*n))*1.3
        dim(r_npo) = c(m,maxPd,n)
        r_npo[1,1,] = 1  
    }

    if (missing(crk)) {
       crk = (1:n)/(1:n)
    }


    alpha = list()
    beta  = list()

    if (anyNA(d)) 
        d = 1
    Pmax = max(p[, 1:2])
    P = max(p)
    if (!anyNA(X)) 
        k = dim(X)[2]
    if (anyNA(Bo)) {
        Bo = (1:(m * m * Pmax * n)) * 0
        dim(Bo) = c(m, m, Pmax, n)
        #for (i in 1:n) {
        #        r_np = c(1:(m * p[i,1])) * 0
        #        dim(r_np) = c(m,p[i,1])
        #        r_np = r_npo[,1:p[i,1],i]
        #    VARD = VARData(m, p[i, 1], T,r_np=r_np)
        #    Bo[, , 1:p[i, 1], i] = VARD$B
        #    #r_npo[, 1:p[i, 1], i] = VARD$r_np
        #    alpha[[i]] = B2CIB(VARD$B)[[2]]
        #    beta[[i]]  = B2CIB(VARD$B)[[3]]
        #}
    

        Balphabeta = VARB_commtrend(m,p,T,r_npo,Ncommtrend,n)
        Bo    = Balphabeta[[1]]
        alpha = Balphabeta[[2]]
        beta  = Balphabeta[[3]]
    }


    
    if (anyNA(Sigmao)) {
        Sigmao = matrix(0, m * n, m * n)
        VARD = VARData(m * n, p[1, 1], T)
        Sigmao = VARD$Sigma
    }
    if (anyNA(Uo)) {
        Uo = rnormSIGMA(T, Sigmao)
    }
    if (anyNA(W)) {
        W = matrix(((1:(n * n))/(1:(n * n)) * 1/(n - 1)), n, 
            n)
        for (i in 1:n) W[i, i] = 0
    }
    if (anyNA(Ao)) {
        Ao = (1:(m * m * Pmax * n)) * 0
        dim(Ao) = c(m, m, Pmax, n)
        for (i in 1:n) {
        if (p[i,2] < 2) Ao = Ao  
        if (p[i,2] >= 2) {
                VARD = VARData(m, p=(p[i, 2]-1),T)
            BB = VARD$B/1
            Ao[,,1,i] = BB[,,1]+alpha[[i]]%*%t(beta[[i]])*d
            Ao[,,p[i,2],i] = -BB[,,p[i,2]-1]
            if ((p[i, 2]-1)>=2) { for (L in 2:(p[i, 2]-1)) Ao[, ,L, i] = BB[,,L]-BB[,,L-1]}
        }
        }
    }
    if (anyNA(type)) {
        type = "none"
    }
    if (type == "none") {
        Co = matrix(0, m, n)
        dim(Co) = c(m, 1, n)
        mu = matrix(0, m, n)
        dim(mu) = c(m, 1, n)
    }
    G = (1:(n * m * n * m * Pmax)) * 0
    dim(G) = c(n * m, n * m, Pmax)
    for (i in 1:n) {
        for (j in 1:n) {
            for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + 
                (j - 1) * m):(j * m), L] = Ao[, , L, i] * W[i, 
                j]
        }
    }
    for (i in 1:n) {
        for (L in 1:Pmax) G[(1 + (i - 1) * m):(i * m), (1 + (i - 
            1) * m):(i * m), L] = Bo[, , L, i]
    }
    Ct = Uo * 0
    if (type == "const") {
        if (anyNA(mu)) {
            mu = matrix(rnorm(n * m), m, n)
            dim(mu) = c(m, 1, n)
        }
        if (anyNA(Co)) {
            Co = mu
            muV = as.vector(mu)
            CoV = muV
            for (L in 1:Pmax) CoV = CoV - G[, , L] %*% muV
            Co = CoV; dim(Co) = c(m, 1, n)
        }
        else {
            H = diag(n * m)
            for (L in 1:Pmax) H = H - G[, , L]
            CoV = as.vector(Co)
            if ( min(abs(eigen(H)$values)) > 0.001 )  {
                muV = solve(H) %*% CoV
                mu = muV
                dim(mu) = c(m, 1, n)
              }   else   {
                mu = NA

            }

              
        }
        Ct = matrix(1, T, 1) %*% t(CoV)
    }
    if (type == "exog0") {
        if (anyNA(Co)) {
            Co = matrix(rnorm(m * n * (k + 1)), m * (k + 1), 
                n)
            dim(Co) = c(m, k + 1, n)
            Co[, 1, ] = (1:m) * 0
            for (i in 1:n) if (p[i, 3] < k) 
                Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i, 3] + 
                  2):(k + 1), i] * 0
        }
        DMCo = dim(Co[, -1, ])
        if (length(DMCo) < 3) 
            DMCo = c(dim(Co)[1], 1, dim(Co)[3])
        if (!(DMCo[2] == dim(X)[2]) | !(DMCo[1] == m) | !(DMCo[3] == 
            n)) {
            print("dimension problem")
            return("dimension")
        }
        CoV = matrix(0, dim(X)[2], m * n)
        for (i in 1:dim(X)[2]) CoV[i, ] = as.vector(Co[, 1 + 
            i, ])
        Ct = matrix(0, dim(X)[1], m * n)
        for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * m + 
            m)] = as.matrix(X[, , i]) %*% CoV[, ((i - 1) * m + 
            1):((i - 1) * m + m)]
        mu = NA
    }
    else {
        if (type == "exog1") {
            if (anyNA(Co)) {
                Co = matrix(rnorm(m * n * (k + 1)), m * (k + 
                  1), n)
                dim(Co) = c(m, k + 1, n)
                for (i in 1:n) if (p[i, 3] < k) 
                  Co[, (p[i, 3] + 2):(k + 1), i] = Co[, (p[i, 
                    3] + 2):(k + 1), i] * 0
            }
            DMCo = dim(Co[, -1, ])
            if (!(DMCo[2] == dim(X)[2]) | !(DMCo[1] == m) | !(DMCo[3] == 
                n)) {
                print("dimension problem")
                return("dimension")
            }
            CoV = matrix(0, dim(X)[2] + 1, m * n)
            for (i in 1:(1 + dim(X)[2])) CoV[i, ] = as.vector(Co[, 
                i, ])
            Ct = matrix(0, dim(X)[1], m * n)
            for (i in 1:n) Ct[, ((i - 1) * m + 1):((i - 1) * 
                m + m)] = cbind(matrix(1, dim(X)[1], 1), X[, 
                , i]) %*% CoV[, ((i - 1) * m + 1):((i - 1) * 
                m + m)]
            mu = NA
        }
        else {
            X = NA
        }
    }
    Y = Uo + Ct
    for (t in (P + 1):T) {
        for (L in 1:Pmax) Y[t, ] = Y[t, ] + Y[t - L, ] %*% t(G[, 
            , L])
    }
    check = max(abs(Y))
    result = list(Y, X, Uo, G, C, Sigmao, r_npo, Ao, Bo, Co, 
        W, m, n, p, mu, check, type,crk)
    names(result) = c("Y", "X", "Uo", "G", "C", "Sigmao", "r_npo", 
        "Ao", "Bo", "Co", "W", "m", "n", "p", "mu", "check", 
        "type","crk")
    return(result)
}



################

 
#VARB_commtrend(m,p,T,r_npo,Ncommtrend=1,n)


#' @export
VARB_commtrend = function(m,p,T,r_npo, Ncommtrend,n) {  
   alpha = list()
   beta  = list()

  if (missing(Ncommtrend))  {Ncommtrend=NA}
  Pmax = max(p[,1:2])
  Pmin = min(p[,1:2][p[,1:2]>0])
  Bo = (1:(m * m * Pmax * n)) * 0
  dim(Bo) = c(m, m, Pmax, n)
  
  if (anyNA(Ncommtrend)) {
      for (i in 1:n) {
                r_np = c(1:(m * p[i,1]))*0   
                dim(r_np) = c(m,p[i,1])
                r_np = r_npo[,1:p[i,1],i]
            VARD = VARData(m, p[i, 1], T,r_np=r_np)
            Bo[, , 1:p[i, 1], i] = VARD$B
            #r_npo[, 1:p[i, 1], i] = VARD$r_np
            alpha[[i]] = B2CIB(VARD$B)[[2]]
            beta[[i]]  = B2CIB(VARD$B)[[3]]
      }
  }
  if (!anyNA(Ncommtrend)) {
      #Pmin = min(p[,1:2])   # put in the top
      ## B transform the block diagonal underlying coefficients
      B =  matrix(rnorm(m*m),m,m)                 ### jizhu 0 nadiao
      Nnocommtrend = m - Ncommtrend
      B[(Nnocommtrend+1):m,] = 0 
      B = B + diag(m)
     	r_np = c(1:(Ncommtrend*Pmin))/c(1:(Ncommtrend*Pmin))*1.5
      dim(r_np) = c(Ncommtrend,Pmin)
      r_np[1:Ncommtrend,1]  = 1     
      VAR_IONE  = VARData(n=Ncommtrend,p=Pmin,T=100,r_np=r_np)

     
      for (i in 1:n) {
         	B =  matrix(rnorm(m*m),m,m)                 ### jizhu 0 nadiao
      	B[(Nnocommtrend+1):m,] = 0 
      	B = B + diag(m)

 
            Nnocommtrend = m - Ncommtrend
      	r_np = c(1:(Nnocommtrend * p[i,1]))*0      
            dim(r_np) = c(Nnocommtrend,p[i,1])
            r_np <- r_npo[(Ncommtrend+1):m,1:p[i,1],i]
            dim(r_np) = c(Nnocommtrend,p[i,1])

            VAR_IZero = VARData(Nnocommtrend, p[i, 1],T=100,r_np=r_np)
            for (j in 1: p[i,1]) {
               Dblock    = matrix(0,m,m)
               Dblock[1:Nnocommtrend,1:Nnocommtrend]     = VAR_IZero$B[,,j]
               if (dim(VAR_IONE$B)[3]>=j) Dblock[(Nnocommtrend+1):m,(Nnocommtrend+1):m] = VAR_IONE$B[,,j]
               Bo[, , j, i] = B%*%Dblock%*%solve(B)
      	}
            #alpha[[i]] = B2CIB(VARD$B)[[2]]
            #beta[[i]]  = B2CIB(VARD$B)[[3]]
            alpha[[i]] = B2CIB(Bo[, , 1:p[i, 1], i])[[2]]
            beta[[i]]  = B2CIB(Bo[, , 1:p[i, 1], i])[[3]]
      }
  }
  return(list(Bo,alpha,beta))
}



























#' Transform coefficient matrices in of CIGVAR to the coefficients matrix of error correction form.  
#' 
#' @param  B  : array of dimension (m,m,L) containing the lag coefficients of the level VAR.
#' @return      A list containing three components.
#' @field CIB        : the coefficients matrix of the error correction form
#' @field alpha	  : the adjustment coefficients 
#' @field beta       : the cointegration vectors
#' @field sm         : the discrepance between the selected alpha beta and CIB[,,1].
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

#' Transform coefficient matrices of the error correction form to the coefficients matrix of CIGVAR in level.  
#' 
#' @param  tst  : an output of CIGVARest
#' @return      A list containing three components.
#' B            : the coefficient matrices of the domestic variables
#' A	        : the coefficient matrices of the foreign variables
#' C            : the coefficient matrices of the deterministic components. 
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


#' Transfor the coefficients in VEC form to the coefficients of the VAR in level 
#' @param CIB of dimension (m,m,L) of the coefficients matrices of a VECM with lag Dy_L-1.
#' @return B of dimenstion (m,mL) of the coefficients of a VAL with lag L
#' @export
CIB3B = function(CIB) { 
	BB = CIB*0
	L = dim(CIB)[3]
      n = dim(BB)[1]
	if (L==1) BB = CIB + diag(4)
	if (L==2) {BB[,,1] = CIB[,,1]+CIB[,,2]+diag(n); BB[,,2] = - CIB[,,2]}
	if (L>2)  {
     		BB[,,1] = CIB[,,1]+CIB[,,2]+diag(n)
     		for (i in 2:(L-1))  BB[,,i] = CIB[,,i+1]-CIB[,,i]
     		BB[,,L] = - CIB[,,L]
	}
return(BB)
}

### not deleted



#' Calculation of the model selection values of for CIGVAR models  
#' 
#' @param  res  : an object generated from CIGVARData or estimated from CIGVARest.
#' @param  L_V  : a 2 components vector containing the maxima of the domestic lag and the foreign lag respectively. 
#' @return      : A matrix with different lag specifications and the model selection criteria.
#' @examples 
#' 
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3; 
#' res_d = GVARData(m=2,n=5,p=p,T=1000,type="const")
#' r_npo = (1:(2 * 3 * 5))*0; dim(r_npo) = c(2,3,5)
#' r_npo[,,1] = matrix(c(1,1.1,1.3,1.2,1.3,1.3),2,3)
#' r_npo[,,2] = matrix(c(1,1.1,1.3,1.2,1.3,1.3),2,3)
#' r_npo[,,3] = matrix(c(1,1.1,1.3,1.2,1.3,1.3),2,3)
#' r_npo[,,4] = matrix(c(1,1.1,1.3,1.2,1.3,1.3),2,3)
#' r_npo[,,5] = matrix(c(1,1.1,1.3,1.2,1.3,1.3),2,3)
#' 
#' res_d = CIGVARData(m=2,n=5,p=p,T=500,r_npo=r_npo,type="const",d=0)
#' res_e = CIGVARest(res_d)
#' I = 1
#'
#' L_V = c(3,3)
#' 
#' res = res_d
#' res=res_e
#' 
#' CIGVARSelecte = CIGVAR_Select(res=res_e,L_V=c(4,4))
#'
#' CIGVARSelecte[which.min(CIGVARSelecte[,17]),]
#' 
#' CIGVARSelecte[which.min(CIGVARSelecte[,19]),]
#' @export
CIGVAR_Select = function(res,L_V=L_v)  {
Criteria = matrix(0,n*(L_V[1]-1)*(L_V[2]-1),4+n*3)
idx = 0
res_dd   = res
idx = 0
n   = res$n
for (I in 1:n )           {
for (l_d in 2: L_V[1] )   {
for (l_f in 2: L_V[2] )   {

      idx = idx + 1
      res_dd$p[I,1] = l_d 
      res_dd$p[I,2] = l_f       
	

    	res_s = CIGVARest(res_dd)
      Criteria[idx,] = c(as.vector(res_dd$p),res_s$AIC_g,res_s$BIC_g,sum(res_s$AIC),sum(res_s$BIC))
      #colnames(Criteria) = c("l_d","l_f","AIC_g","BIC_g","AIC","BIC")
    }
 }
}
return(Criteria)
}


#' Impulse Response Functions of GVAR 
#'
#' This function generates impulse response functions of an estimated GVAR with confidence bands 
#'
#' @param  res  : a list of the output of GVARest
#' @param  nstep: length of the impulse response functions
#' @param  comb : an mn vector specifying combined impulse such as global shocks, reginal shocks, or concerted actions.
#' @param  type : a list of the output of GVARest
#' @param  irf  : types of the impulse response irf=c("gen","chol","chol1","gen1","comb1"), gen for GIRF, gen1 for GIRF with unit impulse, chol for Cholezky decomposition, chol1 for Cholezky decomposition with unit impulse, comb1 for concerted action with unit impulse.
#' @param  runs : number of bootstrap runs to generate the confidence bands   
#' @param  conf :   
#' @return a matrix of (mn,mn,nstep,3) as the IRF colunms respesenting the impulse rows the responses. 
#' @examples 
#'
#'
#' n = 5
#' p = (1:15)*0; dim(p) = c(5,3)
#' p[,1] = 3; p[,2]=3; 
#' r_npo = (1:(2 * 3 * 5))*0; dim(r_npo) = c(2,3,5)
#' r_npo[,,1] = matrix(c(1,2,3,2,3,3),2,3)
#' r_npo[,,2] = matrix(c(1,2,3,2,3,3),2,3)
#' r_npo[,,3] = matrix(c(1,2,3,2,3,3),2,3)
#' r_npo[,,4] = matrix(c(1,2,3,2,3,3),2,3)
#' r_npo[,,5] = matrix(c(1,2,3,2,3,3),2,3)
#' 
#' res_d = CIGVARData(m=2,n=5,p=p,T=500,r_npo=r_npo,type="const",d=0)
#' max(res_d$Y)
#' plot(ts(res_d$Y))
#' res_e = CIGVARest(res_d)
#'
#' ### Impulse response function
#' 
#' IRF_CB = irf_CIGVAR_CB(res=res_e,nstep=20,comb=NA,irf="gen1",runs=200,conf=c(0.05,0.95))
#' 
#' par(mfrow=c(3,3))
#' plott(IRF_CB,1,1) 
#' plott(IRF_CB,1,2) 
#' plott(IRF_CB,1,3) 
#' plott(IRF_CB,1,4) 
#' plott(IRF_CB,1,5) 
#' plott(IRF_CB,1,6) 
#' plott(IRF_CB,1,7) 
#' plott(IRF_CB,1,8) 
#' plott(IRF_CB,1,9) 
#'
#'
#' @export
irf_CIGVAR_CB = function (res, nstep, comb, irf = c("gen", "chol", "chol1", "gen1", "comb1"), runs = 200, conf = c(0.05, 0.95)) 
{
    m = res$m
    n = res$n
    p = res$p
    T = dim(res$Y)[1]
    W = res$W
    Ao = res$Ao
    Bo = res$Bo
    Go = res$G
    Co = res$Co
    type  = res$type
    X     = res$X
    mu    = res$mu
    B     = res$G
    neq   = dim(B)[1]
    nvar  = dim(B)[2]
    sigma = res$Sigmao
    response <- array(0, dim = c(neq, nvar, nstep, length(conf) + 1))
    response[, , , 1] <- irf_GVAR(res, nstep, comb, irf)
    responseR <- array(0, dim = c(neq, nvar, nstep, runs))
    for (i in 1:runs) {
        Uo_run    = rnormSIGMA(T, sigma)
        res_run   = CIGVARData(m, n, p, T, W, r_npo = NA, Ao, Bo, Co, Uo = Uo_run, Sigmao = NA, type, X, mu)
        res_e     = CIGVARest(res_run)
        B_run     = res_e$G
        sigma_run = res_e$Sigmao
        responseR[, , , i] <- irf_GVAR(res_e, nstep, comb, irf)
    }
    for (tt in 1:(nstep)) {
        for (i in 1:neq) {
            for (j in 1:nvar) {
                Ao
                response[i, j, tt, -1] = quantile(responseR[i, 
                  j, tt, ], conf)
            }
        }
    }
    return(response)
}





