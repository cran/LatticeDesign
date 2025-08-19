#### Interleaved lattice-based maximin distance designs with continous, ordinal, and binary variables
#### Author: Blinded

######################## Algorithm 1 ########################################

Alg1_DeterFkSk <- function(L_ordinal) {
  
  # Algorithm A : for each ordinal dimension k, for each s_k = 2,...,length(ordinal),
  #               find all possible f, i.e. find s_k points from k-d ordinal list
  # input:
  # L_ordinal : list of ordinal variable;
  # p_ordinal : the number of ordinal dimensions
  
  e <- 10^-14
  p_ordinal <- dim(L_ordinal)[2]
  # min distance for discrete dimensions
  Di_min <- list()
  Di_2min <- list()
  De_good <- list()
  S_ordinal <- list()
  
  for(j in 1:p_ordinal){
    q <- L_ordinal[,j]
    len_q <- length(na.exclude(q))
    q <- L_ordinal[1:len_q,j]
    Di_min[[j]] <- list()
    Di_2min[[j]] <- list()
    De_good[[j]] <- list()
    s_select <- seq(2,len_q,1)
    
    for(i in s_select){
      d_min_0 <- vector(mode="numeric",length=0)  # the minimun of distance of any two adjacent points in one dimension
      d_2min_0 <- vector(mode="numeric",length=0)  # the minimun of distance of any two points separated by a point in one dimension
      De_good[[j]][[i]] <- list()  # store the f of j-th ordinal dimension, and s=s_select[i],
      # De_good's s_vector is correct, but Di_min and Di_2min's s_vector is not corresponce, need minus 1.
      if(i==2){
        d_min_0 <- q[len_q]
        d_2min_0  <- Inf
        De_good[[j]][[i]] <- c(De_good[[j]][[i]],list(c(0,q[len_q])))
      }
      if(i==3){
        dd <- abs(q-0.5)
        ddd <- q[which(dd==min(dd))][1]
        d_min_0 <- min(ddd,q[len_q]-ddd)
        d_2min_0 <- q[len_q]
        De_good[[j]][[i]] <-  c(De_good[[j]][[i]],list(c(0,ddd,q[len_q])))
      }
      
      if(i > 3){
        tobreak2 <- 0
        d_plus_1  <- 0
        while(tobreak2 == 0){
          d_star <- Alg2_dstar(q,s=i,d_plus_1)[[1]]
          d_plus_result <- Alg3_dplus(q,s=i,d_star)
          d_plus <- d_plus_result[[1]]
          if(d_plus < d_plus_1+e) {tobreak2 <- 1; break}
          d_plus_1 <- d_plus
          d_min_0 <- c(d_min_0,d_star)
          d_2min_0 <- c(d_2min_0,d_plus)
          De_good[[j]][[i]] <-  c(De_good[[j]][[i]],list(d_plus_result[[2]]))
        }
        
      }
      
      Di_min[[j]] <- c(Di_min[[j]],list(d_min_0))
      Di_2min[[j]] <- c(Di_2min[[j]],list(d_2min_0))
    } # end for S_select
    S_ordinal  <- c(S_ordinal ,list(s_select ))      # S_ordinal[[k]] stores all possible s in k-dimension
  } # end for ordinal column: jj
  
  return(list(De_good,Di_min,Di_2min,S_ordinal))
  
}

######################## Algorithm 2 ########################################

Alg2_dstar <- function(q,s,d_plus){
  # find the minimum distance btw 2 layers give the
  # q: the ordinal sequence
  # s: the number of unique values
  # d_plus: the minimum dist btw 3 layers
  
  q <- na.exclude(q)
  e<- 10^-14
  r<- rep(0,s)
  de <-q
  
  if(s==2) {d_star = q[length(q)]; r <- c(0,q[length(q)])}
  if(s==3) {
    dd <- abs(q-0.5)
    r[2] <- q[which(dd==min(dd))][1]
    r[3] <- q[length(q)]
    de <- r
    d_star <- min(diff(de))
  }
  if(s>3) {
    d <- min(diff(q))
    dobreak <- 0
    while(dobreak ==0){
      r[2] <- min(q[q>d])
      for(i in 3:s){
        q_down <- max(r[i-1]+d,r[i-2]+d_plus)
        r[i] <- min(q[q>q_down+e],2) # when q[q>q_down+e] is empty, it will be numeric(0), there will only left 2, then 2 will be the min, any other number bigger than 1 will be ok
        if(r[i]==2) break
      }
      if(r[i]==2){ d_star <- d;  dobreak<-1; break; }else{
        r[s] <- q[length(q)]
        de <- r
        d <- min(diff(de))
      }
    }
  }
  return(list(d_star,de))
  
}

######################## Algorithm 3 ########################################

Alg3_dplus <- function(q,s,d_star){
  # find the maximum distance of 3 layers in a sequcence r with s points
  # q: all possible value of the ordinal variable
  # s: the number of different values in the sequence
  # d_star: the maximum dist of 2 layers of r
  
  q <- na.exclude(q)
  e <- 10^-13
  e1 <- 10^-14
  r <- rep(0,s)
  de <-q
  
  if(s==2) {d_plus = Inf; de <- c(0,q[length(q)])}
  if(s==3) {
    dd <- abs(q-0.5)
    r[2] <- q[which(dd==min(dd))][1]
    r[3] <- q[length(q)]
    de <- r
    d_plus <- q[length(q)]
  }
  if(s>3) {
    d <- Inf
    for(i in 3:s) d <- min(d, q[i]-q[i-2])
    r[2] <- min(q[q>d_star-e])
    
    dobreak <- 0
    while(dobreak ==0){
      for(i in 3:s){
        q_down <- max(r[i-1] + d_star-e, r[i-2]+d)
        r[i] <- min(q[q>q_down+e1],2)
        if(r[i]==2) break
      }
      if(r[i]==2){ d_plus <- d;  dobreak<-1; break; }else{
        r[s] <- q[length(q)]
        de <- r
        d <-Inf
        for(i in 3:s) d <- min(d,de[i]-de[i-2])
      }
    }
  }
  return(list(d_plus,de))
}



######################## Algorithm 4 ########################################

Alg4_FindY <- function(discrete_dims,L01alt,is.rep,s_vector,weight,Di_min,Di_2min){
  #### given L&s, find f best
  #### discrete_dims: the index of discret dimensions
  
  p <- length(s_vector)
  
  len_odn <- length(discrete_dims)
  if (len_odn>0){
    con_dims <-setdiff(c(1:p), discrete_dims)
  }else{
    con_dims <-c(1:p)
  }
  
  # initialize v (consider weight)
  v <- rep(1,p)   # v is the vector of index for choosing which f in the f_1,...f_n.
  v[con_dims] <- weight[con_dims]/(s_vector[con_dims]-1)
  
  
  # initialize d_star & d_plus (don't consider weight)
  D_minn <- rep(0,p)
  D_2minn <- rep(0,p)
  
  if (len_odn>0){
    for(i in 1:len_odn){
      k <- discrete_dims[i]
      D_minn[k] <- Di_min[[i]][[s_vector[k]-1]][v[k]]
      D_2minn[k] <- Di_2min[[i]][[s_vector[k]-1]][v[k]]
    }
  }
  
  len_con <- p-len_odn
  if (len_con>0){
    for (i in 1:len_con){
      k <- con_dims[i]
      D_minn[k]<- 1/(s_vector[k]-1)
      if (s_vector[k]>2){
        D_2minn[k]<- 2/(s_vector[k]-1)
      }else{
        D_2minn[k]<- Inf
      }
    }
  }
  
  tobreak1 <- 0
  Sep_star <- 0
  #### find f best   ## Alg_B
  while(tobreak1==0){
    
    Sep_all <- Sep_Known_al1_ordinal(L01alt,is.rep,s_vector,D_minn,D_2minn,weight)
    Sep <- Sep_all[[1]]
    d1_min <- Sep_all[[2]]  #the first term in RHS of (4)
    d2_min <- Sep_all[[3]]   #the 2nd term in RHS of (4)
    
    if(Sep > Sep_star){
      Sep_star <- Sep
      v_star <- v  ## equivalent to f_T
    }else{
      tobreak1 <- 1; break
    }
    
    
    ### determin A
    B <- rep(0,p)
    for (j in 1:p) {
      ## find the dimensions which equals to the d_+_min
      if(Sep_star == D_2minn[j]*weight[j]) {B[j] <- 1}
    }
    A <- which(B==1)
    
    Ac <- intersect(A,con_dims)
    Ao <- setdiff(A,Ac)
    
    #### check odd s_k
    check1 <-0
    if (length(Ac)){ for(j in Ac){if (s_vector[j]%%2 !=0){check1 <- 1;break} } }
    if (check1) {tobreak1 <- 1; break}
    
    ### check avialibel f_k in A_o
    check3 <-0
    if (length(Ao)){
      for (k in 1:length(Ao)){
        jj <- Ao[k]
        kk <- which(discrete_dims==jj)
        if(is.na(Di_min[[kk]][[s_vector[jj]-1]][v[jj]+1])==0){
          v[jj] <- v[jj] +1
        }else{ check3 <- 1; break}
      }
    }
    if (check3) {tobreak1 <- 1; break}
    
    ### update dmin in K_o
    if (length(Ao)>0){
      for(i in 1:len_odn){
        k <- discrete_dims[i]
        D_minn[k] <- Di_min[[i]][[s_vector[k]-1]][v[k]]
        D_2minn[k] <- Di_2min[[i]][[s_vector[k]-1]][v[k]]
      }
    }
    
    
    #### if there are continue dimensions that are critical dimensions, do following,
    #### if not, skip the following.
    if (length(Ac)){
      
      ## compute new separation distance, rho_T
      Sep_all <- Sep_Known_al1_ordinal(L01alt,is.rep,s_vector,D_minn,D_2minn,weight)
      rho_t <- Sep_all[[1]]
      
      if(rho_t < Sep_star){
        tobreak1 <- 1; break
      }
      
      ### find all x in {0,1}^p
      # all_X = L0 + e_k \in L
      if (is.null(L01alt)){
        all_X <- diag(p)
      }else{
        all_X <- matrix(NA, dim(L01alt)[1], p)
        all_X[,is.rep != 1] = L01alt
        all_X[,is.rep == 1] = 0
        tmpid <- which(is.rep == 1)
        for(i in tmpid){
          tmpx <- rep(0,p)
          tmpx[i] <- 1
          all_X <- rbind(all_X, tmpx)
        }
      }
      
      #### remove the 0_p
      postion_nz <- rowSums(all_X)
      all_X_nz <- all_X[postion_nz!=0, drop= F]
      all_X_nz <- matrix(all_X_nz, ncol = p)
      
      ### determin X1 & X0
      rho_0 <- Inf
      X1 <- NULL
      X0 <- NULL
      for (j in 1:dim(all_X_nz)[1]){
        x <- all_X_nz[j,]
        check4 <- intersect(Ac,which(x==1))
        
        if (length(check4)){
          X1 <- rbind(X1,x)
        }else{
          X0 <- rbind(X0,x)
          rho_0 <- min(rho_0,sum((x*D_minn*weight)^2))
        }
      }
      
      X1 <- matrix(X1,ncol = p)
      
      ### compute rho_0
      Ac_bar <- setdiff(c(1:p),Ac)
      for(j in Ac_bar) if(is.rep[j]==0) if(s_vector[j]>2) rho_0 <- min(rho_0,(weight[j]*D_2minn[j])^2)
      rho_0 <- sqrt(rho_0)
      
      ### There are no values in the feasible domain
      up_alpha <- min(rho_0,min(2*weight[Ac]/(s_vector[Ac]-2)))
      lb_alpha <- rho_t
      if (up_alpha<lb_alpha){tobreak1 <- 1; break }
      
      #### solve the optimial problem
      lb <- c(lb_alpha, lb_alpha)
      ub <- c(up_alpha, up_alpha)
      
      # set the objective function & constraints & their grad
      fn <- function(x) return(-x[2])
      
      fn_grad <- function(x) return(c(0,-1))
      
      hin1 <- function(x) {
        
        consts <- rep(0,dim(X1)[1]+1)
        consts[1] <- x[2]-x[1]
        for(k in 2:(dim(X1)[1]+1))  {
          y <- X1[k-1,]
          consts[k] = x[2]^2 - sum(y[Ac]*(weight[Ac]-(s_vector[Ac]-2)/2*x[1])^2) - sum((y[Ac_bar]*D_minn[Ac_bar]*weight[Ac_bar])^2)
        }
        
        return(consts)
      }
      
      hin1_grad <- function(x) {
        consts_grad <- matrix(1,(dim(X1)[1]),2)
        for(k in 1:(dim(X1)[1])){
          y <- X1[k,]
          consts_grad[k,1] = sum( y[Ac]*(s_vector[Ac]-2)*(weight[Ac]-(s_vector[Ac]-2)/2*x[1]) )
        }
        consts_grad[,2] = 2*x[2]
        consts_grad <- rbind(c(-1,1),consts_grad)
        return(consts_grad)
      }
      
      x0 <-c(Sep_star,Sep_star)
      rests = nloptr::nloptr( x0,
                              eval_f=fn,
                              eval_grad_f=fn_grad,
                              lb = lb,
                              ub = ub,
                              eval_g_ineq = hin1,
                              eval_jac_g_ineq = hin1_grad,
                              opts = list("algorithm" = "NLOPT_LD_SLSQP",
                                          "xtol_rel"=1.0e-8) )#,
      
      opt_solution = rests$solution
      
      
      for (i in 1:length(Ac)){
        k <- Ac[i]
        v[k] <- weight[k]-(s_vector[k]-2)/2*opt_solution[1]
        D_minn[k]<- v[k]/weight[k]   #  1-(s_vector[k]-2)/2*opt_solution[1]
        D_2minn[k]<- opt_solution[1]/weight[k]
      }
    }
    
    
  } # end while
  
  return(list(Sep_star,v_star))
}


Sep_Known_al1_ordinal<- function(L01alt,is.rep,s_vector,D_minn,D_2minn,weight) {
  
  # Compute separation distance based on known G (and thus known L01alt).
  # v is the vector of the index for choosing which f in the f_1,...f_n.
  
  p <- length(s_vector)
  SepSquare <- sum(weight^2)
  
  d2_min <-min(weight*D_2minn)^2
  d1_min <- min((weight*D_minn))^2
  
  if(sum(is.rep)>0) SepSquare <- min(weight[is.rep==1]*D_minn[is.rep==1])^2
  if(sum(is.rep)<p) {
    for(j in 1:dim(L01alt)[2])  L01alt[,j] <- L01alt[,j] * weight[is.rep==0][j]
    for(i in 2:dim(L01alt)[1]) SepSquare <- min(SepSquare,sum((L01alt[i,]*D_minn[is.rep==0])^2))
    d1_min <- SepSquare
    d2_min <- Inf
    for(j in 1:p) if(is.rep[j]==0) if(s_vector[j]>2) d2_min <- min(d2_min,(weight[j]*D_2minn[j])^2)
    SepSquare <- min(d1_min,d2_min)
  }
  
  return(list(sqrt(SepSquare),sqrt(d1_min),sqrt(d2_min)))
}


######################## Algorithm 5 ########################################

Alg5_closetoN <- function(save_results,n){
  
  ### consider cosets, to find the design with the closest size to n
  ## compare all saved results,
  m_saved <- rep(0,length(save_results))
  is_n <- rep(0,length(save_results))
  for (i in 1:length(save_results)){
    m_saved[i]<-save_results[[i]][["m"]]
    is_n[i] = m_saved[i]==n
  }
  
  p <- length(save_results[[1]][["s_vector"]])
  
  flag <- sum(is_n)
  if (flag>0) {
    t_best <- which(is_n==1)[1]
    id_cosetBest <- 1
    n_best <- n
    L01altBest <- save_results[[t_best]][["L01alt"]]
    palt <- ncol(L01altBest)
    u_best <- rep(0,p)
    if(palt>0) u_best <- rep(0,palt)
    
  }else{
    
    dm0 <- Inf
    
    for (t in 1:length(save_results)){
      s_vectortmp <- save_results[[t]][["s_vector"]]
      p1tmp <-save_results[[t]][["p1"]]
      if(p1tmp==p){
        id_cosettmp <-  1
        mtmp <- save_results[[t]][["m"]]
        L01alttmp <- save_results[[t]][["L01alt"]]
        u_tmp <- rep(0,p)
      }else{
        L01alttmp <- save_results[[t]][["L01alt"]]
        is.reptmp <-save_results[[t]][["is.rep"]]
        tmps <- my_consider_coset(n,s_vectortmp,p-p1tmp,L01alttmp,is.reptmp)
        id_cosettmp <-  tmps[["id_cosetBest"]]
        mtmp <- tmps[["n_best"]]
        L01alttmp <- tmps[["L01altBest"]]
        u_tmp <- tmps[["u_best"]]
      }
      
      dmtmp <- mtmp-n
      if (dmtmp<0)  dmtmp <- Inf
      if (dmtmp < dm0){
        dm0 <- dmtmp
        t_best <- t
        L01altBest <- L01alttmp
        id_cosetBest <- id_cosettmp
        n_best <- mtmp
        u_best <- u_tmp
        
      }
    }
  }
  
  return(list(t_best=t_best,L01altBest=L01altBest,id_cosetBest=id_cosetBest,n_best=n_best,u_best=u_best))
  
}

my_consider_coset <- function(n,s_vector,pp,L01alttmp,is.reptmp){
  
  # index each point
  L01alttmp_idx <- cbind(L01alttmp,rep(0,dim(L01alttmp)[1]))
  for(i in 1:dim(L01alttmp)[1]){
    tmp <- 0
    for(j in pp:1 ){
      tmp <- tmp + L01alttmp[i,j]*2^{pp-j}
    }
    L01alttmp_idx[i,pp+1] <- tmp
  }
  
  # find all points in {0,1}^p
  is = rep(0,pp)
  allpoints <- matrix(0,2^pp,pp+1)
  for(i in 1:(2^pp)) {
    allpoints[i,] = c(is,i-1)
    is[pp]=is[pp]+1
    if(pp>1) for(j in pp:2) if(is[j]==2) { is[j]=0; is[j-1]=is[j-1]+1; }
  }
  
  difidx <- setdiff(allpoints[,pp+1], L01alttmp_idx[,pp+1]) # find the points not in design
  ini_idx <- L01alttmp_idx[,pp+1] # find the points in design as initialization
  
  # construct COEFD
  COEFD <- cbind(L01alttmp,rep(0,dim(L01alttmp)[1]))
  t <- 1
  while (length(difidx)>0) {
    tmp <- allpoints[difidx[1]+1,1:pp]
    L01alttmp_coset <- L01alttmp + matrix(rep(tmp,dim(L01alttmp)[1]),dim(L01alttmp)[1],pp,byrow = T)
    L01alttmp_coset[L01alttmp_coset==2]=0
    L01alttmp_coset <- unique(L01alttmp_coset)       ## remove repeated
    tmp2 <- cbind(L01alttmp_coset,rep(t,dim(L01alttmp_coset)[1]))
    COEFD <- rbind(COEFD,tmp2)
    
    L01alttmp_coset_idx <- rep(0,dim(L01alttmp_coset)[1])
    for(i in 1:dim(L01alttmp_coset)[1]){
      tmp <- 0
      for(j in pp:1 ){
        tmp <- tmp + L01alttmp_coset[i,j]*2^{pp-j}
      }
      L01alttmp_coset_idx[i] <- tmp
    }
    
    new_idx <- union(L01alttmp_coset_idx,ini_idx)
    difidx <- setdiff(allpoints[,pp+1], new_idx)
    ini_idx <- new_idx
    
    t <- t+1
  }
  
  # obtain COEFC
  p3 <- log2(t)
  p <- length(s_vector)
  p1 <- p - pp
  ps <- rep(0,3)
  ps[1] <- p1
  ps[2] <- p - p1 - p3
  ps[3] <- p3
  COEFC = supp2(COEFD,ps)
  
  # compute size
  lsss <- s_vector[is.reptmp==0]
  counts <- tnp1(lsss,COEFC)
  counts2 <- counts*prod(s_vector[is.reptmp==1]) #size
  near_n <- counts2-n
  near_n[which(near_n<0)] <- Inf
  id_coset<- which(near_n==min(near_n))
  id_cosetBest <- id_coset[1]     # the closest
  n_best <- counts2[id_cosetBest]
  u_best <- COEFD[which(COEFD[,pp+1]==id_cosetBest-1)[1],1:pp]
  
  L01alttmp <- COEFD[which(COEFD[,pp+1]==id_cosetBest-1),1:pp]
  
  return(list(id_cosetBest=id_cosetBest,u_best=u_best,n_best=n_best,L01altBest=L01alttmp))
  
  
}

supp2 = function(coefD,ps) {  # Supplement generating information.
  coefC=matrix(0,2^(ps[2]+ps[3]),2^ps[3])
  for(i in 1:(2^(ps[2]+ps[3]))) {
    corres=0
    for(j in (ps[2]+ps[3]):1) { if(coefD[i,j]==0) corres=c(corres,corres+2^(ps[2]+ps[3]-j)); }
    coefC[ corres+1, coefD[i,ps[2]+ps[3]+1]+1 ] = coefC[ corres+1, coefD[i,ps[2]+ps[3]+1]+1]+1
  }
  
  coefC
  
}
######################## Algorithm 6 ########################################
#### To generate the design in p dimensions and at least n points using the weight vector w and Algorithm 6 (p=2,3,4,5 only),
# run: ILMmDMixVarsAlg6(p,n,weight,discrete_dims,ordinal_levels) .


ILMmDMixVarsAlg6 <- function(p,n,discrete_dims,ordinal_levels,weight=rep(1,p)) {
  # THIS algorithm is to obtain the design with the highest separation distance when p=2,3,4,5.
  # p: the number of dimensions
  # n: the number of runs/points
  # weight: prior weight of variables
  # discrete_dims: the index of discrete variables in 1:p
  # ordinal_levels: a list conting the allowable levels of ordinal variables
  # NOTE: for continuous \gamma_k =[0,1], for ordinal \gamma_k = [l_{k,1}=0,...,l_{k,j_k}=1]
  
  ### load generator matrix
	if(p==2) GsName <- data(GeneratorMatrices2, envir=environment())
	if(p==3) GsName <- data(GeneratorMatrices3, envir=environment())
	if(p==4) GsName <- data(GeneratorMatrices4, envir=environment())
	if(p==5) GsName <- data(GeneratorMatrices5, envir=environment())
	Gs <- get(GsName)

  if(!p>=2|!p<=5|!p==round(p)) stop("p must be an integer greater than one and no greater than five.") 
  if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
  if(!length(weight)==p) stop("weight must be a vector with length p.") 
  
  
  ### load ordinal levels
  len_odn <- length(discrete_dims)
  if (len_odn){
    #load the ordinal list as a matrix
    if(is.character(ordinal_levels)) L_ordinal <- read.csv(file=paste(ordinal_levels,".csv",sep=''),header=T)
    if(is.list(ordinal_levels)){
      my_list = ordinal_levels
      length_diff <- max(sapply(my_list, length))
      my_list_diff_fixed <- lapply(my_list, function(x) {
        c(x, rep(NA, length_diff - length(x)))
      })
      L_ordinal <- data.frame(my_list_diff_fixed)
    }
    if(is.data.frame(ordinal_levels)) L_ordinal <- ordinal_levels
  }else{
    L_ordinal <- matrix(0,0,0)
  }
  
  
  Gsn <- dim(Gs)[1]/p
  SepBest <- 0
  GindexBest <- 0
  s_vectorBest <- rep(0,p)
  mBest <- 0
  e <- 10^-6
  MAX_row <- dim(L_ordinal)[1]
  len_odn <- length(discrete_dims)
  con_dims <- setdiff(c(1:p),discrete_dims)
  d1min <- 1e-6
  
  # the maximum of each s_k
  # if no ordinal variable
  s_lim <- rep(Inf,p)
  if (len_odn>0){
    for(i in 1:len_odn) {
      k <- discrete_dims[i]
      s_lim[k] <- MAX_row-sum(is.na(L_ordinal[,i]))
    }
  }
  
  # all are ordinal variable
  if(len_odn==p){
    m1 <- 1
    for(i in 1:p) m1 <- m1*s_lim[i]    # the max experiment runs
    if(m1 < n) { return(cat("The trials number beyonds the maximum design size!!!"))}
  }
  
  # run the algorithm A, find all possible combinations of ordinal values
  if (len_odn>0){
    ordinal_reultus <- Alg1_DeterFkSk(L_ordinal)
    De_good <- ordinal_reultus[[1]]  #layers
    Di_min <- ordinal_reultus[[2]]   #dstar
    Di_2min <- ordinal_reultus[[3]]  #dplus
    S_ordinal  <- ordinal_reultus[[4]] #S
    
    ## the max dplus of each dimension
    dplus_max_each <- list()
    for(i in 1:len_odn) {
      len_dplus <- length(Di_2min[[i]])
      dplus_max_each[[i]] <- matrix(0,len_dplus,1)
      for (j in 1:len_dplus){
        last_idx <- length(Di_2min[[i]][[j]])
        dplus_max_each[[i]][j,1]<- Di_2min[[i]][[j]][last_idx]
      }
    }
  }else{
    De_good <- list()
    Di_min <- list()
    Di_2min <- list()
    S_ordinal  <- list()
    dplus_max_each <- list()
  }
  
  
  # ergodic all generater matrixes
  show_G_sep_s_v <- rep(0,1,2+2*p)
  for(Gindex in 1:Gsn) {
    G <- Gs[(Gindex*p-p+1):(Gindex*p),]
    temp <- r_v(G)
    r_vector <- temp[[1]]
    #is.rep <- temp[[2]]
    is.rep <- unlist(temp[[2]], use.names = F)
    L01alt <- temp[[3]]
    p3     <- sum(G==2)
    s_vector <- rep(2,p)
    
    
    # the last dimension
    tmp <- max(c(ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, 2 ))
    if (tmp <= s_lim[p] ){
      s_vector[p] = tmp
    }else{
      s_vector[p] = 2
      for (k in (p-1):1) {if(s_vector[k] +1 <= s_lim[k]){break} }
      s_vector[k] = s_vector[k] +1
    }
    
    
    tobreak <- 0
    # ergodic all possible s
    
    while(tobreak==0) {
      
      m <- prod(ceiling(s_vector[is.rep==0]/2)*2)*prod(s_vector[is.rep==1])/2^p3  # A quick upper bound of m
      if(m>=n) m <- m_exact_known(G,r_vector,s_vector)
      if(m<n) {
        ##  find s which min(m(s)>n)
        ## check
        s_available <- rep(1,p)
        s_available[which(s_vector>=s_lim)] <- 0
        # if there are no s available, break this G
        if (sum(s_available)==0){ tobreak<-1; break;}
        s_avlb_idx <- rev(which(s_available==1))
        
        if( s_avlb_idx[1] ==p ) {
          s_vector[p] <- s_vector[p]+1
        }else{
          inc_idx = s_avlb_idx[1]
          s_vector[inc_idx] <- s_vector[inc_idx] +1
          s_vector[(inc_idx+1):p] = 2
        }
      }
      
      if(m>=n) {
        
        # Algrithm B: find the best separation distance & f given L & s
        tmp_rsults <- Alg4_FindY(discrete_dims,L01alt,is.rep,s_vector,weight,Di_min,Di_2min)
        Sep_T <- tmp_rsults[[1]]
        v_T <- tmp_rsults[[2]]
        show_G_sep_s_v <- rbind(show_G_sep_s_v,c(Gindex,Sep_T,s_vector,v_T))
        
        if(abs(Sep_T-SepBest)<1e-6) {
          SepBest <- Sep_T
          GindexBest <- Gindex
          s_vectorBest <- s_vector
          mBest <- m
          p3Best <- p3
          p1Best <- sum(is.rep)
          
          if(is.null(L01alt)) L01alt <- matrix(0,1,0)
          L01altBest <- L01alt
          is.repBest <- is.rep
          idbest <- v_T
          ### save this results
          t_save <- t_save + 1
          save_results[[t_save]] <- list(Sep=SepBest,L01alt=L01altBest,m=mBest,Gidx = GindexBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,id=idbest, is.rep= is.repBest)
          
        }
        
        ## update
        if(Sep_T>=SepBest+1e-6) {
          SepBest <- Sep_T
          GindexBest <- Gindex
          s_vectorBest <- s_vector
          mBest <- m
          p3Best <- p3
          p1Best <- sum(is.rep)
          if(is.null(L01alt)) L01alt <- matrix(0,1,0)
          L01altBest <- L01alt
          is.repBest <- is.rep
          idbest <- v_T
          ### save this results
          save_results <- list()
          t_save <- 1
          save_results[[t_save]] <- list(Sep=SepBest,L01alt=L01altBest,m=mBest,Gidx = GindexBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,id=idbest, is.rep= is.repBest)
          
          ## update s_lim
          if (len_odn>0){
            for (i in 1:len_odn){
              k <- discrete_dims[i]
              s_abdn <- which(dplus_max_each[[i]]<SepBest) #
              if (length(s_abdn) !=0){s_lim[k] <- min(s_abdn)}
            }
          }
          
          if (len_odn<p){
            
            for (k in con_dims) {
              tmp <- floor(2*(weight[k])/SepBest + 1)
              s_lim[k] <- tmp
              if (floor(tmp/2)*2!=tmp){
                tmp2 <- 2*(weight[k]-d1min)/(tmp+1-2)
                if(tmp2>=SepBest) {
                  s_lim[k] <- tmp+1
                }
              }
            }
            
          }
          
        }  # end if(Sep_T>SepBest)
        
        
        ## gradually increase s_vector
        s_vector[p] <- s_lim[p]
        
        s_available <- rep(1,p)
        s_available[which(s_vector>=s_lim)] <- 0
        
        # if there are no possible s, change to next generator
        if (sum(s_available)==0){tobreak<-1; break;} else{
          for(kk in p:1) if(s_vector[kk]<s_lim[kk])   break;
        }
        
        s_vector[kk] <- s_vector[kk]+1;
        s_vector[(kk+1):p] <- 2
        
      }  # end if m>n
      
    } # end while
  } # end all generator matrixes
  
  # write.xlsx(show_G_sep_s_v, file = "exp1.xlsx", row.names = FALSE, sheetName = "Sheet1")
  
  
  bestresults <- Alg5_closetoN(save_results,n)
  t_best <- bestresults[["t_best"]]
  id_cosetBest <- bestresults[["id_cosetBest"]]
  n_best <- bestresults[["n_best"]]
  L01altBest <- bestresults[["L01altBest"]]
  uuBest <- bestresults[["u_best"]]
  
  
  ## choose the best one
  SepBest <- save_results[[t_best]][["Sep"]]
  GindexBest <- save_results[[t_best]][["Gidx"]]
  s_vectorBest <- save_results[[t_best]][["s_vector"]]
  mBest <- n_best
  p3Best <- save_results[[t_best]][["p3"]]
  p1Best <- save_results[[t_best]][["p1"]]
  is.repBest <- save_results[[t_best]][["is.rep"]]
  idbest <- save_results[[t_best]][["id"]]
  
  uBest <- rep(0,p)
  uBest[is.repBest==0] <- uuBest
  
  #generate Design
  #use L01 to generate the design
  L01Best <- matrix(0,dim(L01altBest)[1],p)
  L01Best[,is.repBest==0] <- L01altBest
  Design <- L01Best
  L01BestS <- L01Best
  for(j in p:1) {
    k <- j
    if(is.repBest[k]==0) if(s_vectorBest[j]>3) for(i in 2:floor(s_vectorBest[j]/2)) {
      DesignAdd <- L01Best
      DesignAdd[,j] <- DesignAdd[,j]+(i-1)*2
      Design <- rbind(Design,DesignAdd)
    }
    if(is.repBest[k]==0) if(floor(s_vectorBest[j]/2)*2 != s_vectorBest[j]) {
      DesignAdd <- L01Best[L01Best[,j]==0,,drop=FALSE]
      DesignAdd[,j] <- DesignAdd[,j] + s_vectorBest[j]-1
      Design <- rbind(Design,DesignAdd)
    }
    if(is.repBest[k]==1) for(i in 2:s_vectorBest[j]) {
      DesignAdd <- L01Best
      DesignAdd[,j] <- DesignAdd[,j]+(i-1)
      Design <- rbind(Design,DesignAdd)
    }
    L01Best <- Design
  }
  
  # debug 2024-03-31:discrete_dims is 0-1 sequence, and now it is a index set
  h=1
  for(j in discrete_dims){
    D_ordinal <- c(De_good[[h]][[s_vectorBest[j]]][[idbest[j]]])
    for(ss in 0:(s_vectorBest[j]-1)) Design[which(Design[,j]==ss),j] <- D_ordinal[ss+1]
    h=h+1
  }
  for(j in con_dims) Design[,j] <- Design[,j]/(s_vectorBest[j]-1)
  DesignTransformed <- Design
  for(j in 1:p) DesignTransformed[,j] <- DesignTransformed[,j]*weight[j]
  
  #return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,Gidx = GindexBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,uBest = uBest,L01=L01BestS,id=idbest,cosetId =id_cosetBest, Nbest = n_best))
  return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,s_vector=s_vectorBest,L01=L01BestS,d_star=idbest,uBest = uBest,ordinal_levels = ordinal_levels))
  
}


######################## Algorithm 7 ########################################

Alg7_construcL <- function(p,p3,VectorInConsideration,VectorListCannot,wdsOrdered,is.rep,Sep){
  
  p1 <- sum(is.rep) 
  
  ## construct L0
  YesAlready <- 0
  VectorListNo <- VectorListCannot[,!is.rep,drop=FALSE] ## select the non rep column
  for(i in dim(VectorListNo)[1]:1) if(sum(VectorListNo[i,])==0) VectorListNo <- VectorListNo[-i,,drop=FALSE]##remove the row with all 0
  VectorListYes <- matrix(0,0,p-p1)  ##L0L
  VectorGeneratorAlt <- VectorListYes
  VectorListInfo <- VectorInConsideration
  for(j in p:1) if(is.rep[j]==1) VectorListInfo <- VectorListInfo[,-j,drop=FALSE]## remove the kth column if ek==1
  for(i in dim(VectorListInfo)[1]:1) if(sum(VectorListInfo[i,1:(p-p1)])==0) VectorListInfo <- VectorListInfo[-i,,drop=FALSE]#remove the row with all 0
  while(dim(VectorListYes)[1]<2^(p-p1-p3)-1) {
    VectorToNoIndex <- order(VectorListInfo[,p-p1+1])[1]
    SepN <- min(VectorListInfo[,p-p1+1])
    VectorToNo <- VectorListInfo[VectorToNoIndex,1:(p-p1)]
    for(kkk in 1:(p-p1)) if(VectorToNo[kkk]>0) break
    if(kkk-1>0) for(i in 1:(kkk-1)) {
      TempVector <- VectorToNo
      TempVector[i] <- TempVector[i]+1
      TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered[!is.rep])^2)))
      VectorListInfo <- rbind(VectorListInfo,TempVector)
    }
    VectorListInfo <- VectorListInfo[-VectorToNoIndex,,drop=FALSE]
    NeedTo <- 1
    if(dim(VectorListYes)[1]>0) for(i in 1:dim(VectorListYes)[1]) if(sum(VectorToNo==VectorListYes[i,])==p-p1) NeedTo <- 0
    if(dim(VectorListYes)[1]==0) for(i in 1:dim(VectorListNo)[1]) if(sum(VectorToNo==VectorListNo[i,])==p-p1) NeedTo <- 0
    if(dim(VectorListYes)[1]>0) for(i in 1:dim(VectorListNo)[1]) for(j in 1:dim(VectorListYes)[1]) if(sum(VectorToNo==BAdd(VectorListNo[i,],VectorListYes[j,]))==p-p1) NeedTo <- 0
    if(NeedTo==1) {
      CanNo <- ( CheckArrange_old(p-p1,p3,rbind(VectorListNo,VectorToNo),VectorListYes) >0 )
      if(CanNo) VectorListNo <- rbind(VectorListNo,VectorToNo)
      if(!CanNo) {
        VectorListYes <- AddToGroup_old(VectorListYes,VectorToNo)
        VectorGeneratorAlt <- rbind(VectorGeneratorAlt,VectorListYes)
        if(YesAlready==0) {
          Sep <- min(c( Sep, SepN ))  # the first vector added from Gamma to L0
          YesAlready <- 1
        }
      }
    }
  }
  
  L01alt <- rbind(rep(0,p-p1),VectorListYes)  ## VectorGenerator and L01 not updated.
  
  
  return(list(L01alt = L01alt, VectorListCannot=VectorListCannot,Sep=Sep, VectorGeneratorAlt=VectorGeneratorAlt))
}


CheckArrange_old <- function(p,p3,VectorListCannot,VectorListYes)  {
  # Return 0 if cannot arrange; return number of (nonzero) cosets if can arrange.
  # VectorListYes must have group structure; VectorListCannot needs not.
  # VectorListYes can have zero rows; VectorListCannot needs at least one row.
  StatusList <- rep(-1,2^p-1)
  if(dim(VectorListCannot)[1]>0) for(i in 1:dim(VectorListCannot)[1]) StatusList[SiteBinary(VectorListCannot[i,])] <- 2^p3
  VectorListArrange <- rep(list(matrix(0,0,p)),2^p3-1)
  if(dim(VectorListCannot)[1]>0) for(i in 1:dim(VectorListCannot)[1]) {
    VectorToArrange <- VectorListCannot[i,]
    if(StatusList[SiteBinary(VectorToArrange)]==0) return(0)
    if(StatusList[SiteBinary(VectorToArrange)]!=2^p3) next
    j=1; while(j<=2^p3-1) {
      AddHere <- 1
      if(dim(VectorListArrange[[j]])[1]==0) break
      VectorListToYes <- BAdds_old(VectorListArrange[[j]],VectorToArrange)
      StatusListNew <- StatusList
      VectorListNewYes <- VectorListYes
      for(ii in 1:dim(VectorListToYes)[1]) {
        temp <- AddToGroup_old(VectorListNewYes,VectorListToYes[ii,])
        if(dim(temp)[1]>dim(VectorListNewYes)[1]) {
          VectorListNewYes <- temp    # H0
        }
      }
      for(iii in 1:dim(VectorListNewYes)[1]) {
        if(StatusListNew[SiteBinary(VectorListNewYes[iii,])]>0) { AddHere <- 0; break; }
        StatusListNew[SiteBinary(VectorListNewYes[iii,])] = 0
      }
      for(jj in 1:(2^p3-1)) if(dim(VectorListArrange[[jj]])[1]>0) {
        if(AddHere==0) break
        VectorListToCoset <- BAddss_old(VectorListArrange[[jj]],VectorListNewYes)
        for(iii in 1:dim(VectorListToCoset)[1]) {
          OldStatus <- StatusListNew[SiteBinary(VectorListToCoset[iii,])]
          if(OldStatus!=jj & OldStatus!=2^p3 & OldStatus!=-1) { AddHere <- 0; break; }
          StatusListNew[SiteBinary(VectorListToCoset[iii,])] <- jj
        }
      }
      if(AddHere==1) break
      j=j+1
    }
    if(j==2^p3) return(0)
    VectorListArrange[[j]] <- rbind(VectorListArrange[[j]],VectorToArrange)
    StatusList[SiteBinary(VectorToArrange)] <- j
    if(dim(VectorListArrange[[j]])[1]>1) {
      VectorListYes <- VectorListNewYes
      StatusList <- StatusListNew
      for(jj in 1:(2^p3-1)) if(dim(VectorListArrange[[jj]])[1]>0) {
        VectorListToCoset <- BAddss_old(VectorListArrange[[jj]],VectorListYes)
        for(iii in 1:dim(VectorListToCoset)[1]) {
          if(StatusList[SiteBinary(VectorListToCoset[iii,])]!=jj) {
            VectorListArrange[[jj]] <- rbind(VectorListArrange[[jj]],VectorListToCoset[iii,])
            StatusList[SiteBinary(VectorListToCoset[iii,])] <- jj
          }
        }
      }
    }
  }
  for(j in 1:(2^p3-1)) if(dim(VectorListArrange[[j]])[1]==0) return(j-1)
  return(2^p3-1)
}

######################## Algorithm 8 ########################################


#### To generate the design in p dimensions and at least n points using the weight vector w and Algorithm 8 (p=6,7,8 only),
# run: ILMmDMixVarsAlg8(p,n,weight,discrete_dims,ordinal_levels,ReturnAll = 0) .

ILMmDMixVarsAlg8 <- function(p,n,discrete_dims,ordinal_levels,weight=rep(1,p),N = 1,ReturnAll = 0) {
  # THIS algorithm is to obtain the design with the highest separation distance when p=6,7,8.
  # p: the number of dimensions
  # n: the number of runs/points
  # weight: prior weight of variables
  # discrete_dims: the index of discrete variables in 1:p
  # ordinal_levels: a list conting the allowable levels of ordinal variables
  # N: the top N results will be saved
  # ReturnAll: if ReturnAll = True, it will output all designs with the highest separation distance, not considering n.
  # NOTE: for continuous \gamma_k =[0,1], for ordinal \gamma_k = [l_{k,1}=0,...,l_{k,j_k}=1]
  
  if(!p>=2|!p==round(p)) stop("p must be an integer greater than one.") 
  if(!p<=8) warning("The program can be extremely slow for p>8.")
  if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
  if(!length(weight)==p) stop("weight must be a vector with length p.") 
  
  
  
  len_odn <- length(discrete_dims)
  if (len_odn){
    #load the ordinal list as a matrix
    if(is.character(ordinal_levels)) L_ordinal <- read.csv(file=paste(ordinal_levels,".csv",sep=''),header=T)
    if(is.list(ordinal_levels)){
      my_list = ordinal_levels
      length_diff <- max(sapply(my_list, length))
      my_list_diff_fixed <- lapply(my_list, function(x) {
        c(x, rep(NA, length_diff - length(x)))
      })
      L_ordinal <- data.frame(my_list_diff_fixed)
    }
    if(is.data.frame(ordinal_levels)) L_ordinal <- ordinal_levels
  }else{
    L_ordinal <- matrix(0,0,0)
  }
  
  # initialize some variables
  d1min <- 1e-6
  MAX_row <- dim(L_ordinal)[1]
  slmt_tmp <- rep(MAX_row,len_odn) - colSums(is.na(L_ordinal))
  s_lim <- rep(Inf,p)
  s_lim[discrete_dims] <- slmt_tmp
  is_odn <- rep(0,p)
  is_odn[discrete_dims] <-1
  con_dims <- setdiff(c(1:p),discrete_dims)
  SepBest <- 0
  
  # run the algorithm A, find all possible combinations of ordinal values
  if (len_odn>0){
    ordinal_reultus <- Alg1_DeterFkSk(L_ordinal)
    De_good <- ordinal_reultus[[1]]
    Di_min <- ordinal_reultus[[2]]
    Di_2min <- ordinal_reultus[[3]]
    S_ordinal  <- ordinal_reultus[[4]]
    
    ## the max dplus of each dimension
    dplus_max_each <- list()
    for(i in 1:len_odn) {
      len_dplus <- length(Di_2min[[i]])
      dplus_max_each[[i]] <- matrix(0,len_dplus,1)
      for (j in 1:len_dplus){
        last_idx <- length(Di_2min[[i]][[j]])
        dplus_max_each[[i]][j,1]<- Di_2min[[i]][[j]][last_idx]
      }
    }
  }else{
    De_good <- list()
    Di_min <- list()
    Di_2min <- list()
    S_ordinal  <- list()
    dplus_max_each <- list()
  }
  
  
  ######### save the first 20th seperation distance results
  #N <- 20
  save_results_bestN <- list()
  bestNth <- 0
  sep_bestN <- rep(-Inf,N)
  sep_bestN_lowest <- 0
  lowest_indx <- NA
  #isfull_bestN <- 0
  
  ##### Initialize separation distance ##########
  ##### Fast solution for p3=1, all alt., not guaranteed best s_vector.
  ## all layers are alternative layers
  
  p3 <- 1   # q = p-p3
  
  s_vector=rep(2,p)
  is_lmt <- s_vector==s_lim
  m <- ceiling(prod(s_vector)/2)
  
  SepEachNext <- rep(0,p)
  while(m<n && sum(is_lmt)!=p) {
    ##  find s which min(m(s)>n)
    SepEachNext[con_dims] <- weight[con_dims]/(s_vector[con_dims])
    if (len_odn>0){
      for(i in 1:len_odn){
        odn_idx <- discrete_dims[i]
        if (is_lmt[odn_idx]==0) {
          SepEachNext[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]]][1]
        }
      }
    }
    
    SepEachNext[is_lmt] <- 0
    s_vector[order(-SepEachNext)[1]] <- s_vector[order(-SepEachNext)[1]]+1
    m <- ceiling(prod(s_vector)/2)
    is_lmt <- s_vector==s_lim
  }
  
  if(m>=n){
    
    # construct L01alt
    L01alt <- matrix(0,1,p)
    for(j in p:1) { ## construct p*p matrix with 0, 1
      L01altAdd <- L01alt; L01altAdd[,j]<-L01altAdd[,j]+1; L01alt <- rbind(L01alt,L01altAdd); }
    temp <- rep(0,2^p)
    for(j in 1:p) temp <- temp + L01alt[,j] #sum row
    L01alt <- L01alt[floor(temp/2)*2==temp,] ##if sum row = even, to be L0
    is.rep <- rep(0,p)
    
    # compute d^* and d^+
    v_T <- rep(1,p)
    SepEach <- rep(0,p)
    SepEachPlus <- rep(0,p)
    SepEach[con_dims] <- weight[con_dims]/(s_vector[con_dims]-1) # d^*
    SepEachPlus[con_dims] <- 2*weight[con_dims]/(s_vector[con_dims]-1) # d^+
    v_T[con_dims] <- SepEach[con_dims]
    
    if (len_odn>0){
      for(i in 1:len_odn){
        odn_idx <- discrete_dims[i]
        SepEach[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]-1]][1]
        SepEachPlus[odn_idx] <- weight[odn_idx]*Di_2min[[i]][[s_vector[odn_idx]-1]][1]
      }
    }
    
    # compute sep
    # Algrithm B: find the best separation distance & f given L & s
    tmp_rsults <- Alg4_FindY(discrete_dims,L01alt,is.rep,s_vector,weight,Di_min,Di_2min)
    SepUpperBound <- tmp_rsults[[1]]
    v_T <- tmp_rsults[[2]]
    
    if(SepUpperBound>SepBest) {
      SepBest <- SepUpperBound
      dmBest <- m-n
      mBest <- m
      p3Best <- p3
      p1Best <- 0
      s_vectorBest <- s_vector
      NCannotBest <- sum(SepEach<SepUpperBound)
      DimOrderBest <- 1:p
      is.repBest <- is.rep
      L01altBest <- L01alt
      VectorGeneratorAltBest <- cbind(rep(1,p-1),diag(p-1))
      idbest <- v_T
      ### save this results
      L01Best <- matrix(0,dim(L01altBest)[1],p)
      L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
      save_results <- list()
      t_save <- 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    }
    
    # save
    bestNth = bestNth + 1
    save_results_bestN[[bestNth]] <-list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    sep_bestN[bestNth] <- SepUpperBound
    sep_bestN_lowest <- min(sep_bestN[is.finite(sep_bestN)])
    lowest_indx <- which(sep_bestN_lowest==sep_bestN)[1]
    # isfull_bestN <- bestNth==N
    
    
  }
  
  
  ##### Fast solution for p3=p-1, all alt., not guaranteed best s_vector.
  p3 <- p-1   #q=1
  
  s_vector=rep(2,p)
  is_lmt <- s_vector==s_lim
  
  L01alt <- matrix(0,1,p)
  L01alt <- rbind(L01alt,rep(1,p))
  is.rep <- rep(0,p)
  m <- m_exact_fromL01(L01alt,is.rep,s_vector)
  
  
  SepEachNext <- rep(0,p)
  while(m<n && sum(is_lmt)!=p) {
    ##  find s which min(m(s)>n)
    SepEachNext[con_dims] <- weight[con_dims]/(s_vector[con_dims])
    if (len_odn>0){
      for(i in 1:len_odn){
        odn_idx <- discrete_dims[i]
        if (is_lmt[odn_idx]==0) SepEachNext[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]]][1]
      }
    }
    SepEachNext[is_lmt] <- 0
    s_vector[order(-SepEachNext)[1]] <- s_vector[order(-SepEachNext)[1]]+1
    m <- m_exact_fromL01(L01alt,is.rep,s_vector)
    is_lmt <- s_vector==s_lim
  }
  
  
  if(m>=n){
    v_T <- rep(1,p)
    SepEach <- rep(0,p)
    SepEachPlus <- rep(0,p)
    SepEach[con_dims] <- weight[con_dims]/(s_vector[con_dims]-1)
    SepEachPlus[con_dims] <- 2*weight[con_dims]/(s_vector[con_dims]-1)
    v_T[con_dims] <- SepEach[con_dims]
    
    if (len_odn>0){
      for(i in 1:len_odn){
        odn_idx <- discrete_dims[i]
        SepEach[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]-1]][1]
        SepEachPlus[odn_idx] <- weight[odn_idx]*Di_2min[[i]][[s_vector[odn_idx]-1]][1]
      }
    }
    
    # Algrithm B: find the best separation distance & f given L & s
    tmp_rsults <- Alg4_FindY(discrete_dims,L01alt,is.rep,s_vector,weight,Di_min,Di_2min)
    SepUpperBound <- tmp_rsults[[1]]
    v_T <- tmp_rsults[[2]]
    
    if(abs(SepUpperBound-SepBest)<1e-4) {
      SepBest <- SepUpperBound
      dmBest <- m-n
      mBest <- m
      p3Best <- p3
      p1Best <- 0
      s_vectorBest <- s_vector
      NCannotBest <- sum(SepEach<SepUpperBound)
      DimOrderBest <- 1:p
      is.repBest <- is.rep
      L01altBest <- L01alt
      VectorGeneratorAltBest <- rep(1,p)
      idbest <- v_T
      ### save this results
      L01Best <- matrix(0,dim(L01altBest)[1],p)
      L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
      t_save <- t_save + 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    }
    
    if(SepUpperBound>=SepBest+1e-6) {
      SepBest <- SepUpperBound
      dmBest <- m-n
      mBest <- m
      p3Best <- p3
      p1Best <- 0
      s_vectorBest <- s_vector
      NCannotBest <- sum(SepEach<SepUpperBound)
      DimOrderBest <- 1:p
      is.repBest <- is.rep
      L01altBest <- L01alt
      VectorGeneratorAltBest <- rep(1,p)
      idbest <- v_T
      L01Best <- matrix(0,dim(L01altBest)[1],p)
      L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
      ### save this results
      save_results <- list()
      t_save <- 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    }
    
    # save
    bestNth = bestNth + 1
    save_results_bestN[[bestNth]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    sep_bestN[bestNth] <- SepUpperBound
    sep_bestN_lowest <- min(sep_bestN[is.finite(sep_bestN)])
    lowest_indx <- which(sep_bestN_lowest==sep_bestN)[1]
    #isfull_bestN <- bestNth==N
    
  }
  
  ##### Exact solution for p3=0. all repeated layer
  p3 <- 0
  s_vector <- rep(2,p)
  is_lmt <- s_vector==s_lim
  m <- prod(s_vector)
  
  SepEachNext <- rep(0,p)
  while(m<n && sum(is_lmt)!=p) {
    SepEachNext[con_dims] <- weight[con_dims]/(s_vector[con_dims])
    if (len_odn>0){
      for(i in 1:len_odn){
        odn_idx <- discrete_dims[i]
        if (is_lmt[odn_idx]==0) SepEachNext[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]]][1]
      }
    }
    # find the k maximizes d*_k,sk+1 max, let s_k = s_k +1
    SepEachNext[is_lmt] <- 0
    s_vector[order(-SepEachNext)[1]] <- s_vector[order(-SepEachNext)[1]]+1
    m <- prod(s_vector)
    is_lmt <- s_vector==s_lim
  }
  
  if(m>=n){
    
    SepEach <- rep(0,p)
    SepEach[con_dims] <- weight[con_dims]/(s_vector[con_dims]-1)
    if (len_odn>0){
      for(i in 1:len_odn){
        odn_idx <- discrete_dims[i]
        SepEach[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]-1]][1]
      }
    }
    SepUpperBound <- min(SepEach)
    
    
    if(abs(SepUpperBound-SepBest)<1e-4) {
      SepBest <- SepUpperBound
      dmBest <- m-n
      mBest <- m
      p3Best <- p3
      p1Best <- p
      s_vectorBest <- s_vector
      NCannotBest <- 0
      DimOrderBest <- 1:p
      is.repBest <- rep(1,p)
      L01altBest <- matrix(0,1,0)
      VectorGeneratorAltBest <- matrix(0,0,0)
      idbest <- rep(1,p)
      L01Best <- matrix(0,dim(L01altBest)[1],p)
      L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
      ### save this results
      t_save <- t_save + 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
      
    }
    
    if(SepUpperBound>SepBest) {
      SepBest <- SepUpperBound
      dmBest <- m-n
      mBest <- m
      p3Best <- p3
      p1Best <- p
      s_vectorBest <- s_vector
      NCannotBest <- 0
      DimOrderBest <- 1:p
      is.repBest <- rep(1,p)
      L01altBest <- matrix(0,1,0)
      VectorGeneratorAltBest <- matrix(0,0,0)
      idbest <- rep(1,p)
      L01Best <- matrix(0,dim(L01altBest)[1],p)
      L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
      ### save this results
      save_results <- list()
      t_save <- 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    }
    
    # save
    bestNth = bestNth + 1
    save_results_bestN[[bestNth]] <-list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01alt=L01altBest,L01=L01Best,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
    sep_bestN[bestNth] <- SepUpperBound
    sep_bestN_lowest <- min(sep_bestN[is.finite(sep_bestN)])
    lowest_indx <- which(sep_bestN_lowest==sep_bestN)[1]
    #isfull_bestN <- bestNth==N
    
  }
  
  ##########   End initialize separation distance
  
  ## update s_limit by rho_B
  ### s_k such that d+(s_k) > rho_B
  v_lim <- rep(0,p)
  if (len_odn>0){
    for (i in 1:len_odn){
      k <- discrete_dims[i]
      s_abdn <- which(weight[k]*dplus_max_each[[i]]<SepBest)
      if (length(s_abdn) !=0){
        s_lim[k] <- min(s_abdn)
      }
      tmp <- which(unlist(weight[k]*Di_2min[[i]][[s_lim[k]-1]])>=SepBest)
      v_lim[k] <- tmp[1]
    }
  }
  
  if (len_odn<p){
    for (k in con_dims) {
      tmp <- floor(2*(weight[k])/SepBest + 1)
      s_lim[k] <- tmp
      v_lim[k] <- (weight[k])/(s_lim[k]-1)
      if (floor(tmp/2)*2!=tmp){
        tmp2 <- 2*(weight[k]-weight[k]*d1min)/(tmp+1-2)
        if(tmp2>=SepBest) {
          s_lim[k] <- tmp+1
          v_lim[k] <- max(c(weight[k]*d1min,weight[k] - (SepBest)*(s_lim[k]-2)/2))
        }
      }
    }
  }
  
  ###########################################################################################
  ## Search from p3=1. Only consider best Sep for given p3 and s_vector.                 ###
  ###########################################################################################
  changeL <- 0
  checkloop <- 1
  p3_svector <- matrix(0,0,p+1)
  
  for(p3 in 1:(p-1)){
    # choose s
    s_vector <- rep(2,p)
    temp =  min(c(ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, s_lim[p] )) # A quick upper bound of m: prod(ceiling(s_vector[-p]/2)*2)/2^p3
    s_vector[p] = max(c(temp, 2 ))
    is_lmt <- s_vector==s_lim
    
    tobreak <- 0
    while(tobreak==0) {  # m<n: increase last s; IsChance==0: increase other s and set last as 2.
      
      # determine y
      v_T <- rep(1,p)
      v_T[con_dims] <- weight[con_dims]/(s_vector[con_dims]-1)
      if (1){  # =0, only consider the layers with the max d*
        if (sum(is_lmt)){
          if(sum(is_lmt[discrete_dims])){
            tmpid <- discrete_dims[is_lmt[discrete_dims]==1]
            v_T[tmpid] <- v_lim[tmpid]
          }
          if(sum(is_lmt[con_dims])){
            tmpid <- con_dims[is_lmt[con_dims]==1]
            v_T[tmpid] <- v_lim[tmpid]
          }
        }
      }
      
      # determine d^*
      SepEach <- rep(0,p)
      SepEach[con_dims] <-v_T[con_dims]
      if (len_odn>0){
        for(i in 1:len_odn){
          odn_idx <- discrete_dims[i]
          SepEach[odn_idx] <- weight[odn_idx]*Di_min[[i]][[s_vector[odn_idx]-1]][v_T[odn_idx]]
        }
      }
      
      
      DimMustAlt <- as.vector(SepEach<SepBest) ## non-rep
      if (sum(DimMustAlt)<=p3) {
        # Force some dimensions to be Alt . no-rep dimensions should larger than p3
        # if s_k is odd and SepEach_k > SepBest, make it be alt layer
        DimSeq <- s_vector + (floor(s_vector/2)*2==s_vector)*max(s_vector) + (SepEach<SepBest)*max(s_vector)*2
        DimMustAlt[order(DimSeq)[1:(p3+1-sum(DimMustAlt))]] = TRUE
      }
      
      # find a upper bound on m, if it smaler than n, then increase s
      m <- ceiling(prod(s_vector[DimMustAlt])/2)
      if(p3>=2) m <- m - prod(floor(s_vector[DimMustAlt]/2)*2)/2 + prod(floor(s_vector[DimMustAlt]/2)*2)/2^p3
      if((!is.na(sum(DimMustAlt)<=p)) && (sum(DimMustAlt)<=p)) if(sum(DimMustAlt)>0) m <- m * prod(s_vector[DimMustAlt==0])  # Another upper bound of m, good for low p3. p3=1 exact; p3=2 very good.
      if(m<n) {
        ##  find s which min(m(s)>n)
        # if there are no s available, try next p3
        
        if (sum(is_lmt)==p){tobreak<-1; break;}
        s_avlb_idx <- rev(which(is_lmt==0))
        
        if( s_avlb_idx[1] ==p ) {
          s_vector[p] <- s_vector[p]+1
        }else{
          inc_idx = s_avlb_idx[1]
          s_vector[inc_idx] <- s_vector[inc_idx] +1
          s_vector[(inc_idx+1):p] = 2
          temp =  min(c(ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, s_lim[p] ))  # A quick upper bound of m: prod(ceiling(s_vector[-p]/2)*2)/2^p3
          s_vector[p] = max(c(temp, 2 ))
        }
        
        is_lmt <- s_vector==s_lim
        next;
      }
      
      
      
      IsChance <- 1  # Is the choice of p3 and s_vector has a chance to be SepBest? 0 stands for no chance.
      # give a sep upper bound, if upper bound< sep, no need to do
      SepUpperBound <- sqrt(sum(weight^2))
      SepEachPlusmax <- max(weight)
      #Find the largest possible d+and d*
      if(max(s_vector)>2) {
        # ordinal & s_k>2
        tmpidx <- which((s_vector>2)+is_odn==2)
        if(length(tmpidx)>0){
          for(i in 1:length(tmpidx)) {
            k <- which(discrete_dims==tmpidx[i])
            tmpid <- tmpidx[i]
            SepEachPlusmax <-min(c(SepEachPlusmax, weight[tmpid]*Di_2min[[k]][[s_vector[tmpid]-1]][v_T[tmpid]]))
          }
        }
        # continuous
        tmpidx <- which((s_vector>2)+(!is_odn)==2)
        # continuous & even
        tmpidx2 <- tmpidx[floor(s_vector[tmpidx]/2)*2==s_vector[tmpidx]]
        if(length(tmpidx2)) {
          SepEachPlusmax <- min(c(SepEachPlusmax,2*(weight[tmpidx2]-v_T[tmpidx2])/(s_vector[tmpidx2]-2)))
        }
        # continuous & odd
        tmpidx3 <- setdiff(tmpidx,tmpidx2)
        if(length(tmpidx3)) SepEachPlusmax <- min(c(SepEachPlusmax,SepEach[tmpidx3]*2))
        
        SepUpperBound <- SepEachPlusmax
      }
      SepUpperBound <- min(c(SepUpperBound,sqrt(sum(SepEach^2))))
      if(SepUpperBound<SepBest) IsChance <- 0
      
      
      if(IsChance==1) {
        DimOrder <- order(SepEach)
        wdsOrdered <- SepEach[DimOrder] # d* ordered.
        VectorInConsideration <- cbind(diag(p),wdsOrdered)## add ek, k = 1,2,,..p
        colnames(VectorInConsideration) <- NULL
        VectorListCannot <- matrix(0,0,p) ## no need vector
        
        # Add vectors that must be added.
        VectorToCannotIndex <- order(VectorInConsideration[,p+1])[1] ## consider the vector with the smallest wds
        VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
        while(sqrt(sum((VectorToCannot*wdsOrdered)^2))<SepBest) {  # can be equal
          for(kkk in 1:p) if(VectorToCannot[kkk]>0) break ## find the first neq to 0
          if(kkk-1>0) for(i in 1:(kkk-1)) {
            TempVector <- VectorToCannot
            TempVector[i] <- TempVector[i]+1
            TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered)^2)))
            VectorInConsideration <- rbind(VectorInConsideration,TempVector)
          }
          VectorInConsideration <- VectorInConsideration[-VectorToCannotIndex,,drop=FALSE]## remove vector with the smallest wds
          VectorListCannot <- rbind(VectorListCannot,VectorToCannot) ### put it in cannot set
          VectorToCannotIndex <- order(VectorInConsideration[,p+1])[1] ## consider another vector with the smallest wds
          VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
          if(dim(VectorInConsideration)[1]==0)  { IsChance <- 0;  break; }
          
        }
      }
      
      # Increase nonrep dimensions if necessary.
      if(IsChance==1) {
        is.rep <- rep(1,p)  # assume all rep
        for(j in 1:p) is.rep[j] <- (sum(VectorListCannot[,j])==0) ##check ek
        while(sum(is.rep)>p-p3-1) { # if too much ek, put it into cannot set according to d*
          VectorToCannotIndex <- 1
          VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
          for(kkk in 1:p) if(VectorToCannot[kkk]>0) break
          if(kkk-1>0) for(i in 1:(kkk-1)) {
            TempVector <- VectorToCannot
            TempVector[i] <- TempVector[i]+1
            TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered)^2)))
            VectorInConsideration <- rbind(VectorInConsideration,TempVector)
          }
          VectorInConsideration <- VectorInConsideration[-VectorToCannotIndex,,drop=FALSE]
          VectorListCannot <- rbind(VectorListCannot,VectorToCannot)
          if(dim(VectorInConsideration)[1]==0)  { IsChance <- 0;  break; }
          
          for(j in 1:p) is.rep[j] <- (sum(VectorListCannot[,j])==0)
        }
      }
      
      # Check if the vectors can be arranged (check if u is large enough).
      if(IsChance==1) {
        VectorListMust <- matrix(0,0,p)
        ##if ncosets=0, q(L) \neq q
        # if ncosets=0, there is a vector in VectorListCannot could not be put in any cosets.
        ncosets <- CheckArrange_old(p,p3,VectorListCannot,VectorListMust)
        IsChance <- (ncosets>0)
      }
      
      # Now check subproblems for each p1
      if(IsChance==1) {
        Op1 <- sum(is.rep) ## r, r <  p -p3
        for(p1aim in Op1:0) {
          
          if(p1aim!=Op1) {	# Increase nonrep dimensions by one.
            VectorToCannotIndex <- 1
            VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
            for(kkk in 1:p) if(VectorToCannot[kkk]>0) break## the first >0
            if(kkk-1>0) for(i in 1:(kkk-1)) {##
              TempVector <- VectorToCannot
              TempVector[i] <- TempVector[i]+1
              TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered)^2)))
              VectorInConsideration <- rbind(VectorInConsideration,TempVector)
            }
            VectorInConsideration <- VectorInConsideration[-VectorToCannotIndex,,drop=FALSE]
            VectorListCannot <- rbind(VectorListCannot,VectorToCannot)
            if(dim(VectorInConsideration)[1]==0)  { IsChance <- 0;  break; }
            
            for(j in 1:p) is.rep[j] <- (sum(VectorListCannot[,j])==0)## record ek
          }
          p1 <- sum(is.rep) # r
          
          Sep <- SepUpperBound
          if(p1>0) Sep <- min(c(Sep, sort(SepEach)[p+1-p1] )) #
          
          Ltmps <- Alg7_construcL(p,p3,VectorInConsideration,VectorListCannot,wdsOrdered,is.rep,Sep)
          
          VectorListCannot <- Ltmps[["VectorListCannot"]]
          Sep <- Ltmps[["Sep"]]
          L01alt <- Ltmps[["L01alt"]]
          VectorGeneratorAlt <- Ltmps[["VectorGeneratorAlt"]]
          
          m <- m_exact_fromL01(L01alt,is.rep,s_vector[DimOrder])
          
          if(m>=n){
            # compute sep
            # Algrithm B: find the best separation distance & f given L & s
            # discrete_dimsOrd <- match(discrete_dims,DimOrder)
            
            # order the discrete_dims
            if (!is.null(discrete_dims)){
              dimorder_tmp<-sapply(discrete_dims, function(x){which(DimOrder==x)})
              new_discrete_dims <- sort(dimorder_tmp)
              new_discrete_dims_order <- order(dimorder_tmp)
            }else{
              new_discrete_dims <- discrete_dims
              new_discrete_dims_order <- c()
              
            }
            
            
            tmp_rsults <- Alg4_FindY(new_discrete_dims,L01alt,is.rep,s_vector[DimOrder],weight[DimOrder],Di_min[new_discrete_dims_order],Di_2min[new_discrete_dims_order])
            Sep_B <- tmp_rsults[[1]]
            
            
            if ((Sep-Sep_B)>10^-4) {
              why <- 1
              print('separation distance does not increase')
              #  stop("CHECK WHY!")  # debug
            }
            v_B <- tmp_rsults[[2]]
            sordord <- sort(DimOrder,index.return = T)
            dm <- m-n
            
            if((SepBest-Sep_B)>1e-6) {
              print('separation distance will decrease')
              print(s_vector)
              IsChance <- 0
              # stop("CHECK WHY!")  # debug
              
            }
            
            if(abs(Sep_B-SepBest)<1e-6) {
              changeL <- changeL +1
              dmBest <- dm
              mBest <- m
              p3Best <- p3
              p1Best <- p1
              s_vectorBest <- s_vector
              NCannotBest <- dim(VectorListCannot)[1]
              DimOrderBest <- DimOrder
              is.repBest <- is.rep
              L01altBest <- L01alt
              VectorGeneratorAltBest <- VectorGeneratorAlt
              idbest <- v_B[sordord$ix]
              ### save this results
              L01Best <- matrix(0,dim(L01altBest)[1],p)
              L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
              t_save <- t_save + 1
              save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01=L01Best,L01alt=L01altBest,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest,VectorGeneratorAlt=VectorGeneratorAltBest)
              
              if(p1==0) IsChance <- 0
              
            }
            
            if(Sep_B>=SepBest+1e-6) {
              SepBest <- Sep_B
              dmBest <- m-n
              mBest <- m
              p3Best <- p3
              p1Best <- p1
              s_vectorBest <- s_vector
              NCannotBest <- dim(VectorListCannot)[1]
              DimOrderBest <- DimOrder
              is.repBest <- is.rep
              L01altBest <- L01alt
              VectorGeneratorAltBest <- VectorGeneratorAlt
              idbest <- v_B[sordord$ix]
              changeL <- 0
              
              ### save this results
              L01Best <- matrix(0,dim(L01altBest)[1],p)
              L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
              save_results <- list()
              t_save <-  1
              save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01=L01Best,L01alt=L01altBest,id=idbest, is.rep= is.repBest,DimOrder = DimOrderBest)
              
              if(p1==0) IsChance <- 0
              
              ## update s_limit by rho_B
              v_lim <- rep(0,p)
              if (len_odn>0){
                for (i in 1:len_odn){
                  k <- discrete_dims[i]
                  s_abdn <- which(weight[k]*dplus_max_each[[i]]<SepBest)
                  if (length(s_abdn) !=0){
                    s_lim[k] <- min(s_abdn)
                    tmp <- which(unlist(weight[k]*Di_2min[[i]][[s_lim[k]-1]])>=SepBest)
                    v_lim[k] <- tmp[1]
                  }else{
                    tmp <- which(unlist(weight[k]*Di_2min[[i]][[s_lim[k]-1]])>=SepBest)
                    v_lim[k] <- tmp[1]
                  }
                }
              }
              
              if (len_odn<p){
                for (k in con_dims) {
                  tmp <- floor(2*(weight[k])/SepBest + 1)
                  s_lim[k] <- tmp
                  v_lim[k] <- (weight[k])/(s_lim[k]-1)
                  if (floor(tmp/2)*2!=tmp){
                    tmp2 <- 2*(weight[k]-d1min)/(tmp+1-2)
                    if(tmp2>=SepBest) {
                      s_lim[k] <- tmp+1
                      v_lim[k] <- max(c(weight[k]*d1min,weight[k] - (SepBest)*(s_lim[k]-2)/2))
                    }
                  }
                }
              } # END FOR IF
            }
            
            if(bestNth<N) {
              
              bestNth = bestNth + 1
              L01here <- matrix(0,dim(L01alt)[1],p)
              L01here[,DimOrder[is.rep==0]] <- L01alt
              save_results_bestN[[bestNth]] <-list(Sep=Sep_B,m=m,p3=p3,p1=p1,s_vector=s_vector,L01=L01here,L01alt=L01alt,id=v_B[sordord$ix],is.rep= is.rep,DimOrder=DimOrder)
              sep_bestN[bestNth] <- Sep_B
              sep_bestN_lowest <- min(sep_bestN[is.finite(sep_bestN)])
              lowest_indx <- which(sep_bestN_lowest==sep_bestN)[1]
              # isfull_bestN <- bestNth==N
              
            }else{
              
              if(Sep_B>sep_bestN_lowest){
                L01here <- matrix(0,dim(L01alt)[1],p)
                L01here[,DimOrder[is.rep==0]] <- L01alt
                save_results_bestN[[lowest_indx]] <- list(Sep=Sep_B,m=m,p3=p3,p1=p1,s_vector=s_vector,L01=L01here,L01alt=L01alt,id=v_B[sordord$ix],is.rep= is.rep,DimOrder=DimOrder)
                sep_bestN[lowest_indx] <- Sep_B
                sep_bestN_lowest <- min(sep_bestN[is.finite(sep_bestN)])
                lowest_indx <- which(sep_bestN_lowest==sep_bestN)[1]
                #isfull_bestN <- bestNth==N
              }
            }
            
            
          } # end for m>=n
          
          if(m<n) {
            
            if (sum(is_lmt)==p){tobreak<-1; break;}
            s_avlb_idx <- rev(which(is_lmt==0))
            
            if( s_avlb_idx[1] ==p ) {
              s_vector[p] <- s_vector[p]+1
            }else{
              inc_idx = s_avlb_idx[1]
              s_vector[inc_idx] <- s_vector[inc_idx] +1
              s_vector[(inc_idx+1):p] = 2
              temp =  min(c(ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, s_lim[p] ))  # A quick upper bound of m: prod(ceiling(s_vector[-p]/2)*2)/2^p3
              s_vector[p] = max(c(temp, 2 ))
            }
            
            is_lmt <- s_vector==s_lim
            
            break;  # break r loop}
          }
        }
      }
      
      if(IsChance==0) {
        # change s
        if (0){
          s_vector[p] <- s_lim[p]
          s_available <- rep(1,p)
          s_available[which(s_vector>=s_lim)] <- 0
          
          # if there are no possible s, change to next generator
          if (sum(s_available)==0){tobreak<-1; break;} else{
            for(kk in p:1) if(s_vector[kk]<s_lim[kk])   break;
          }
          
          s_vector[kk] <- s_vector[kk]+1;
          s_vector[(kk+1):p] <- 2
          temp =  min(c(ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, s_lim[p] ))
          s_vector[p] = max(c(temp, 2 ))
          is_lmt <- s_vector==s_lim
        }else{
          
          for(kk in p:1) if(s_vector[kk]>2) break;
          if(kk==1) { tobreak<-1; break; }
          kk <- kk-1
          
          s_available <- rep(1,p)
          s_available[which(s_vector>=s_lim-1e-6)] <- 0
          
          while(kk>0) {   # Check if the kk-1 dimension has chance
            if(!s_available[kk]){kk <- kk-1; next;}
            tmp <- s_vector[kk]+1
            if(kk %in% con_dims){ # continuous dimension
              
              if(floor(tmp/2)*2==tmp){ #  continuous & s_kk +1 is even
                if(weight[kk]/(tmp-2)*2>=SepBest) break
              }else{       #  continuous & s_kk +1 is odd
                if(weight[kk]/(tmp-1)*2>=SepBest) break
              }
              
            }else{ # ordinal
              kkk <- which(discrete_dims==kk)
              if(weight[kk]*dplus_max_each[[kkk]][tmp-1]>=SepBest) break
            }
            kk <- kk-1;
          }
          if(kk==0) { tobreak<-1; break; }
          
          s_vector[kk] <- s_vector[kk]+1  ## kk >= 1
          s_vector[(kk+1):p] <- 2
          temp =  min(c(ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, s_lim[p] ))
          s_vector[p] = max(c(temp, 2 ))
          is_lmt <- s_vector==s_lim
        }
        
      } # end this IsChance
      
    } # end for while
  } # end for p3
  
  
  #if (ReturnAll) return(save_results)
  if (ReturnAll) return(save_results_bestN)
  
  
  ### consider cosets, to find the design with the closest size to n
  ## compare all saved results,
  bestresults <- Alg5_closetoN(save_results,n)
  t_best <- bestresults[["t_best"]]
  id_cosetBest <- bestresults[["id_cosetBest"]]
  n_best <- bestresults[["n_best"]]
  L01altBest <- bestresults[["L01altBest"]]
  u_best <- bestresults[["u_best"]]
  
  
  ## choose the best one
  SepBest <- save_results[[t_best]][["Sep"]]
  DimOrderBest <- save_results[[t_best]][["DimOrder"]]
  s_vectorBest <- save_results[[t_best]][["s_vector"]]
  mBest <- save_results[[t_best]][["m"]]
  idbest <- save_results[[t_best]][["id"]]
  p3Best <- save_results[[t_best]][["p3"]]
  p1Best <- save_results[[t_best]][["p1"]]
  is.repBest <- save_results[[t_best]][["is.rep"]]
  
  ## Generate the best design
  
  L01Best <- matrix(0,dim(L01altBest)[1],p)
  L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
  
  
  Design <- L01Best
  L01BestS <- L01Best
  for(j in p:1) {
    for(jj in 1:p) if(DimOrderBest[jj]==j) break
    if(is.repBest[jj]==0) if(s_vectorBest[j]>3) for(i in 2:floor(s_vectorBest[j]/2)) {
      DesignAdd <- L01Best
      DesignAdd[,j] <- DesignAdd[,j]+(i-1)*2
      Design <- rbind(Design,DesignAdd)
    }
    if(is.repBest[jj]==0) if(floor(s_vectorBest[j]/2)*2 != s_vectorBest[j]) {
      DesignAdd <- L01Best[L01Best[,j]==0,,drop=FALSE]
      DesignAdd[,j] <- DesignAdd[,j]+s_vectorBest[j]-1
      Design <- rbind(Design,DesignAdd)
    }
    if(is.repBest[jj]==1) for(i in 2:s_vectorBest[j]) {
      DesignAdd <- L01Best
      DesignAdd[,j] <- DesignAdd[,j]+(i-1)
      Design <- rbind(Design,DesignAdd)
    }
    L01Best <- Design
  }
  
  h=1
  for(j in discrete_dims){
    D_ordinal <- c(De_good[[h]][[s_vectorBest[j]]][[idbest[j]]])
    for(ss in 0:(s_vectorBest[j]-1)) Design[which(Design[,j]==ss),j] <- D_ordinal[ss+1]
    h=h+1
  }
  for(j in con_dims) Design[,j] <- Design[,j]/(s_vectorBest[j]-1)
  DesignTransformed <- Design
  for(j in 1:p) DesignTransformed[,j] <- DesignTransformed[,j]*weight[j]
  
  
  #return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,p3=p3Best,p1=p1Best,s_vector=s_vectorBest,L01=L01BestS,id=idbest,cosetId =id_cosetBest,u_best=u_best))
  return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,s_vector=s_vectorBest,L01=L01BestS,d_star=idbest,uBest = u_best,ordinal_levels = ordinal_levels))
  
}


######################## Algorithm 9 ########################################
#### To generate the design in p dimensions and at least n points using the weight vector w and Algorithm 9 (p>8 ),
# run: ILMmDMixVarsAlg9(p,n,weight,discrete_dims,ordinal_levels) .


ILMmDMixVarsAlg9 <- function(p,n,discrete_dims,ordinal_levels,weight=rep(1,p),pfrom=8) {
  # THIS algorithm is to obtain the design with the highest separation distance when p>8.
  # p: the number of dimensions
  # n: the number of runs/points
  # weight: prior weight of variables
  # discrete_dims: the index of discrete variables in 1:p
  # ordinal_levels: a list conting the allowable levels of ordinal variables
  # NOTE: for continuous \gamma_k =[0,1], for ordinal \gamma_k = [l_{k,1}=0,...,l_{k,j_k}=1]
  
  if(!p>pfrom|!p==round(p)) stop("p must be an integer greater than pfrom.") 
  if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
  if(!length(weight)==p) stop("weight must be a vector with length p.") 
  
  
  ### load ordinal levels
  len_odn <- length(discrete_dims)
  if (len_odn){
    #load the ordinal list as a matrix
    if(is.character(ordinal_levels)) L_ordinal <- read.csv(file=paste(ordinal_levels,".csv",sep=''),header=T)
    if(is.list(ordinal_levels)){
      my_list = ordinal_levels
      length_diff <- max(sapply(my_list, length))
      my_list_diff_fixed <- lapply(my_list, function(x) {
        c(x, rep(NA, length_diff - length(x)))
      })
      L_ordinal <- data.frame(my_list_diff_fixed)
    }
    if(is.data.frame(ordinal_levels)) L_ordinal <- ordinal_levels
  }else{
    L_ordinal <- matrix(0,0,0)
  }
  
  ### initialize variables
  MAX_row <- dim(L_ordinal)[1]
  slmt_tmp <- rep(MAX_row,len_odn) - sum(is.na(L_ordinal))
  s_lim <- rep(Inf,p)
  s_lim[discrete_dims] <- slmt_tmp
  is_odn <- rep(0,p)
  is_odn[discrete_dims] <-1
  con_dims <-setdiff(c(1:p),discrete_dims)
  
  # run the algorithm A, find all possible combination of ordinal values
  if (len_odn>0){
    ordinal_reultus <- Alg1_DeterFkSk(L_ordinal)
    De_good <- ordinal_reultus[[1]]
    Di_min <- ordinal_reultus[[2]]
    Di_2min <- ordinal_reultus[[3]]
    S_ordinal  <- ordinal_reultus[[4]]
    
  }else{
    De_good <- list()
    Di_min <- list()
    Di_2min <- list()
    S_ordinal  <- list()
    dplus_max_each <- list()
  }
  
  
  Id2max <- pfrom
  WeightOrder <- order(-weight)  # sort the variable by importance
  weightOrdered <- weight[WeightOrder]
  discrete_dims_Ordered <- c()
  discrete_dims2 <- c()  # sort the discrete_dims by importance
  for (i in 1:Id2max) {
    dimtmp <- WeightOrder[i]
    if (dimtmp %in% discrete_dims){
      timpID <- which(discrete_dims==dimtmp)
      discrete_dims_Ordered <- c(discrete_dims_Ordered,timpID)
      discrete_dims2 <- c(discrete_dims2,i)
    }
  }
  
  L_ordinal2  <- L_ordinal[,discrete_dims_Ordered]
  # Run alg 2 to construct designs with p=pfrom and m>=n
  all_Dpart <- ILMmDMixVarsAlg8(Id2max,n,discrete_dims2,L_ordinal2,weightOrdered[1:Id2max],N=20,ReturnAll = "True")
  SepBest <- 0
  dmBest <- Inf
  
  # for loop all designs generated from Alg.2
  for (it in 1:length(all_Dpart)){
    
    Dpart <- all_Dpart[[it]]
    
    s_vectorOrdered <- c(Dpart$s_vector,rep(2,p-Id2max))
    L01 <- Dpart$L01
    for(j in Id2max:1) if(max(L01[,j])==0) { L01Add <- L01; L01Add[,j] <- 1; L01 <- rbind(L01,L01Add); }
    
    for(pp in Id2max:(p-1)){
      L01List <- cbind(L01,rep(-1,dim(L01)[1]),rep(0,dim(L01)[1]))
      # pp+2: the distance from x to 0
      for(j in 1:pp) L01List[,pp+2] <- L01List[,pp+2] + (L01List[,j]*weightOrdered[j]/(s_vectorOrdered[j]-1))^2
      L01List[,pp+2] <- sqrt(L01List[,pp+2])
      VectorListYes <- matrix(0,0,pp)
      VectorOrder <- order(L01List[,pp+2])
      #### L_1:0, L_2:1
      L01List[VectorOrder[1],pp+1] <- 0 # 0_j put in L_1
      L01List[VectorOrder[2],pp+1] <- 1 # put the point with the smallest distance to 0 to L_2
      if(dim(L01List)[1]>2) for(i in 3:dim(L01List)[1]) if(L01List[VectorOrder[i],pp+1]==-1) {
        L01List[VectorOrder[i],pp+1] <- 1
        #If x and y are two vectors in L2, then x + y must be in L1, because otherwise 0j must belong to L2.
        VectorListYes <- AddToGroup_old(VectorListYes,BAdd(L01List[VectorOrder[i],1:pp],L01List[VectorOrder[2],1:pp]))
        for(ii in 1:dim(L01List)[1]) for(iii in 1:dim(VectorListYes)[1])
          if(sum(L01List[ii,1:pp]==VectorListYes[iii,])==pp)
            L01List[ii,pp+1] <- 0
      }
      L01 <- L01List[,1:(pp+1)]
    }
    
    # best L01 alternative
    L01BestSS <- L01
    for(j in 1:p) L01BestSS[,WeightOrder[j]] <- L01[,j] # order by importance
    p1 <- 0
    is.rep <- rep(0,p)  # compute is.rep
    for(j in 1:p) {
      for(i in 2:dim(L01BestSS)[1]) if(sum(L01BestSS[i,-j])==0) {
        p1 <- p1+1
        is.rep[j] <- 1
        L01BestSS <- L01BestSS[L01BestSS[,j]==0,,drop=FALSE]
        break
      }
    }  #L01BestSS doesn't contain e_k, so ncol<=p
    
    s_vector <- rep(0,p)
    for(j in 1:p) s_vector[WeightOrder[j]] <- s_vectorOrdered[j]   # order by importance
    
    # determine d^* index
    v_Best <- rep(NA,p)
    idbest8 <- Dpart[["id"]] # v, the index of d^*
    for(j in 1:Id2max)  v_Best[WeightOrder[j]] <- idbest8[j]  # order by importance
    checkna <- which(is.na(v_Best))
    
    for (i in checkna){
      if (i %in% discrete_dims) v_Best[i] <- 1
      if (i %in% con_dims) v_Best[i] <- weight[i]/(s_vector[i]-1)
    }
    
    # determine y
    SepEach <- rep(0,p)  # didn't change the order
    Sep2Each <- rep(0,p)
    for( i in 1:p){
      if (i %in% discrete_dims){
        k <- which(discrete_dims==i)
        SepEach[i] <- weight[i]*Di_min[[k]][[s_vector[i]-1]][v_Best[i]]
        Sep2Each[i] <- weight[i]*Di_2min[[k]][[s_vector[i]-1]][v_Best[i]]
      }
      
      if (i %in% con_dims) {
        SepEach[i] <- v_Best[i]
        Sep2Each[i] <- weight[i]
        if (s_vector[i] > 2) Sep2Each[i] <- 2*(weight[i]-v_Best[i])/(s_vector[i]-2)
      }
    }
    
    
    # compute the distance of x in {0,1}^p to 0
    SepEachOrdered <- SepEach[WeightOrder]
    
    L01List[,p+1] <- 0
    for(j in 1:p) L01List[,p+1] <- L01List[,p+1] + (L01List[,j]*SepEachOrdered[j])^2
    L01List[,p+1] <- sqrt(L01List[,p+1])  ## p+1: the distance of x in {0,1}^p to 0
    
    # compute separation distance
    Sep_B <- min(L01List[-1,p+1])
    if (max(s_vector)>2) Sep_B <- min(c(Sep_B,Sep2Each[s_vector>2]))
    if (sum(s_vector>2)>0 && (Sep_B==min(Sep2Each[s_vector>2])) ){
      L01alt <- L01BestSS
      tmp_rsults <- Alg4_FindY(discrete_dims,L01alt,is.rep,s_vector,weight,Di_min,Di_2min)
      Sep1 <- tmp_rsults[[1]]
      v_Best <- tmp_rsults[[2]]
      if (Sep_B>=Sep1+1e-6){
        print("Sep becomes lower!!!! CHECK!!!!")
        # stop("CHECK WHY!") # debug
      }
      Sep_B <- Sep1
    }
    idbest <- v_Best
    
    m <- Dpart$m
    L01tmp <- matrix(0,dim(L01)[1],p) # ncol = p
    for(j in 1:p) L01tmp[,WeightOrder[j]] <- L01[,j]
    dm <- m-n
    
    
    if(abs(Sep_B-SepBest)<1e-6) {
      SepBest <- Sep_B
      dmBest <- dm
      mBest <- m
      p1Best <- p1
      s_vectorBest <- s_vector
      is.repBest <- is.rep
      L01Best <- L01tmp
      L01altBest <- L01BestSS
      idbest <- v_Best
      ### save this results
      t_save <- t_save + 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p1=p1Best,s_vector=s_vectorBest,L01=L01Best,L01alt=L01altBest,id=idbest, is.rep= is.repBest)
    }
    
    
    if(Sep_B>=SepBest+1e-6) {
      SepBest <- Sep_B
      dmBest <- dm
      mBest <- m
      p1Best <- p1
      s_vectorBest <- s_vector
      is.repBest <- is.rep
      L01Best <- L01tmp
      L01altBest <- L01BestSS  #L01BestS
      idbest <- v_Best
      ### save this results
      save_results <- list()
      t_save <- 1
      save_results[[t_save]] <- list(Sep=SepBest,m=mBest,p1=p1Best,s_vector=s_vectorBest,L01=L01Best,L01alt=L01altBest,id=idbest, is.rep= is.repBest)
      
    }
  }
  
  
  
  
  ### consider cosets
  ## compare all saved results,
  
  bestresults <- Alg5_closetoN(save_results,n)
  t_best <- bestresults[["t_best"]]
  id_cosetBest <- bestresults[["id_cosetBest"]]
  n_best <- bestresults[["n_best"]]
  L01altBest <- bestresults[["L01altBest"]]
  u_best <- bestresults[["u_best"]]
  
  
  p1Best <- save_results[[t_best]][["p1"]]
  p3Best <- save_results[[t_best]][["p3"]]
  is.repBest <- save_results[[t_best]][["is.rep"]]
  
  ## Generate the best design
  L01Best <- matrix(0,dim(L01altBest)[1],p)
  L01Best[,is.repBest==0]=L01altBest
  
  
  Design <- L01Best
  for(j in p:1) {
    if(s_vectorBest[j]>3) for(i in 2:floor(s_vectorBest[j]/2)) {
      DesignAdd <- L01Best
      DesignAdd[,j] <- DesignAdd[,j]+(i-1)*2
      Design <- rbind(Design,DesignAdd)
    }
    if(floor(s_vectorBest[j]/2)*2 != s_vectorBest[j]) {
      DesignAdd <- L01Best[L01Best[,j]==0,,drop=FALSE]
      DesignAdd[,j] <- DesignAdd[,j]+s_vectorBest[j]-1
      Design <- rbind(Design,DesignAdd)
    }
    L01Best <- Design
  }
  h=1
  for(j in discrete_dims){
    D_ordinal <- c(De_good[[h]][[s_vectorBest[j]]][[idbest[j]]])
    for(ss in 0:(s_vectorBest[j]-1)) Design[which(Design[,j]==ss),j] <- D_ordinal[ss+1]
    h=h+1
  }
  for(j in con_dims) Design[,j] <- Design[,j]/(s_vectorBest[j]-1)
  DesignTransformed <- Design
  for(j in 1:p) DesignTransformed[,j] <- DesignTransformed[,j]*weight[j]
  
  
  #return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,p1=p1Best,p3=p3Best,s_vector=s_vectorBest,L01=L01altBest,id=idbest,cosetId =id_cosetBest,u_best=u_best))
  return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,s_vector=s_vectorBest,L01=L01altBest,d_star=idbest,uBest = u_best,ordinal_levels = ordinal_levels))
  
}


InterleavedMaximinDMixVars <- function(p,n,discrete_dims,ordinal_levels,weight = rep(1, p)) {
  
  if(!p>=2|!p==round(p)) stop("p must be an integer greater than one and no greater than five.")
  if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.")
  if(!length(weight)==p) stop("weight must be a vector with length p.")
  
  if(p<=5) return(ILMmDMixVarsAlg6(p,n,discrete_dims,ordinal_levels,weight = rep(1, p)))
  if(p<=8) return(ILMmDMixVarsAlg8(p,n,discrete_dims,ordinal_levels,weight = rep(1, p), N = 1, ReturnAll = 0))
  return(ILMmDMixVarsAlg9(p,n,discrete_dims,ordinal_levels,weight = rep(1, p)))
}




############################### codes from other package





tnp1 <- function(lsss,coefC) {
  oneortwo = (lsss>2)
  counts = rep(0,dim(coefC)[2])
  if(sum(oneortwo)==0)  counts = coefC[ sum(c(2-lsss)*2^((length(lsss)-1):0))+1 ,]
  if(sum(oneortwo)!=0)  {
    for(j in 1:length(lsss)) if(lsss[j]>2) break
    thedim = j
    thetwos = floor((lsss[j]-1)/2)
    lsss1 = lsss
    lsss1[thedim]=2
    lsss2 = lsss
    lsss2[thedim]=lsss[thedim]-thetwos*2
    counts = thetwos*tnp1(lsss1,coefC) + tnp1(lsss2,coefC)
  }
  
  counts
}



