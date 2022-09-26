gradient.new = function(w,A,Y){
  n = length(w)
  rowA = nrow(A)
  R0 = matrix(exp(w),n,rowA)
  R0 = t(R0)
  RR0 = rowSums(A*R0)
  RR1 = R0/RR0
  RR1 = RR1*A
  G1 = (YY0-RR1)*A
  G = colSums(G1)
  return(G)
}
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
MLE.new = function(A,Y){
  w0 = rnorm(n); w0 = w0 - mean(w0)
  Gw = gradient.new(w0,A,Y)
  V = max(abs(Gw))
  w = w0
  k0 = 1
  while(V > 0.00001 && k0 <= 10000){
    w = w + 0.01*Gw
    print(mean(w))
    ##w = w - mean(w)
    Gw = gradient.new(w,A,Y)
    V = max(abs(Gw))
    print(V)
    k0 = k0 + 1
  }
  return(w)
}
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
variance.new = function(m,w,A){
  n = length(w)
  rowA = nrow(A)
  R0 = matrix(exp(w),n,rowA)
  R0 = t(R0)
  RR0 = rowSums(A*R0)
  RR1 = R0/RR0
  RR1 = RR1*A
  Idx.m = (A[,m]==1)
  tmp.sub.A.m = A[Idx.m,m]
  tmp.sub.R.m = RR1[Idx.m,m]
  V1.tmp = tmp.sub.R.m*(1-tmp.sub.R.m)*tmp.sub.A.m
  VV = sum(V1.tmp)
  return(VV)
}
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
F.new = function(m,w,A,Y.A,L){
  n = length(w)
  rowA = nrow(A)
  R0 = matrix(exp(w),n,rowA)
  R0 = t(R0)
  RR0 = rowSums(A*R0)
  RR1 = R0/RR0
  RR1 = RR1*A
  Idx.m = (A[,m]==1)
  tmp.sub.R.m = RR1[Idx.m,m]
  Y.A.tmp = Y.A[Idx.m,m,]-tmp.sub.R.m
  FF = sum(Y.A.tmp)
  return(FF)
}
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
xx1 = read.csv('~/Desktop/jester_1.csv')
xx1 = xx1[,-1]
xx2 = read.csv('~/Desktop/jester_2.csv')
xx2 = xx2[,-1]
xxyy = rbind(xx1,xx2)   ## dim 14116*100
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
set.seed(2022)                       
MM0 = 2; n = 100; p = 0.05; L = 80; B = 1000
seq.C = 1#seq(0,1,length=10)
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
edge.NN = choose(n,MM0)
graph.NN = rbinom(1,edge.NN,p)
AB = matrix(0,graph.NN+1000,n)
AB[,1:MM0] = 1
AB = apply(AB,1,sample)
AB = t(AB)
A00.tmp = unique(AB)
A = A00.tmp[1:graph.NN,]
YY0 = matrix(0,graph.NN,n)
YYY0 = array(0,c(graph.NN,n,L))
for(o in 1:graph.NN){
  tmp.A.v = (A[o,]==1)
  IdxL = sample(1:14116,L)
  for(oo in 1:L){
    Y.value = as.numeric(xxyy[IdxL[oo],tmp.A.v])
    Y.max = which.max(Y.value)
    YY0.tmp = numeric(MM0)
    YY0.tmp[Y.max] = 1
    YYY0[o,tmp.A.v,oo] = YY0.tmp
    YY0[o,tmp.A.v] = YY0[o,tmp.A.v] + YY0.tmp/L
  }
}
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
w.est = MLE.new(A,YY0)                                                         ## MLE estimator
w.sort = sort(w.est)
MM = 15                                                                        ## MM
w.tmp = tail(w.sort,MM)
w.idx = numeric(MM)                                                            ## w.idx is the index of largest MM values of theta
for(pp in 1:MM){
  w.idx[pp] = which(w.est == w.tmp[pp])
}
w.idx
w.idx_new=c(10,30,50,70,90)
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
var.V = numeric(n)                                                             
F.V = matrix(0,n,L)
for(kk in c(1:n)){
  var.V[kk] = variance.new(kk,w.est,A)
  F.V[kk,] = F.new(kk,w.est,A,YYY0,L)
}                                                                              ## compute f and g
FG.V = F.V/var.V
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
B.M = matrix(rnorm(L*B),L,B)
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
FF = FG.V[1,] - t(FG.V[-1,])
SDD = sqrt(1/var.V[1]+1/var.V[-1])
FF.N = t(FF)/SDD
for(mm in 2:n){
  FF.tmp = FG.V[mm,] - t(FG.V[-mm,])
  SDD.tmp = sqrt(1/var.V[mm]+1/var.V[-mm])
  FF.N.tmp = t(FF.tmp)/SDD.tmp
  FF.N = rbind(FF.N,FF.N.tmp)
}
FF.B  = FF.N%*%B.M
cut.FF = apply(FF.B,2,max,na.rm=TRUE)                                                      ## uniform normalized one-sided CI
cut.v.FF = quantile(cut.FF,0.95,na.rm=TRUE)/sqrt(L)                                        ## R.one.uniform 
R.one.uniform = numeric(n)
for(mmm in 1:n){
  w.V.uniform = w.est[mmm]-w.est[-mmm]
  SDD.tmp.uniform = sqrt(1/var.V[mm]+1/var.V[-mm])
  w.V.N.uniform = w.V.uniform*sqrt(L)/SDD.tmp.uniform
  R.one.uniform[mmm] = 1 + sum(1*(-w.V.N.uniform>cut.v.FF))
}
R.one.uniform
write.csv(R.one.uniform,"R.one.uniform.csv")
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
FFF = FG.V[1,] - t(FG.V[-1,])
FFF.N = t(FFF)
for(mm in 2:n){
  FFF.tmp = FG.V[mm,] - t(FG.V[-mm,])
  FFF.N.tmp = t(FFF.tmp)
  FFF.N = rbind(FFF.N,FFF.N.tmp)
}
FFF.B  = FFF.N%*%B.M
cut.FFF = apply(FFF.B,2,max,na.rm=TRUE)                                                      ## uniform one-sided CI
cut.v.FFF = quantile(cut.FFF,0.95,na.rm=TRUE)/sqrt(L)                                        ## RR.one.uniform
RR.one.uniform = numeric(n)
for(mmm in 1:n){
  w.V.u = w.est[mmm]-w.est[-mmm]
  w.V.N.u = w.V.u*sqrt(L)
  RR.one.uniform[mmm] = 1 + sum(1*(-w.V.N.u>cut.v.FFF))
}
RR.one.uniform
write.csv(RR.one.uniform,"RR.one.uniform.csv")
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
R.left_vec<-c()
R.right_vec<-c()
R.one_vec<-c()
R.left.C_vec<-c()
R.right.C_vec<-c()
R.left.N_vec<-c()
R.right.N_vec<-c()
for (i in c(1:length(w.idx_new))){
  #m = 53  ## item of interest
  m=w.idx_new[i]
  ##---------------------------------------------------------------------------------##
  ##---------------------------------------------------------------------------------##
  F.M = FG.V[m,] - t(FG.V[-m,])
  sd.M = sqrt(1/var.V[m]+1/var.V[-m])
  F.M.N = t(F.M)/sd.M
  B.M.P = F.M.N%*%B.M
  cut.tmp = apply(abs(B.M.P),2,max,na.rm=TRUE)
  cut.v = quantile(cut.tmp,0.95,na.rm=TRUE)/sqrt(L)
  cut.tmp.one = apply(B.M.P,2,max,na.rm=TRUE)
  cut.v.one = quantile(cut.tmp.one,0.95,na.rm=TRUE)/sqrt(L)
  w.V = w.est[m]-w.est[-m]
  w.V.N = w.V*sqrt(L)/sd.M
  R.left = 1 + sum(1*(-w.V.N>cut.v))                                             ## [R.left,R.right] CI for item m
  R.right = n - sum(1*(w.V.N>cut.v))
  R.length = R.right-R.left
  R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
  R.one = 1 + sum(1*(-w.V.N>cut.v.one))                                          ## [R.one,n] one-sided CI for item m
  cutZ = qnorm(0.975)
  C.len = length(seq.C)
  cut.C = numeric(C.len)
  R.left.C = numeric(C.len); R.right.C = numeric(C.len)
  R.CI.C = numeric(C.len)
  R.length.C = numeric(C.len)
  R.C.I = numeric(C.len)
  R.left.N = numeric(C.len); R.right.N = numeric(C.len)
  R.CI.N = numeric(C.len)
  R.length.N = numeric(C.len)
  R.N.I = numeric(C.len)
  for(uu in 1:C.len){
    tmp.C = seq.C[uu]
    sd.M.C = cutZ/sqrt(var.V[m])/(1+tmp.C)/sqrt(2*log(n))+1/sqrt(var.V[-m])
    F.M.N.C = t(F.M)/sd.M.C
    B.M.P.C = F.M.N.C%*%B.M
    cut.tmp.C = apply(abs(B.M.P.C),2,max,na.rm=TRUE)
    cut.C[uu] = quantile(cut.tmp.C,0.95,na.rm=TRUE)/sqrt(L)
    w.V.N.C = w.V*sqrt(L)/sd.M.C
    R.left.C[uu] = 1 + sum(1*(-w.V.N.C>cut.C[uu]))
    R.right.C[uu] = n - sum(1*(w.V.N.C>cut.C[uu]))
    R.length.C[uu] = R.right.C[uu]-R.left.C[uu]
    tmp.N.cut = (1+tmp.C)*sqrt(2*log(n))
    R.left.N[uu] = 1 + sum(1*(-w.V.N.C>tmp.N.cut))
    R.right.N[uu] = n - sum(1*(w.V.N.C>tmp.N.cut))
    R.length.N[uu] = R.right.N[uu]-R.left.N[uu]
    R.CI.C[uu] = sum(1*(m<R.left.C[uu]))+sum(1*(m>R.right.C[uu]))
    R.CI.N[uu] = sum(1*(m<R.left.N[uu]))+sum(1*(m>R.right.N[uu]))
  }
  #####################################################################################
  R.left_vec<-c(R.left_vec, R.left)                                                                        ## [R.left,R.right] is our two-sided CI      
  write.csv(R.left_vec,"R.left_vec.csv")
  R.right_vec<-c(R.right_vec,R.right)   
  write.csv(R.right_vec,"R.right_vec.csv")
  R.one_vec<-c(R.one_vec,R.one)                                                                                ## [R.left,n] is our one-sided CI
  write.csv(R.one_vec,"R.one_vec.csv.csv")
  R.left.C_vec<-c(R.left.C_vec,R.left.C)  
  write.csv(R.left.C_vec,"R.left.C_vec.csv")
  R.right.C_vec<-c(R.right.C_vec,R.right.C) 
  write.csv(R.right.C_vec,"R.right.C_vec.csv")
  R.left.N_vec<-c(R.left.N_vec,R.left.N) 
  write.csv(R.left.N_vec,"R.left.N_vec.csv")                                                                ## [R.left.N,R.right.N] is Chao Gao CI       
  R.right.N_vec<-c(R.right.N_vec,R.right.N) 
  write.csv(R.right.N_vec,"R.right.N_vec.csv")}
  #####################################################################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  