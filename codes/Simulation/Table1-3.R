gradient = function(w,A,Y1,Y2,Y3){
	n = length(w)
	R1 = array(exp(w),c(n,n,n))
	R2 = aperm(R1,c(2,3,1))
	R3 = aperm(R1,c(3,1,2))
	RR1 = R1/(R1+R2+R3)
	RR2 = R3/(R1+R2+R3)
	RR3 = R2/(R1+R2+R3)
	G1 = A*(Y1-RR1); G2 = A*(Y2-RR2); G3 = A*(Y3-RR3)
	G1.tmp = aperm(G1,c(2,3,1))
	G1.gradient = colSums(colSums(G1.tmp))
	G2.tmp = aperm(G2,c(3,1,2))
	G2.gradient = colSums(colSums(G2.tmp))
	G3.tmp = G3
	G3.gradient = colSums(colSums(G3.tmp))
	G = G1.gradient+G2.gradient+G3.gradient 
	return(G)
}
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
MLE = function(A,Y1,Y2,Y3){
	w0 = rnorm(n); w0 = w0 - mean(w0)
	Gw = gradient(w0,A,Y1,Y2,Y3)
	V = max(abs(Gw))
	w = w0
	k0 = 1
	while(V > 0.0000001 && k0 <= 10000){
		w = w + 0.01*Gw
		##w = w - mean(w)
		Gw = gradient(w,A,Y1,Y2,Y3)
		V = max(abs(Gw))
		##print(V)
		k0 = k0 + 1
	}
	return(w)
}
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
variance = function(m,w,A){
	n = length(w)
	R1 = array(exp(w),c(n,n,n))
	R2 = aperm(R1,c(2,3,1))
	R3 = aperm(R1,c(3,1,2))
	RR1 = R1/(R1+R2+R3)
	RR2 = R3/(R1+R2+R3)
	RR3 = R2/(R1+R2+R3)
	V1.tmp = A*RR1*(1-RR1); V2.tmp = A*RR2*(1-RR2); V3.tmp = A*RR3*(1-RR3)
	VV = sum(V1.tmp[m,,]+V2.tmp[,m,]+V3.tmp[,,m])
	return(VV)
}
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
F = function(m,w,A,YY1,YY2,YY3,L){
	FF = numeric(L)
	n = length(w)
	R1 = array(exp(w),c(n,n,n))
	R2 = aperm(R1,c(2,3,1))
	R3 = aperm(R1,c(3,1,2))
	RR1 = R1/(R1+R2+R3)
	RR2 = R3/(R1+R2+R3)
	RR3 = R2/(R1+R2+R3)
	for(ffff in 1:L){
		GG1 = A*(YY1[,,,ffff]-RR1); GG2 = A*(YY2[,,,ffff]-RR2); GG3 = A*(YY3[,,,ffff]-RR3)
		FF[ffff] = sum(GG3[,,m])+sum(GG2[,m,])+sum(GG1[m,,])
	}
	return(FF)
}
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
GA.CI = function(n,m,A,Y1,Y2,Y3,YY1,YY2,YY3,L,B,W,seq.C){
	w.est = MLE(A,Y1,Y2,Y3)
	var.V = numeric(n)
	F.V = matrix(0,n,L)
	for(kk in 1:n){
		var.V[kk] = variance(kk,w.est,A)
		F.V[kk,] = F(kk,w.est,A,YY1,YY2,YY3,L)
	}
	FG.V = F.V/var.V
	F.M = FG.V[m,] - t(FG.V[-m,])
	sd.M = sqrt(1/var.V[m]+1/var.V[-m])
	F.M.N = t(F.M)/sd.M
	B.M = matrix(rnorm(L*B),L,B)
	B.M.P = F.M.N%*%B.M
	cut.tmp = apply(abs(B.M.P),2,max)
	cut.v = quantile(cut.tmp,0.95)/sqrt(L)
#####################################################################################
#####################################################################################
	w.V = w.est[m]-w.est[-m]
	w.V.N = w.V*sqrt(L)/sd.M
	R.left = 1 + sum(1*(-w.V.N>cut.v))
	R.right = n - sum(1*(w.V.N>cut.v))
	R.length = R.right-R.left
	R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
#####################################################################################
#####################################################################################
	ww.V = w.est[m]-w.est[-m]-(W[m]-W[-m])
	stat.V = max(abs(ww.V/sd.M))*sqrt(L)
	vv = 1*(stat.V>cut.v)
#####################################################################################
#####################################################################################
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
		##sd.M.C = sqrt(cutZ/var.V[m]/(1+tmp.C)/sqrt(2*log(n))+1/var.V[-m])
		F.M.N.C = t(F.M)/sd.M.C
		B.M.P.C = F.M.N.C%*%B.M
		cut.tmp.C = apply(abs(B.M.P.C),2,max)
		cut.C[uu] = quantile(cut.tmp.C,0.95)/sqrt(L)
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
		stat.V.C = max(abs(ww.V/sd.M.C))*sqrt(L)
		R.C.I[uu] = 1*(stat.V.C>cut.C[uu])
		R.N.I[uu] = 1*(stat.V.C>tmp.N.cut)
	}
#####################################################################################
#####################################################################################
	result = matrix(0,3,(1+2*C.len))
	result[1,] = c(R.CI,R.CI.C,R.CI.N)
	result[2,] = c(R.length,R.length.C,R.length.N)
	result[3,] = c(vv,R.C.I,R.N.I)
	return(result)
}
#####################################################################################





		














































































#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
n = 80; p = 0.05; L = 60; B = 1000; NN = 500; m = 10                                  ## item m 
seq.C = seq(0,1,length=10)                                                           ## sequence of C
C.len = length(seq.C)
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
W.est.M = matrix(0,n,NN)                          ## MLE estimator 
R.one.uniform.M = matrix(0,n,NN)                  ## normalized one sided CI for all items 
R.one.uniform.V = numeric(NN)                     ## normalized one sided CI length for all items    
RR.one.uniform.M = matrix(0,n,NN)                 ## not normalized 
RR.one.uniform.V = numeric(NN)
R.two.V = numeric(NN)                             ## empirical size of two sided CI (theta)
R.two.length.V = numeric(NN)                      ## two sided CI length 
R.CI.V = numeric(NN)                              ## empirical size of two sided CI (rank)
R.one.V = numeric(NN)                             ## empirical size of one sided CI (theta)   
##################################################################################### ## compare with Chao Gao
R.CI.C.M = matrix(0,C.len,NN)                     ## empirical size of two sided CI (rank) based on GA 
R.length.C.M = matrix(0,C.len,NN)                 ## length CI based on GA 
R.C.I.M = matrix(0,C.len,NN)                      ## empirical size CI (theta)
R.CI.N.M = matrix(0,C.len,NN)                     ## empirical size of two sided CI (rank) Chao Gao  
R.length.N.M = matrix(0,C.len,NN)                 ## length CI Chao Gao 
R.N.I.M = matrix(0,C.len,NN)                      ## empirical size CI (theta) Chao gao 
#####################################################################################
#####################################################################################
K = 10
seq.m = (K-5):(K+5)                            ## item index 
#####################################################################################
len.m = length(seq.m)
K.test.M = matrix(0,len.m,NN)                     ## empirical size of one-sided CI (theta)
K.test.rank.M = matrix(0,len.m,NN)                ## empirical size of one-sided CI (rank)
#####################################################################################
#####################################################################################
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
KK.seq = c(5,10,15)                                                                ##--------## sequence of KK
KK.len = length(KK.seq)                                                          
KK.set.R = matrix(0,KK.len,NN)                                                     ##--------## set cardinality uniform normalized CI 
KK.power.R = matrix(0,KK.len,NN)                                                   ##--------## uniform normalized CI empirical power 
KK.set.RR = matrix(0,KK.len,NN)                                                    ##--------## set cardinality uniform CI 
KK.power.RR = matrix(0,KK.len,NN)                                                  ##--------## uniform CI empirical power
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
#####################################################################################
#####################################################################################
#####################################################################################
for(hh in 1:NN){
print(hh)
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
YY1 = array(0,c(n,n,n,L))
YY2 = array(0,c(n,n,n,L))
YY3 = array(0,c(n,n,n,L))
##W = runif(n,2,4); W = W - mean(W)  ## true theta
W = seq(4,2,length=n); W = W - mean(W)                                               ## true theta
for(o in 1:(n-2)){
	for(oo in (o+1):(n-1)){
		for(ooo in (oo+1):n){
			A[o,oo,ooo] = rbinom(1,1,p)
			p1.tmp = exp(W[o])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p2.tmp = exp(W[oo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p3.tmp = exp(W[ooo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p.tmp = c(p1.tmp,p2.tmp,p3.tmp)
			for(oooo in 1:L){		
				Y.tmp = rmultinom(1,1,p.tmp)
				YY1[o,oo,ooo,oooo] = Y.tmp[1]
				YY2[o,oo,ooo,oooo] = Y.tmp[2]
				YY3[o,oo,ooo,oooo] = Y.tmp[3]
				Y1[o,oo,ooo] = Y1[o,oo,ooo] + YY1[o,oo,ooo,oooo]/L
				Y2[o,oo,ooo] = Y2[o,oo,ooo] + YY2[o,oo,ooo,oooo]/L          ##rm() and gc() to save as sparse tensors 
				Y3[o,oo,ooo] = Y3[o,oo,ooo] + YY3[o,oo,ooo,oooo]/L
			}
		}
	}
}
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
	w.est = MLE(A,Y1,Y2,Y3)                                                        ## MLE estimator
	W.est.M[,hh] = w.est
#####################################################################################
	var.V = numeric(n)                                                             
	F.V = matrix(0,n,L)
	for(kk in 1:n){
		var.V[kk] = variance(kk,w.est,A)
		F.V[kk,] = F(kk,w.est,A,YY1,YY2,YY3,L)
	}                                                                              ## compute f and g
	FG.V = F.V/var.V
#####################################################################################
	B.M = matrix(rnorm(L*B),L,B)                                                   ## Bootstrap matrix
#####################################################################################
#####################################################################################
#####################################################################################
	FF = FG.V[1,] - t(FG.V[-1,])
	SDD = sqrt(1/var.V[1]+1/var.V[-1])
	FF.N = t(FF)/SDD
	ww.V.R = w.est[1]-w.est[-1]-(W[1]-W[-1])
	ww.V.R.N = ww.V.R/SDD
	for(mm in 2:n){
		FF.tmp = FG.V[mm,] - t(FG.V[-mm,])
		SDD.tmp = sqrt(1/var.V[mm]+1/var.V[-mm])
		FF.N.tmp = t(FF.tmp)/SDD.tmp
		FF.N = rbind(FF.N,FF.N.tmp)
		ww.V.R.tmp = w.est[mm]-w.est[-mm]-(W[mm]-W[-mm])
		ww.V.R.N.tmp = ww.V.R.tmp/SDD.tmp
		ww.V.R.N = c(ww.V.R.N,ww.V.R.N.tmp)
	}
	FF.B  = FF.N%*%B.M*(-1)
	cut.FF = apply(FF.B,2,max)                                                      ## uniform normalized one-sided CI
	cut.v.FF = quantile(cut.FF,0.95)/sqrt(L)                                        ## R.one.uniform 
	stat.V.R = max(-ww.V.R.N)*sqrt(L)
	vv.V.R = 1*(stat.V.R>cut.v.FF)
	R.one.uniform.V[hh] = vv.V.R
	R.one.uniform = numeric(n)
	for(mmm in 1:n){
		w.V.uniform = w.est[mmm]-w.est[-mmm]
		SDD.tmp.uniform = sqrt(1/var.V[mmm]+1/var.V[-mmm])
		w.V.N.uniform = w.V.uniform*sqrt(L)/SDD.tmp.uniform
		R.one.uniform[mmm] = 1 + sum(1*(-w.V.N.uniform>cut.v.FF))
	}
	R.one.uniform.M[,hh] = R.one.uniform
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
	FFF = FG.V[1,] - t(FG.V[-1,])
	FFF.N = t(FFF)
	ww.V.RR = w.est[1]-w.est[-1]-(W[1]-W[-1])
	for(mm in 2:n){
		FFF.tmp = FG.V[mm,] - t(FG.V[-mm,])
		FFF.N.tmp = t(FFF.tmp)
		FFF.N = rbind(FFF.N,FFF.N.tmp)
		ww.V.RR.tmp = w.est[mm]-w.est[-mm]-(W[mm]-W[-mm])
		ww.V.RR = c(ww.V.RR,ww.V.RR.tmp)
	}
	FFF.B  = FFF.N%*%B.M*(-1)
	cut.FFF = apply(FFF.B,2,max)                                                      ## uniform one-sided CI
	cut.v.FFF = quantile(cut.FFF,0.95)/sqrt(L)                                        ## RR.one.uniform
	stat.V.RR = max(-ww.V.RR)*sqrt(L)
	vv.V.RR = 1*(stat.V.RR>cut.v.FFF)
	RR.one.uniform.V[hh] = vv.V.RR
	RR.one.uniform = numeric(n)
	for(mmm in 1:n){
		w.V.u = w.est[mmm]-w.est[-mmm]
		w.V.N.u = w.V.u*sqrt(L)
		RR.one.uniform[mmm] = 1 + sum(1*(-w.V.N.u>cut.v.FFF))
	}
	RR.one.uniform.M[,hh] = RR.one.uniform
#####################################################################################
#####################################################################################
#####################################################################################
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
for(kkkk in 1:KK.len){
	KK.tmp = KK.seq[kkkk]
	KK.vector = 1:KK.tmp
	KK.set.R[kkkk,hh] = n - sum(1*(R.one.uniform>KK.tmp))
	KK.set.RR[kkkk,hh] = n - sum(1*(RR.one.uniform>KK.tmp))
	KK.seq.tmp.R = c(1*(R.one.uniform[KK.vector]>KK.tmp))
	KK.seq.tmp.RR = c(1*(RR.one.uniform[KK.vector]>KK.tmp))
	KK.power.tmp.R = 1*(sum(KK.seq.tmp.R)>0)
	KK.power.tmp.RR = 1*(sum(KK.seq.tmp.RR)>0)
	KK.power.R[kkkk,hh] = KK.power.tmp.R
	KK.power.RR[kkkk,hh] = KK.power.tmp.RR
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
#####################################################################################
#####################################################################################
	F.M = FG.V[m,] - t(FG.V[-m,])
	sd.M = sqrt(1/var.V[m]+1/var.V[-m])
	F.M.N = t(F.M)/sd.M
	B.M.P = F.M.N%*%B.M
	cut.tmp = apply(abs(B.M.P),2,max)
	cut.v = quantile(cut.tmp,0.95)/sqrt(L)
#####################################################################################
#####################################################################################
	cut.tmp.one = apply(-B.M.P,2,max)
	cut.v.one = quantile(cut.tmp.one,0.95)/sqrt(L)
	w.V = w.est[m]-w.est[-m]
	w.V.N = w.V*sqrt(L)/sd.M
#####################################################################################
#####################################################################################
	ww.V = w.est[m]-w.est[-m]-(W[m]-W[-m])
	stat.V = max(abs(ww.V/sd.M))*sqrt(L)
	vv = 1*(stat.V>cut.v)
	R.two.V[hh] = vv
	stat.V.one = max(-ww.V/sd.M)*sqrt(L)
	vv.one = 1*(stat.V.one>cut.v.one)
	R.one.V[hh] = vv.one
#####################################################################################
#####################################################################################
	R.left = 1 + sum(1*(-w.V.N>cut.v))                                             ## [R.left,R.right] CI for item m
	R.right = n - sum(1*(w.V.N>cut.v))
	R.length = R.right-R.left
	R.two.length.V[hh] = R.length
	R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
	R.CI.V[hh] = R.CI
	R.one = 1 + sum(1*(-w.V.N>cut.v.one))                                          ## [R.one,n] one-sided CI for item m
#####################################################################################
#####################################################################################
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
		cut.tmp.C = apply(abs(B.M.P.C),2,max)
		cut.C[uu] = quantile(cut.tmp.C,0.95)/sqrt(L)
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
		stat.V.C = max(abs(ww.V/sd.M.C))*sqrt(L)
		print(stat.V.C)
		R.C.I[uu] = 1*(stat.V.C>cut.C[uu])
		R.N.I[uu] = 1*(stat.V.C>tmp.N.cut)
	}
	R.CI.C.M[,hh] = R.CI.C
	R.length.C.M[,hh] = R.length.C
	R.C.I.M[,hh] = R.C.I
	R.CI.N.M[,hh] = R.CI.N
	R.length.N.M[,hh] = R.length.N
	R.N.I.M[,hh] = R.N.I
#####################################################################################
	for(hhhh in 1:len.m){
		tmp.mm = seq.m[hhhh]
		F.M = FG.V[tmp.mm,] - t(FG.V[-tmp.mm,])
		sd.M = sqrt(1/var.V[tmp.mm]+1/var.V[-tmp.mm])
		F.M.N = t(F.M)/sd.M
		B.M.P = F.M.N%*%B.M
#####################################################################################
#####################################################################################
		cut.tmp.one = apply(-B.M.P,2,max)
		cut.v.one = quantile(cut.tmp.one,0.95)/sqrt(L)
		w.V = w.est[tmp.mm]-w.est[-tmp.mm]
		w.V.N = w.V*sqrt(L)/sd.M
#####################################################################################
#####################################################################################
		ww.V = w.est[tmp.mm]-w.est[-tmp.mm]-(W[tmp.mm]-W[-tmp.mm])
		stat.V.one = max(-ww.V/sd.M)*sqrt(L)
		vv.one.K = 1*(stat.V.one>cut.v.one)
		K.test.M[hhhh,hh] = vv.one.K
#####################################################################################
#####################################################################################
		R.one = 1 + sum(1*(-w.V.N>cut.v.one))
		K.test.rank.M[hhhh,hh] = 1*(R.one>K)
	}
}

		           			










#W.est.M = matrix(0,n,NN)                          ## MLE estimator 
write.csv(W.est.M,"Westm.csv")
#R.one.uniform.M = matrix(0,n,NN)                  ## normalized one sided CI for all items 
write.csv(R.one.uniform.M,"Roneuniform.csv")
#R.one.uniform.V = numeric(NN)                     ## normalized one sided CI length for all items   
write.csv(R.one.uniform.V,"Roneuniformv.csv")
#RR.one.uniform.M = matrix(0,n,NN)                 ## not normalized 
write.csv(RR.one.uniform.M,"RRoneuniform.csv")
#RR.one.uniform.V = numeric(NN)
write.csv(RR.one.uniform.V,"RRoneuniform.csv")
#R.two.V = numeric(NN)                             ## empirical size of two sided CI (theta)
write.csv(R.two.V,"Rtwov.csv")
#R.two.length.V = numeric(NN)                      ## two sided CI length 
write.csv(R.two.length.V,"Rtwolengthv.csv")
#R.CI.V = numeric(NN)                              ## empirical size of two sided CI (rank)
write.csv(R.CI.V,"RCIV.csv")
#R.one.V = numeric(NN)                             ## empirical size of one sided CI (theta)   
write.csv(R.one.V,"Ronev.csv")
##################################################################################### ## compare with Chao Gao
#R.CI.C.M = matrix(0,C.len,NN)                     ## empirical size of two sided CI (rank) based on GA 
write.csv(R.CI.C.M,"Rcicm.csv")
#R.length.C.M = matrix(0,C.len,NN)                 ## length CI based on GA 
write.csv(R.length.C.M,"Rlengthcm.csv")
#R.C.I.M = matrix(0,C.len,NN)                      ## empirical size CI (theta)
write.csv(R.C.I.M,"Rcim.csv")
#R.CI.N.M = matrix(0,C.len,NN)                     ## empirical size of two sided CI (rank) Chao Gao  
write.csv(R.CI.N.M,"Rcinm.csv")
#R.length.N.M = matrix(0,C.len,NN)                 ## length CI Chao Gao 
write.csv(R.length.N.M,"Rlengthnm.csv")
#R.N.I.M = matrix(0,C.len,NN)                      ## empirical size CI (theta) Chao gao 
write.csv(R.N.I.M,"rnim.csv")

#K.test.M = matrix(0,len.m,NN)                     ## empirical size of one-sided CI (theta)
write.csv(K.test.M,"Ktestm.csv")
#K.test.rank.M = matrix(0,len.m,NN)                ## empirical size of one-sided CI (rank)
write.csv(K.test.rank.M,"Ktestrankm.csv")

#KK.len = length(KK.seq)                                                          
#KK.set.R = matrix(0,KK.len,NN)   ##--------## set cardinality uniform normalized CI 
write.csv(KK.set.R,"KKsetr.csv")
#KK.power.R = matrix(0,KK.len,NN)                                                   ##--------## uniform normalized CI empirical power 
write.csv(KK.power.R,"KKpowerR.csv")
#KK.set.RR = matrix(0,KK.len,NN)                                                    ##--------## set cardinality uniform CI 
write.csv(KK.set.RR,"KKsetRR.csv")
#KK.power.RR = matrix(0,KK.len,NN)                                                  ##--------## uniform CI empirical power
write.csv(KK.power.RR,"KKpowerRR.csv")



















































































































