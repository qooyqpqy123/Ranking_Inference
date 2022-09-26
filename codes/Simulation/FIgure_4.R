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
MLE = function(A,Y1,Y2,Y3){
	w0 = rnorm(n); w0 = w0 - mean(w0)
	Gw = gradient(w0,A,Y1,Y2,Y3)
	V = max(abs(Gw))
	w = w0
	k0 = 1
	while(V > 0.000001 && k0 <= 1000){
		w = w + 0.01*Gw
		##w = w - mean(w)
		Gw = gradient(w,A,Y1,Y2,Y3)
		V = max(abs(Gw))
		##print(V)
		k0 = k0 + 1
	}
	return(w)
}
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
GA = function(m,A,Y1,Y2,Y3,YY1,YY2,YY3,L,B,W,para){
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
	##F.M.N = t(F.M)
	B.M = matrix(rnorm(L*B),L,B)
	B.M.P = F.M.N%*%B.M
	cut.tmp = apply(abs(B.M.P),2,max)
	cut.v = quantile(cut.tmp,para)/sqrt(L)
	w.V = w.est[m]-w.est[-m]-(W[m]-W[-m])
	stat.V = max(abs(w.V/sd.M))*sqrt(L)
	##stat.V = max(abs(w.V))
	vv = 1*(stat.V>cut.v)
	return(vv)
}
#####################################################################################
paraseq=seq(0.05,0.5,0.05)
n = 60; p = 0.05; L = 80; B = 500; NN = 500
m = 1 
vvvv = matrix(NN,length(paraseq))
for(hh in 1:NN){
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
YY1 = array(0,c(n,n,n,L))
YY2 = array(0,c(n,n,n,L))
YY3 = array(0,c(n,n,n,L))
W = runif(n,2,4); W = W - mean(W)  ## true theta
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
				##Y.tmp.mean = rowMeans(Y.tmp)
				YY1[o,oo,ooo,oooo] = Y.tmp[1]
				YY2[o,oo,ooo,oooo] = Y.tmp[2]
				YY3[o,oo,ooo,oooo] = Y.tmp[3]
				Y1[o,oo,ooo] = Y1[o,oo,ooo] + YY1[o,oo,ooo,oooo]/L
				Y2[o,oo,ooo] = Y2[o,oo,ooo] + YY2[o,oo,ooo,oooo]/L
				Y3[o,oo,ooo] = Y3[o,oo,ooo] + YY3[o,oo,ooo,oooo]/L
			}
		}
	}
}
for (i in length(paraseq)){
vvvv[hh,i] = GA(m,A,Y1,Y2,Y3,YY1,YY2,YY3,L,B,W,paraseq[i])
print(hh)}
write.csv(vvvv,'~/Desktop/gauss_app.csv')
}
#print(mean(vvvv))


