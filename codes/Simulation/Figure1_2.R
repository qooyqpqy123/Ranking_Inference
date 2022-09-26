gradient = function(w,A,Y1,Y2,Y3){ #Gradient Function
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
MLE = function(A,Y1,Y2,Y3){  ##MLE function
	w0 = rnorm(n); w0 = w0 - mean(w0)
	Gw = gradient(w0,A,Y1,Y2,Y3)
	V = max(abs(Gw))
	w = w0
	while(V > 0.00001){
		w = w + 0.01*Gw
		##w = w - mean(w)
		Gw = gradient(w,A,Y1,Y2,Y3)
		V = max(abs(Gw))
		#print(V)
	}
	return(w)
}
#####################################################################################
n = 60; L = 30; NN = 200
W = runif(n,2,4)                                                ## true theta 
W = W - mean(W) 
r.v = seq(1.4*log(60)/sqrt(n*0.2*L),1.4*log(60)/sqrt(n*0.01*L),length=8)        ## ratio sequence  
EEE = matrix(0,NN,8)    ## p from 0.01 to 0.2 
EEE2=matrix(0,NN,8) #qqplot
EEE3=matrix(0,NN,8)
for(ee in 1:8){
	p = 1/n/n/L/r.v[ee]/r.v[ee]
for(eee in 1:NN){
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
for(o in 1:(n-2)){
	for(oo in (o+1):(n-1)){
		for(ooo in (oo+1):n){
			A[o,oo,ooo] = rbinom(1,1,p)
			p1.tmp = exp(W[o])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p2.tmp = exp(W[oo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p3.tmp = exp(W[ooo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p.tmp = c(p1.tmp,p2.tmp,p3.tmp)
			Y.tmp = rmultinom(L,1,p.tmp)
			Y.tmp.mean = rowMeans(Y.tmp)
			Y1[o,oo,ooo] = Y.tmp.mean[1]
			Y2[o,oo,ooo] = Y.tmp.mean[2]
			Y3[o,oo,ooo] = Y.tmp.mean[3]
	##print(Y1[o,oo,ooo]+Y2[o,oo,ooo]+Y3[o,oo,ooo])
		}
	}
}
w.est = MLE(A,Y1,Y2,Y3)## MLE estimate
EEE2[eee,ee]=w.est-W
write.csv(EEE2,"qqplot_matrix.csv")
EEE[eee,ee] = max(abs(w.est-W)) ## L infty norm 
EEE3[eee,ee]=sum((w.est-W)^2) #L_2 norm
}}

EEE.mean = colMeans(EEE)
EEE3.mean=colMeans(EEE3)
write.csv(EEE.mean,"consistency_p_infty.csv") #save mean-infty
write.csv(EEE3.mean,"consistency_p_2.csv") #save mean-l2
plot(r.v,EEE.mean,type='b',col='red')
plot(r.v, EEE3.mean,type='b',col='blue')

