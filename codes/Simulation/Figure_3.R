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
	while(V > 0.00001){
		w = w + 0.03*Gw
		##w = w - mean(w)
		Gw = gradient(w,A,Y1,Y2,Y3)
		V = max(abs(Gw))
		#print(V)
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
#####################################################################################
n = 60; L = 10; NN = 500
W = runif(n,2,4)                                                ## true theta 
W = W - mean(W) 
p_list= c(0.015,0.03)
qq=matrix(0,NN,length(p_list))
#r.v = seq(1/sqrt(n*n*0.2*L),1/sqrt(n*n*0.01*L),length=8)        ## ratio sequence  
#EEE = matrix(0,NN,1)    ## p from 0.01 to 0.2 
#EEE2=matrix(0,NN,8) #qqplot
#EEE3=matrix(0,NN,8)
#	p = 1/n/n/L/r.v[ee]/r.v[ee]
for (index in c(1:length(p_list))){
for(eee in 1:NN){
  print(eee)
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
for(o in 1:(n-2)){
	for(oo in (o+1):(n-1)){
		for(ooo in (oo+1):n){
			A[o,oo,ooo] = rbinom(1,1,p_list[index])
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
var=variance(1,w.est,A)
qq[eee,index]=(w.est[1]-W[1])*sqrt(L*var)
#qq=c(qq,w.est[1]-W[1])
write.csv(qq,"~/Desktop/qq.csv")
#print(eee)
#write.csv(EEE2,"qqplot_matrix.csv")
#EEE[eee,ee] = max(abs(w.est-W)) ## L infty norm 
#EEE3[eee,ee]=sum((w.est-W)^2) #L_2 norm
}}
	
data<-read.csv("~/Desktop/qq_L5.csv")
gau=rnorm(1000000,0,1)
qqplot(gau,data[,2])
qqline(gau,data[,2])

par(mfrow=c(3,3))
qqplot(gau,data[,2],col='blue',ylab="",xlab="",main="p=0.008,L=5")
qqline(gau,data[,2])
qqplot(gau,data[,3],col='blue',ylab="",xlab="",main="p=0.015,L=5")
qqline(gau,data[,3])
qqplot(gau,data[,4],col='blue',ylab="",xlab="",main="p=0.030,L=5")
qqline(gau,data[,4])

data2<-read.csv("~/Desktop/qq_L10.csv")
qqplot(gau,data2[,2],col='blue',ylab="",xlab="",main="p=0.008,L=10")
qqline(gau,data2[,2])
qqplot(gau,data2[,3],col='blue',ylab="",xlab="",main="p=0.015,L=10")
qqline(gau,data2[,3])
qqplot(gau,data2[,4],col='blue',ylab="",xlab="",main="p=0.030,L=10")
qqline(gau,data2[,4])


data3<-read.csv("~/Desktop/qq_L201.csv")
qqplot(gau,data3[,2],col='blue',ylab="",xlab="",main="p=0.008,L=20")
qqline(gau,data3[,2])
qqplot(gau,data3[,3],col='blue',ylab="",xlab="",main="p=0.015,L=20")
qqline(gau,data3[,3])
qqplot(gau,data3[,4],col='blue',ylab="",xlab="",main="p=0.030,L=20")
qqline(gau,data3[,4])



par(mfrow=c(3,3))
hist(data[,2],ylab="",xlab="",main="p=0.008,L=5",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)
#curve(dnorm, -3.5, 3.5, lwd=2, axes = FALSE, xlab = "", ylab = "")
qqplot(gau,data[,2],col='blue',ylab="",xlab="",main="p=0.008,L=5")
qqline(gau,data[,2])
qqplot(gau,data[,3],col='blue',ylab="",xlab="",main="p=0.015,L=5")
hist(data[,3],ylab="",xlab="",main="p=0.015,L=5",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)
qqline(gau,data[,3])
qqplot(gau,data[,4],col='blue',ylab="",xlab="",main="p=0.030,L=5")
qqline(gau,data[,4])
hist(data[,4],ylab="",xlab="",main="p=0.030,L=5",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)
#data2<-read.csv("~/Desktop/qq_L10.csv")
qqplot(gau,data2[,2],col='blue',ylab="",xlab="",main="p=0.008,L=10")
qqline(gau,data2[,2])
hist(data2[,2],ylab="",xlab="",main="p=0.008,L=10",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)
qqplot(gau,data2[,3],col='blue',ylab="",xlab="",main="p=0.015,L=10")
qqline(gau,data2[,3])
hist(data2[,3],ylab="",xlab="",main="p=0.015,L=10",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)
qqplot(gau,data2[,4],col='blue',ylab="",xlab="",main="p=0.030,L=10")
qqline(gau,data2[,4])
hist(data2[,4],ylab="",xlab="",main="p=0.030,L=10",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)


#data3<-read.csv("~/Desktop/qq_L20.csv")
qqplot(gau,data3[,2],col='blue',ylab="",xlab="",main="p=0.008,L=20")
qqline(gau,data3[,2])
hist(data3[,2],ylab="",xlab="",main="p=0.008,L=20",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)

qqplot(gau,data3[,3],col='blue',ylab="",xlab="",main="p=0.015,L=20")
qqline(gau,data3[,3])
hist(data3[,3],ylab="",xlab="",main="p=0.015,L=20",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)

qqplot(gau,data3[,4],col='blue',ylab="",xlab="",main="p=0.030,L=20")
qqline(gau,data3[,4])
hist(data3[,4],ylab="",xlab="",main="p=0.030,L=20",freq=FALSE,breaks=seq(-4,4,0.2))
lines(density(gau), col = 4, lwd = 2)

#EEE.mean = colMeans(EEE)
#EEE3.mean=colMeans(EEE3)
#write.csv(EEE.mean,"consistency_p_infty.csv") #save mean-infty
#write.csv(EEE3.mean,"consistency_p_2.csv") #save mean-l2
#plot(r.v,EEE.mean,type='b',col='red')
#plot(r.v, EEE3.mean,type='b',col='blue')


#gau=rnorm(500,0,1)
#theta1=qq
#qqplot((theta1)/sd(theta1),gau)
#qqline((theta1-W[1])/sd(theta1),gau)

#qq<-read.csv("~/Desktop/qq.csv")
#qq<-qq[,2]
#gau=rnorm(500,0,1)
#qqplot(qq/sd(qq),gau,xlab="",ylab="")
#qqline(qq/sd(qq),gau)
