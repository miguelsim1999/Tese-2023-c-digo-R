#g<-x[[2]]
#g2<-y[[2]]
#w1<-xw[[2]]
#<-yw[[2]]
ss<-seq(0,1000,0.1)
s<-seq(0,1300,0.01)

SL1<-function(x){
  SL1<-c()
  for(i in 1:length(x)){
    SL1<-pmax(g-x,0)}
  SL1<-mean(SL1)
  return(SL1)
}

SL2<-function(x){
  SL2<-c()
  for(i in 1:length(x)){
    SL2<-pmax(g2-x,0)}
  SL2<-mean(SL2)
  return(SL2)
}
sl3<-function(x){
  Sl3<-c()
  Sl3<-SL2(x)-SL1(x)
  return(Sl3)
}
SLtransform<-function(x1){
  z<-x1
  U<-sapply(z,SL1)
  return(U)
}
SLtransformy<-function(x1){
  z<-x1
  U<-sapply(z,SL2)
  return(U)
}
SLdif<-function(x1){
  z<-x1
  U<-sapply(z,sl3)
  return(U)
}
SLLX<-function(x){
  SLx<-c()
  for(i in 1:length(x)){
    SLx<-pmax(upa-x,0)}
  SLx<-mean(SLx)
  return(SLx)
}
SLLY<-function(x){
  SLy<-c()
  for(i in 1:length(x)){
    SLy<-pmax(upa2-x,0)}
  SLy<-mean(SLy)
  return(SLy)
}
sldefagreg<-function(x){
  Sl4<-c()
  Sl4<-SLLY(x)-SLLX(x)
  return(Sl4)
}
SLdifagreg<-function(x1){
  z<-x1
  U<-sapply(z,sldefagreg)
  return(U)
}
SLw1<-function(x){
  SLw1<-c()
  for(i in 1:length(x)){
    SLw1<-pmax(w1-x,0)}
  SLw1<-mean(SLw1)
  return(SLw1)
}
SLw2<-function(x){
  SLw2<-c()
  for(i in 1:length(x)){
    SLw2<-pmax(w2-x,0)}
  SLw2<-mean(SLw2)
  return(SLw2)
}
slw3<-function(x){
  Slw3<-c()
  Slw3<-SLw2(x)-SLw1(x)
  return(Slw3)
}
SLdifw<-function(x1){
  z<-x1
  U<-sapply(z,slw3)
  return(U)
}
SLtransfw1<-function(x1){
  z<-x1
  U<-sapply(z,SLw1)
  return(U)
}
SLtransfw2<-function(x1){
  z<-x1
  U<-sapply(z,SLw2)
  return(U)
}
SLagregx<-function(x){
  SL1<-c()
  for(i in 1:length(x)){
    SL1<-pmax(weibx-x,0)}
  SL1<-mean(SL1)
  return(SL1)
}
SLagregy<-function(x){
  SL1<-c()
  for(i in 1:length(x)){
    SL1<-pmax(weiby-x,0)}
  SL1<-mean(SL1)
  return(SL1)
}
SLtransfagregweibx<-function(x1){
  z<-x1
  U<-sapply(z,SLagregx)
  return(U)
}
SLtransfagregweiby<-function(x1){
  z<-x1
  U<-sapply(z,SLagregy)
  return(U)
}
slweibullagreg<-function(x){
  Slw3<-c()
  Slw3<-SLtransfagregweiby(x)-SLtransfagregweibx(x)
  return(Slw3)
}
SLtransfagregweibdiff<-function(x1){
  z<-x1
  U<-sapply(z,slweibullagreg)
  return(U)
}
countador<-function(x){
  y<-x
  for(i in 1:length(y)){
    if(y[i]<400){y[i]=0}
    else{y[i]=y[i]}}
  y = y[y!= 0]
  return(y)
}

genU2<-function(x1,w,security,lamb,shape,scale,ff){
  n<-rpois(x1,lamb)
  s<-c()
  U<-c()
  p<-(1+security)*mean(n)*mean(ff(x1,shape,scale))
  for(i in 1:x1){
    s[i]<-sum(ff(n[i],shape,scale))}
  U[1]<-w+p-s[1]
  for(n in 2:x1) {
    U[n]<-U[n-1]+p-s[n]
  }
  for(j in 1:x1){if(U[j]<0){
    U[j+1]= 0+p-s[j+1]
  }}
  U<-U[!is.na(U)]
  return(U)
}
genli<-function(x1,w,security,lamb,shape,scale,ff){
  n<-rpois(x1,lamb)
  s<-c()
  U<-c()
  p<-(1+security)*mean(n)*mean(ff(x1,shape,scale))
  for(i in 1:x1){
    s[i]<-sum(ff(n[i],shape,scale))}
  U[1]<-w+p-s[1]
  for(n in 2:x1) {
    U[n]<-U[n-1]+p-s[n]
  }
  for(j in 1:x1){if(U[j]<0){
    U[j+1]= 0+p-s[j+1]
  }}
  U<-U[!is.na(U)]
  Li<-c()
  for(m in 1:x1){if(U[m]<0){
    Li[m]= U[m]}else{Li[m]=0}}
  Li = Li[Li!= 0]
  li=c()
  li[1]= 0-Li[1]
  for(r in 2:length(Li)){
    if(length(Li)>1){
      li[r]=Li[r-1]-Li[r]}else{
        li[r]=0}
  }
  li2=c()
  for(i in 1:length(Li)){
    li2[i]=li[i]
    if (sum(li[i])>w){break}
  }
  leng=length(li2)
  li2=sum(li2)
  #return(list(U,li,li2))
  return(list(li2))
}
sap<-function(x1){
  z<-x1
  #U<-replicate(z,ff1(s,w,tet,lamb,shape,scale,ff))
  U<-replicate(z,genli(10^5,400,0.001,1,4,2,rgamma))
  return(U)
}
sap2<-function(x2){
  #U<-replicate(z,ff1(s,w,tet,lamb,shape,scale,ff))
  U<-replicate(x2,genli(10^5,400,0.001,1,2,1,rgamma))
  return(U)
}
sap3<-function(x2){
  #U<-replicate(z,ff1(s,w,tet,lamb,shape,scale,ff))
  U<-replicate(x2,genli(10^5,400,0.001,1,2,2,rweibull))
  return(U)
}
sap4<-function(x2){
  #U<-replicate(z,ff1(s,w,tet,lamb,shape,scale,ff))
  U<-replicate(x2,genli(10^5,400,0.001,1,2/3,4/3,rweibull))
  return(U)
}

SL<-function(x){
  u<-rgamma(10^5,10,5)
  SL<-c()
  for(i in 1:length(x)){
    SL<-pmax(u-x,0)}
  SL<-mean(SL)
  return(SL)
}

lap<-function(x1){
  z<-x1
  U<-lapply(z,SL)
  return(U)
}
SL2<-function(x){
  u<-rgamma(10^5,20,10)
  SL<-c()
  for(i in 1:length(x)){
    SL<-pmax(u-x,0)}
  SL<-mean(SL)
  return(SL)
}

lap2<-function(x1){
  z<-x1
  U<-lapply(z,SL2)
  return(U)
}  

dif<-function(x){
  u<-rgamma(10^5,100,10)
  SL<-c()
  for(i in 1:length(x)){
    SL<-pmax(u-x,0)}
  SL<-mean(SL)
  u2<-rgamma(10^5,1000,100)
  SL2<-c()
  for(i in 1:length(x)){
    SL2<-pmax(u2-x,0)}
  SL2<-mean(SL2)
  TSL<-c()
  TSL<-SL-SL2
  return(TSL)
}
lap3<-function(x1){
  z<-x1
  U<-lapply(z,dif)
  return(U)
}

H<-function(x){
  f1<-function(z)(z-x)*(dgamma(z,2,1)-dgamma(z,4,2))
  res<-integrate(f1,x,Inf)$value
  return(res)
}
x<-seq(0,40,0.01)
#<-lapply(x,H)
#plot(x,sd, type="l",col="red")


H2<-function(x){
  f1<-function(z)(z-x)*(dweibull(z,2/3,4/3)-dweibull(z,2,2))
  res2<-integrate(f1,x,Inf)$value
  return(res2)
}
ull<-lapply(x,H2)
plot(x,ull, type="l",col="red")

