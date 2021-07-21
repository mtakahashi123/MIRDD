MIdiagRDD<-function(
	y, x, cut, seed=1, M1=100, M2=5, M3=1, p2s1=1, emp=0, bw="mserd",
	ker="triangular", bwidth=1, p1=1, conf=95, upper=1,covs1=NULL,up=NULL,lo=NULL
	){

set.seed(seed)
df1<-na.omit(data.frame(x,y))

if(p1>2){
 p1<-2
}else{
 p1<-p1
}

x1<-df1$x; y3<-df1$y; cut<-cut

d<-NULL
if(upper==1){
 d[x1<cut]<-0
 d[x1>=cut]<-1
}else if(upper==0){
 d[x1<cut]<-1
 d[x1>=cut]<-0
}

y0<-NULL
y0[d==0]<-y3[d==0]
y1<-NULL
y1[d==1]<-y3[d==1]

if(class(covs1)=="NULL"){
 data3<-data.frame(y3,x1,d)
}else{
 data3<-data.frame(y3,x1,d,covs1)
}

x2<-x1^2
if(class(covs1)=="NULL"&p1==1){
 data2<-data.frame(y0,y1,x1)
}else if(class(covs1)!="NULL"&p1==1){
 data2<-data.frame(y0,y1,x1,covs1)
}else if(class(covs1)=="NULL"&p1!=1){
 data2<-data.frame(y0,y1,x1,x2)
}else if(class(covs1)!="NULL"&p1!=1){
 data2<-data.frame(y0,y1,x1,covs1,x2)
}
data1<-data2

#############
#Loading R package Amelia
if(suppressMessages(suppressWarnings(require(Amelia)))){
    print("Loading required package Amelia")
} else {
    print("trying to install Amelia...")
    install.packages("Amelia")
    if(suppressMessages(suppressWarnings(require(Amelia)))){
        print("Amelia installed and loaded")
    } else {
        stop("could not install Amelia")
    }
}

#############
#Loading R package rdrobust
if(suppressMessages(suppressWarnings(require(rdrobust))) ){
    print("Loading required package rdrobust")
} else {
    print("trying to install rdrobust...")
    install.packages("rdrobust")
    if(suppressMessages(suppressWarnings(require(rdrobust))) ){
        print("rdrobust installed and loaded")
    } else {
        stop("could not install rdrobust")
    }
}

#############
#Bandwidth Selection
rbd1<-rdbwselect(data3$y3,data3$x1,c=cut,bwselect=bw,kernel=ker,p=p1,covs=covs1)
h1<-rbd1$bws[1]*bwidth
h2<-rbd1$bws[2]*bwidth
data3a<-data3[(cut-h1<=(data3$x1))&((data3$x1)<=cut+h2),]
rddn<-nrow(data3a)

rbdMI1<-rbd1
hMI1<-h1
hMI2<-h2
data1a<-data1[(cut-hMI1<=(data1$x1))&((data1$x1)<=cut+hMI2),]
data2a<-data2[(cut-hMI1<=(data1$x1))&((data1$x1)<=cut+hMI2),]
MIn<-nrow(data2a)

#############
#Naive Estimator
y1naive<-na.omit(data2a$y1)
y0naive<-na.omit(data2a$y0)
n1naive<-length(y1naive)
n0naive<-length(y0naive)
Naiven<-n1naive+n0naive
se1naive<-sd(y1naive)/sqrt(n1naive)
se0naive<-sd(y0naive)/sqrt(n0naive)
numnaive<-(se1naive+se0naive)^2
denomnaive<-(se1naive^2)/(n1naive-1)+(se0naive^2)/(n0naive-1)
dfna<-numnaive/denomnaive
naive1<-mean(y1naive)-mean(y0naive)
naive2<-se1naive+se0naive

#############
#MI: Amelia

a.out<-amelia(data2a,p2s=p2s1,m=M1,empri=emp*MIn)
tauhata<-matrix(NA,MIn,M1)
tauhat0<-matrix(NA,MIn,M1)
tauhat1<-matrix(NA,MIn,M1)
nMIa1<-NULL
nMIa2<-NULL
for(i in 1:M1){
 y1imp<-a.out$imputations[[i]]$y1
 y0imp<-a.out$imputations[[i]]$y0
 y1imp[y1imp>up]<-up
 y1imp[y1imp<lo]<-lo
 y0imp[y0imp>up]<-up
 y0imp[y0imp<lo]<-lo
 tauhata[,i]<-y1imp-y0imp
 tauhat0[,i]<-y0imp
 tauhat1[,i]<-y1imp
 nMIa1[i]<-mean(tauhata[,i])
 nMIa2[i]<-var(tauhata[,i])/length(tauhata[,i])
}
mia1<-mean(nMIa1)
wmia<-mean(nMIa2)
bmia<-(1/(M1-1))*sum((nMIa1-mean(nMIa1))^2)
tmia<-wmia+(1+1/M1)*bmia
mia2<-sqrt(tmia)

#############
#Regression Discontinuity Design
modelRDD1<-rdrobust(data3$y3, data3$x1, c=cut,all=TRUE,h=c(h1,h2),kernel=ker,p=p1,covs=covs1)
rddb1<-modelRDD1$Estimate[1]
rddb2<-modelRDD1$Estimate[3]

##########################
#Graphs 1 & 2: Use all of the M1 datasets from MI
layout(matrix(1:4,2,2,byrow=TRUE))
minx3<-min(nMIa1,rddb1,naive1)
maxx3<-max(nMIa1,rddb1,naive1)
hist(nMIa1,xlim=c(minx3,maxx3),main="1.MI, RDD, Naive",xlab="LATE at the Cutoff")
abline(v=rddb1,col=2,lwd=3)
abline(v=naive1,col=1,lwd=3)

minx4<-min(nMIa1,rddb1)
maxx4<-max(nMIa1,rddb1)
hist(nMIa1,xlim=c(minx4,maxx4),main="2.MI and RDD",xlab="LATE at the Cutoff")
abline(v=rddb1,col=2,lwd=3)

#############
#Graphs 3 and 4: Densities
maxy0<-NULL
maxy1<-NULL
for(i in 1:M2){
 maxy0[i]<-max(density(tauhat0[,i])$y)
 maxy1[i]<-max(density(tauhat1[,i])$y)
}
maxy<-max(maxy0,maxy1,density(data3a$y3[(data3a$d)==0])$y,density(data3a$y3[(data3a$d)==1])$y)
minx<-min(tauhat0,tauhat1,data3a$y3[(data3a$d)==0],data3a$y3[(data3a$d)==1])
maxx<-max(tauhat0,tauhat1,data3a$y3[(data3a$d)==0],data3a$y3[(data3a$d)==1])

#Control Group
plot(density(tauhat1),xlim=c(minx,maxx),ylim=c(0,maxy),xlab="",ylab="",main="3.Densities (Control)",type="n")
y0n<-data3a$y[data3a$d==0]
par(new=TRUE)
plot(density(y0n),xlim=c(minx,maxx),ylim=c(0,maxy),lwd=3,xlab="",ylab="",main="",lty=1,col=8)
y1n<-data3a$y[data3a$d==1]
par(new=TRUE)
plot(density(y1n),xlim=c(minx,maxx),ylim=c(0,maxy),lwd=3,xlab="",ylab="",main="",lty=1,col=4)
for(i in 1:M2){
 y0mi<-tauhat0[,i][a.out$missMatrix[,1]==TRUE]
 par(new=TRUE)
 plot(density(y0mi),xlim=c(minx,maxx),ylim=c(0,maxy),xlab="",ylab="",main="",col=2,lty=2,lwd=1)
}

#Treatment Group
plot(density(tauhat1),xlim=c(minx,maxx),ylim=c(0,maxy),xlab="",ylab="",main="4.Densities (Treatment)",type="n")
par(new=TRUE)
plot(density(y0n),xlim=c(minx,maxx),ylim=c(0,maxy),lwd=3,xlab="",ylab="",main="",lty=1,col=8)
par(new=TRUE)
plot(density(y1n),xlim=c(minx,maxx),ylim=c(0,maxy),lwd=3,xlab="",ylab="",main="",lty=1,col=4)

for(i in 1:M2){
 y1mi<-tauhat1[,i][a.out$missMatrix[,2]==TRUE]
 par(new=TRUE)
 plot(density(y1mi),xlim=c(minx,maxx),ylim=c(0,maxy),xlab="",ylab="",main="",col=2,lty=2,lwd=1)
}

#############
#Graphs 5 and 6: Scatterplots
x11mi<-data1a$x1[is.na(data1a$y1)==TRUE]
y11mi<-as.numeric(na.omit(data2a$y1))
x01mi<-data1a$x1[is.na(data1a$y0)==TRUE]
y01mi<-as.numeric(na.omit(data2a$y0))

minx<-min(x11mi,x01mi,data1$x1)
maxx<-max(x11mi,x01mi,data1$x1)
miny<-min(y01mi,y11mi,tauhat1,tauhat0,data3$y3)
maxy<-max(y01mi,y11mi,tauhat1,tauhat0,data3$y3)

windows()
layout(matrix(1:4,2,2,byrow=TRUE))
plot(data3$x1[d==0],data3$y3[d==0],xlim=c(minx,maxx),ylim=c(miny,maxy),col=8,main="5.Observed Values",xlab="x",ylab="y")
par(new=TRUE)
plot(data3$x1[d==1],data3$y3[d==1],xlim=c(minx,maxx),ylim=c(miny,maxy),col=4,pch=2,main="",xlab="",ylab="")
minx2<-cut-h1
maxx2<-cut+h2
abline(v=cut,lwd=2)
abline(v=minx2,lty=2,col=3,lwd=2)
abline(v=maxx2,lty=2,col=3,lwd=2)

plot(data3$x1[d==0],data3$y3[d==0],xlim=c(minx,maxx),ylim=c(miny,maxy),col=8,main="6.Observed & Imputed Values",xlab="x",ylab="y")
par(new=TRUE)
plot(data3$x1[d==1],data3$y3[d==1],xlim=c(minx,maxx),ylim=c(miny,maxy),col=4,pch=2,main="",xlab="",ylab="")
minx2<-cut-h1
maxx2<-cut+h2

par(new=TRUE)
a<-NULL
b<-NULL
for(i in 1:M3){
 y0mi<-tauhat0[,i][a.out$missMatrix[,1]==TRUE]
 par(new=TRUE)
 plot(x01mi,y0mi,xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="",ylab="",main="",col=2,lty=2,lwd=1)
}

par(new=TRUE)
a<-NULL
b<-NULL
for(i in 1:M3){
 y1mi<-tauhat1[,i][a.out$missMatrix[,2]==TRUE]
 par(new=TRUE)
 plot(x11mi,y1mi,xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="",ylab="",main="",col=2,lty=2,lwd=1,pch=2)
}
abline(v=cut,lwd=2)
abline(v=minx2,lty=2,col=3,lwd=2)
abline(v=maxx2,lty=2,col=3,lwd=2)

#############
#Graphs 7 and 8: Scatterplots
#Control Group
#par(new=TRUE)
plot(data3$x1[d==0],data3$y3[d==0],xlim=c(minx,maxx),ylim=c(miny,maxy),col=8,main="7.Observed & Imputed (Control)",xlab="x",ylab="y")
 a<-NULL
 b<-NULL
for(i in 1:M3){
 y0mi<-tauhat0[,i][a.out$missMatrix[,1]==TRUE]
 par(new=TRUE)
 plot(x01mi,y0mi,xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="",ylab="",main="",col=2,lty=2,lwd=1)
}
abline(v=cut,lwd=2)
abline(v=minx2,lty=2,col=3,lwd=2)
abline(v=maxx2,lty=2,col=3,lwd=2)

#Treatment Group
#par(new=TRUE)
plot(data3$x1[d==1],data3$y3[d==1],xlim=c(minx,maxx),ylim=c(miny,maxy),col=4,main="8.Observed & Imputed (Treatment)",xlab="x",ylab="y",pch=2)
a<-NULL
b<-NULL
for(i in 1:M3){
 y1mi<-tauhat1[,i][a.out$missMatrix[,2]==TRUE]
 par(new=TRUE)
 plot(x11mi,y1mi,xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="",ylab="",main="",col=2,lty=2,lwd=1,pch=2)
}
abline(v=cut,lwd=2)
abline(v=minx2,lty=2,col=3,lwd=2)
abline(v=maxx2,lty=2,col=3,lwd=2)

#############
#Graphs 9 & 10: Scatterplots
minx<-cut-h1
maxx<-cut+h2
miny0<-min(y01mi,tauhat0)
maxy0<-max(y01mi,tauhat0)
miny1<-min(y11mi,tauhat1)
maxy1<-max(y11mi,tauhat1)

windows()
layout(matrix(1:4,2,2,byrow=TRUE))
#Control Group
plot(x11mi,y01mi,xlim=c(minx,maxx),ylim=c(miny0,maxy0),xlab="",ylab="",main="9.Around Cutoff (Control)",col=8)
a<-NULL
b<-NULL
for(i in 1:M3){
 y0mi<-tauhat0[,i][a.out$missMatrix[,1]==TRUE]
 par(new=TRUE)
 plot(x01mi,y0mi,xlim=c(minx,maxx),ylim=c(miny0,maxy0),xlab="",ylab="",main="",col=2,lty=2,lwd=1)
 
}
for(i in 1:M2){
 y0mi<-tauhat0[,i][a.out$missMatrix[,1]==TRUE]
 model<-lm(y0mi~x01mi)
 a[i]<-as.numeric(model$coefficients[1])
 b[i]<-as.numeric(model$coefficients[2])
 abline(a=a[i],b=b[i],lwd=1)
}

#Treatment Group
plot(x01mi,y11mi,xlim=c(minx,maxx),ylim=c(miny1,maxy1),xlab="",ylab="",main="10.Around Cutoff (Treatment)",col=4,pch=2)
a<-NULL
b<-NULL
for(i in 1:M3){
 y1mi<-tauhat1[,i][a.out$missMatrix[,2]==TRUE]
 par(new=TRUE)
 plot(x11mi,y1mi,xlim=c(minx,maxx),ylim=c(miny1,maxy1),xlab="",ylab="",main="",col=2,lty=2,lwd=1,pch=2)
}
for(i in 1:M2){
 y1mi<-tauhat1[,i][a.out$missMatrix[,2]==TRUE] 
 model<-lm(y1mi~x11mi)
 a[i]<-as.numeric(model$coefficients[1])
 b[i]<-as.numeric(model$coefficients[2])
 abline(a=a[i],b=b[i],lwd=1)
}

#############
#Graphs 11 & 12: Histograms

#Control Group
b0<-NULL
b1<-NULL
for(i in 1:M1){
 y0mi<-tauhat0[,i][a.out$missMatrix[,1]==TRUE]
 model0<-lm(y0mi~x01mi)
 b0[i]<-as.numeric(model0$coefficients[2])
 y1mi<-tauhat1[,i][a.out$missMatrix[,2]==TRUE]
 model1<-lm(y1mi~x11mi)
 b1[i]<-as.numeric(model1$coefficients[2])
}

hist(b0,main="11.Local Slope (Control)",xlab="Coefficient (Slope)")
hist(b1,main="12.Local Slope (Treatment)",xlab="Coefficient (Slope)")

##########################
#Output
C1<-qt((1-(1-conf/100)/2),dfna)
C2<-qt((1-(1-conf/100)/2),rddn-1)
C3<-qt((1-(1-conf/100)/2),MIn-1)
CInaiveUL<-naive1+C1*naive2
CInaiveLL<-naive1-C1*naive2
CIrddUL<-rddb1+C2*rddb2
CIrddLL<-rddb1-C2*rddb2
CImiUL<-mia1+C3*mia2
CImiLL<-mia1-C3*mia2

rownames<-c("MI","RDD","Naive")
Estimate<-round(c(mia1,rddb1,naive1),4)
Std.Error<-round(c(mia2,rddb2,naive2),4)
CI.UL<-round(c(CImiUL,CIrddUL,CInaiveUL),4)
CI.LL<-round(c(CImiLL,CIrddLL,CInaiveLL),4)
size<-c(MIn,rddn,Naiven)
ratio<-round(size/nrow(data3)*100,1)
bandwidth<-round(c(h1,h1,h1),4)

output<-data.frame(rownames,Estimate,Std.Error,CI.LL,CI.UL,size,ratio,bandwidth)
output

}


