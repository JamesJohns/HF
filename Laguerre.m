 function [pt ,wt]= lagptwt2(n,alf)
 format long e
 xn=[0:1:n-1];
 a=2*xn+alf+1; rtb=sqrt(xn.*(xn+alf)); rtb(1)=[];
 J=diag(rtb ,-1)+diag (a)+ diag(rtb ,1);
 [f,lambda]=eig(J); pt= diag(lambda);
size(lambda);
 wt=gamma(alf+1)*f(1,:).^2;
 size(wt);
 ptwt=[pt ,wt'];