 function [pt, wt]=HermQuad(n)
 format long e
 xn=[1:1:n]; rtb=sqrt (xn/2); rtb(n)=[];
 J=diag(rtb ,-1)+diag(rtb ,1);
 [f,lambda]=eig(J); pt= diag(lambda);
 wt=sqrt(pi)*f(1,:).^2;
 ptwt=[pt ,wt'];