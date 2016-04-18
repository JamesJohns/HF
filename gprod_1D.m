function out=gprod_1D(x1,a,ix,x2,b,jx);
%dummy gprod_1D(x1,a,ix,x2,b,jx) Uses the gaussian product rule, and expresses
%the overlap of 2 1D gaussians as a linear combination of hermite gaussians
%centered somewhere on the line x1-x2. Normalization must be done outside
%this routine by the calling function.
%gaussians
%N1 = coefficient of g1, a = exponent of g1, x1 = center of x1
%N2 = coefficient of g2, b = exponent of g2, x2= center of x2
%expand cartesian gaussians in hermite gaussians using McMurchie Davidson
%Algorithm up to a maximum angular momentum of 2 (s-p orbitals only)


tmax=ix+jx;
p= a+b;
q=a*b/p;
A=x1;B=x2;
P=(a*A + b*B)/p;
Q=A-B;
KAB=exp(-q*Q^2);

E=zeros(ix+1, jx+1,tmax+2);
E(1,1,1)=KAB;
i=0;j=0;t=0;
while i < ix
    
    for t=0:i+1+j
        if t==0
            E(i+2,j+1,t+1)=(P-A)*E(i+1,j+1,t+1)+(t+1)*E(i+1,j+1,t+2);
            %E(i+1, j+1, t+1)
            
        else
            E(i+2,j+1,t+1)= 1/(2*p)*E(i+1,j+1,t)+(P-A)*E(i+1,j+1,t+1)+(t+1)*E(i+1,j+1,t+2);
            %E(i+1, j+1, t+1)
            
        end
    end
    i=i+1;
end

while j < jx
          
    for t = 0:i+j+1
        if t== 0;
            E(i+1,j+2,t +1)=(P-B)*E(i+1,j+1,t+1) + (t+1)*E(i+1,j+1,t+2);
        else
            E(i+1,j+2,t+1)=1/(2*p)*E(i+1,j+1,t)+(P-B)*E(i+1,j+1,t+1) + (t+1)*E(i+1,j+1,t+2);
        end
    end
      j=j+1;
end
%out=E
for q=1:i+j+1
    tmp(q)=E(i+1,j+1,q);
end

out = tmp;


% % 
% %varargout.n=tmax;   
% E=zeros(1,tmax+1);
% E000=KAB;
% E100=(P-A)*KAB;
% E010=(P-B)*KAB;
% E101=1/(2*p)*E000;
% E011=1/(2*p)*E000;
% E200=(P-A)*E100 +1*E101;
% E201=.5/p *E100 +(P-A)*E101;
% E202=.5/p *E101;
% E020=(P-B)*E010 + 1*E011;
% E021=.5/p*E010 +(P-B)*E011;
% E022=.5/p*E011;
% E110=(P-A)*E010 + 1*E011;
% E111=.5/p *E010 + (P-A)*E011;
% E112=.5/p*E011;
% 
% if tmax ==0
%     E(1)=KAB;
% end
% if tmax ==1;
%     if ix==1
%         E(1)=E100;E(2)=E101;
%     end
%     if jx==1
%         E(1)=E010;E(2)=E011;
%     end
% end
% if tmax ==2
%     if ix==2
%         E(1)=E200; E(2)=E201; E(3)=E202;
%     end
%     if jx==2
%         E(1)=E020; E(2)=E021; E(3)=E022;
%     end
%     if ix==1
%         E(1)=E110; E(2)=E111; E(3)=E112;
%     end
% end
% out=E;

