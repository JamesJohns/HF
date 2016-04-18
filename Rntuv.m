function output=Rntuv(n,t,u,v,p, P, A);
%Rntuv(n, t,u,v,p,P,A)Determine the helper integral Rntuv for the coulomb
%integral of order n, the t,u,v th Hermite polynomial with exponent p
%centered at [Px Py Pz] and charge centered at location [Ax Ay Az];
nmax=n+t+u+v+2;
PA2=sum((P-A).^2);
R=zeros(nmax, t+1, u+1, v+1);
PA2;
for i=0:nmax-1
    R(i+1,1,1,1)=(-2*p)^i * ChemGammainc(i,p*PA2);
end
%R(3)
 t;
for tt = 1:t
    if tt==1
        for i=1:nmax-1
        R(i,2,1,1)=(P(1)-A(1))*R(i+1,1,1,1);
        end
    else
        for i = 1:nmax-1
               R(i,tt+1,1,1)=(tt-1)*R(i+1,tt-1,1,1) +(P(1)-A(1))*R(i+1,tt,1,1);
        end
    end
end
R;
for uu = 1:u
    if uu==1
        for i=1:nmax-1
        R(i,t+1,2,1)=(P(2)-A(2))*R(i+1,t+1,1,1);
        end
    else
        for i = 1:nmax-1
               R(i,t+1,uu+1,1)=(uu-1)*R(i+1,t+1,uu-1,1) +(P(2)-A(2))*R(i+1,t+1,uu,1);
        end
    end
end


for vv = 1:v
    if vv==1
        for i=1:nmax-1
        R(i,t+1,u+1,2)=(P(3)-A(3))*R(i+1,t+1,u+1,1);
        end
    else
        for i = 1:nmax-1
               R(i,t+1,u+1,vv+1)=(vv-1)*R(i+1,t+1,u+1,vv-1) +(P(3)-A(3))*R(i+1,t+1,u+1,vv);
        end
    end
end
%output=R;    
output=R(n+1,t+1,u+1,v+1);
    
    
    
    
    
    
    
    
    
    
    