function output=coulomb(g1,g2,Ax,Ay,Az,Z);

%coulomb integral is the integral of -Z*phi1*1/(r-rA) * phi2 dr
%Steps:
%Use the gaussian product theorem to replace phi1*phi2 with a new single
%Hermite gaussian located at P with exponent p=a+b;
%Then evaluate the coulomb integral by using the representation of 1/(rA) =
%1/sqrt(pi) * integral ( exp(-t^2*rA^2)dt)
%The resultant integrals are then d/dPX d/dPY d/dPZ * F0(pRpA^2)

a=g1.alpha;b=g2.alpha;

p=a+b;
P=[(a*g1.x0 + b*g2.x0)/p, (a*g1.y0+b*g2.y0)/p, (a*g1.z0 + b*g2.z0)/p];
A=[Ax, Ay, Az];
RPA2=sum((A-P).^2);
%For lx, ly, lz, lx', ly', lz' = 0
%output=sqrt(4*p/pi)*ChemGammainc(0,pRPA2);

Ex=gprod_1D(g1.x0,g1.alpha,g1.lx,g2.x0,g2.alpha,g2.lx);
Ey=gprod_1D(g1.y0,g1.alpha,g1.ly,g2.y0,g2.alpha,g2.ly);
Ez=gprod_1D(g1.z0,g1.alpha,g1.lz,g2.z0,g2.alpha,g2.lz);

tmax=size(Ex,2)-1;umax=size(Ey,2)-1;vmax=size(Ez,2)-1;
V=0;
for tt= 0:tmax
    for uu=0:umax
        for vv=0:vmax
            V=V+2*pi/p * Ez(vv+1)*Ey(uu+1)*Ex(tt+1)*Rntuv(0,tt,uu,vv,p,P,A);
        end
    end
end
V=V*g1.N * g2.N*(-Z);
output=V;
%V_tuv=(d/dPx)^t(d/dPy)^U(d/dPz)^v*\integral(exp(-p(r-p)^2)/(r-A));

%V_tuv=2*pi/p * (d/dPx)^t(d/dPy)^U(d/dPz)^v* F0(p*RPA2);

%R_tuv=(d/dPx)^t(d/dPy)^U(d/dPz)^v F0
%Develop R_tuv through Recursion
%R100=d/dPx(integral(exp(-p*(Ax-Px)^2t^2)dt) =-2p(Ax-Px)*integral(exp(-p*(Ax-Px)^2t^2)*t^2dt) 
% = -2p(Ax-Px)F1(pRPA2)
%R^n(000)=(-2p)^n * F_n(p*RPA2);
