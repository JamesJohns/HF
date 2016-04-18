syms x y z
Z=-4;
for n=1:nb
    a=basis{n}.g.alpha;
    lx=basis{n}.g.lx;
    ly=basis{n}.g.ly;
    lz=basis{n}.g.lz;
    N=basis{n}.g.N
    fn=x^lx*y^ly*z^lz*exp(-a*(x^2+y^2+z^2))*N;
    
    for m=1:nb
        b=basis{m}.g.alpha
        lx=basis{m}.g.lx;
        ly=basis{m}.g.ly;
        lz=basis{m}.g.lz;
        M=basis{m}.g.N;
        fm=x^lx*y^ly*z^lz*exp(-b*(x^2+y^2+z^2))*M;
        gradfm=diff(fm,x, 2) + diff(fm,y,2) + diff(fm,z,2);
        rinv=x^lx*y^ly*z^lz*exp(-b*(x^2+y^2+z^2))*M/(sqrt(x^2 +y^2+z^2))
        SS(n,m)=double(int(int(int(fn*fm,x,-Inf,Inf),y,-Inf,Inf),z,-Inf,Inf));
        TT(n,m)=-1/2*double(int(int(int(fn*gradfm,x,-Inf,Inf),y,-Inf,Inf),z,-Inf,Inf));
        %VV(n,m)=-Z*double(int(int(int(fn*rinv,x,-Inf,Inf),y,-Inf,Inf),z,-Inf,Inf));
    end
end
