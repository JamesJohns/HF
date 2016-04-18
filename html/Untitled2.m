for n=1:10
    clear x;
    syms x;
    a=basis{n}.g(1).alpha;
    f=exp(-a*x*x);
    Ng=basis{n}.g(1).N
    tt(n)=Ng^2*3*pi/a/2*double(int(f*f*(a-2*a^2*x^2),-Inf, Inf));
end
 
for n =1:nb
    b=basis{n}.g(1).alpha
    fxyz=x^(basis{n}.g(1).lx)*y^(basis{n}.g(1).ly)*z^(basis{n}.g(1).lz)*exp(-b*(x^2 +y^2 +z^2));
    for m=1:nb
    clear x y z
    syms x y z
    a=basis{m}.g(1).alpha;
    rxyz=x^(basis{m}.g(1).lx)*y^(basis{m}.g(1).ly)*z^(basis{m}.g(1).lz)*exp(-a*(x^2 +y^2 +z^2));
    gradxyz=diff(rxyz,x,2)+diff(rxyz,y,2)+diff(rxyz,z,2);
    tt(n,m)=-1/2*double(int(int(int(fxyz * gradxyz,z,-Inf,Inf),y,-Inf,Inf),x,-Inf,Inf))*(basis{n}.g(1).N)*basis{m}.g(1).N;
    end;end