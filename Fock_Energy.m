function out=Fock_Energy(D,H0,F);
nb=size(D,1);
E=0;
for n = 1:nb
    for m = 1:nb
        %E=E+ D(n,m)*(F(n,m));
        E=E+ D(n,m)*(H0(n,m) + F(n,m));
    end
end
out=E;