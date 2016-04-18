function out = Build_Fock(H0,D,gabcd,nb);
%out=Build_Fock(H0, D, gabcd, nb); Build a Fock matrix of size nb from a 1 electron
%hamiltonian, a density matrix, and the exchange matrix
for n=1:nb
    for m=1:nb
        F(n,m)=H0(n,m);
        for p=1:nb
            for q=1:nb
                F(n,m)=F(n,m)+D(p,q)*(2*gabcd(n,m,p,q)-gabcd(n,p,u,q));
            end
        end
    end
end
out=F;