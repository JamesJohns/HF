function out=Build_Coulomb_Exchange(D, gabcd);
%build the integrals of the non core hamiltonian from the density matrix
%and the 4center integral matrix gabcd.  gabcd must be of size(nb,nb,nb,nb)
%and D be of size(nb,nb);
nb=size(D,1);
G=zeros(nb,nb);
for a=1:nb
    for b = 1:nb
        for c=1:nb
            for d=1:nb
                G(a,b)=G(a,b)+D(c,d)*(2*gabcd(a,b,c,d)-gabcd(a,c,b,d));
            end
        end
    end
end
out=G;
