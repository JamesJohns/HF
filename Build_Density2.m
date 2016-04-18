function out=Build_Density2(N, eigenvecs)
%Build a density matrix from a list of sorted eigenvectors with N electrons
nb=size(eigenvecs,2);
D=zeros(nb,nb);
for n = 1:nb
    for m = 1:nb
        for i=2:N/2+1
            D(n,m)=D(n,m) + eigenvecs(n,i)*eigenvecs(m,i);
        end
    end
end
out=D;
