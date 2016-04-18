function  [outvectors outeigenvals]=sortEigs(eigvals, eigvecs);
nb=size(eigvals,1);
[sorted_eigvals index]=sort(eigvals);
for n = 1:nb
    sorted_eigvecs(:,n)=eigvecs(:,index(n));
end
outvectors=sorted_eigvecs;
outeigenvals=sorted_eigvals;