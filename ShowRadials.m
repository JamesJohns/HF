function ShowRadials(vecs, energies, basis);
r=0:.001:15;
nb=size(vecs,1);


for i= 1:nb
f=zeros(size(r));
for n=1:nb
f=f+vecs(n,i)*exp(-basis{n}.g.alpha * r.*r)*basis{n}.g.N;
end
figure(1);plot(r,f.*r.*r.*f);
title(strcat('Energy = ', num2str(energies(i))));
pause
end