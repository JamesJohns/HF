function [Eout, Dout, Fout,Gout]=SCF(Din,Fin,H0,X,S, gabcd,N,maxcycle, Nuc);
ncycle=1;
d=1;
nb=size(Fin,1);
D{1}=Din;%Density Corresponding to the Core Hamiltonian
F{1}=Fin;% Fin = Core Hamiltonian
for i = 2:maxcycle
    D{i}=zeros(nb,nb);
    F{i}=zeros(nb,nb);
end;

%e=zeros(nb,nb);
B=0;
while ncycle < maxcycle & d~=0
    %When Starting use solution to core H to build F1;  D1 builds F1; D2
    %builds F2 etc;
     G=Build_Coulomb_Exchange(D{ncycle},gabcd);
     F{ncycle}=H0+G;
    
     
     
%         fprintf(strcat('The Energy difference is ',num2str(E(ncycle)-E(ncycle-1)),'\n'));
%     
%         if abs(E(ncycle) - E(ncycle-1)) < 1d-8
%             d=0;
%         end
%     end
    if ncycle > 1
        F{ncycle}=DIIS(D,F,S,X,ncycle);
    end
     E(ncycle)=Fock_Energy(D{ncycle}, H0, F{ncycle});
    Fprime{ncycle}=X'*F{ncycle}*X;
    [Cprime, epsilon]=eig(Fprime{ncycle}); %[V D]= eig(A);
     C=X'*Cprime;
    [C epsilon]=sortEigs(diag(epsilon),C);
     if ncycle > 2
        if abs(E(ncycle) - E(ncycle-1)) < 1d-10
            d=0;
        end
    
    fprintf(strcat('The Energy for cycle #', num2str(ncycle),' is ',num2str(E(ncycle)),'\n'));
    
    end
    ncycle=ncycle+1;
    
    D{ncycle}=Build_Density(N,C);
%     if ncycle <5
%     D{ncycle}= 0.8*D{ncycle-1} + 0.2 * D{ncycle};
%     end
    
   
end
Eout=E(1: ncycle-1);
Dout=D{ncycle-1};
Fout=F{ncycle-1};
Gout=G;