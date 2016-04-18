%function E=driver();
clear all;
natom=1; %Total number of atoms
%AtomList: Coordinate list of atomic positions
%for NStep=1:20

AL=[0 0 0
        0 0 2];
%         0 -1.906366216 -0.674002238
%        -1.650961567 0.953183118 -0.674002238  
%         1.650961567  0.953183118  -0.674002238];
%   AL=[1.64244187    0.43604651    0.00000000
%         -0.68244187    0.43604651    0.00000000
%         -1.96289646    1.34098234    0.00000000]  *.529;
%List of nuclei with charge Z
    Z=[9,1]%,4];%,1, 1 ];  
%Build the basis functions 
%Starting with uncontracted  6-31G exponents

[basis, N]=Build_Basis(Z, AL);
tic
nb=size(basis,2);

%Compute Overlap Integrals
disp('Building Overlap Integrals');
S=Build_Overlap(basis);

Nuc=Energy_Nuclear(Z,AL);
%Compute Kinetic Energy Integrals
disp('Building Kinetic Energy Integrals');
T=Kinetic(basis);

%Compute Nuclear Attraction
disp('Building Nuclear Attraction Integrals');
vab=nuclear_attraction(basis,AL,Z);
disp('Building Core Hamiltonian');
%Make Independent electron, core Hamiltonian
H0=vab + T;

%Compute 2 electron Integrals (Time consuming part)
gabcd=electron_repulsion_new(basis);

%Start Evaluation fo the energy
maxcycle=3000;
ncycle=1;
%Find Orthogonal basis functions by diagonalizing the S matrix
[V D]=eig(S);
X=V*D^(-.5)*V'; %X is used to orthogonalize the Fock matrix
F=H0;  %Starting guess is the solution to core hamiltonian.  Not great
Fprime=X'*F*X;%Transform F into orthogonal basis functions
[Cprime, epsilon]=eig(Fprime);  %Get Eigenvalues of Fprime
C=X*Cprime; %Convert from orthogonal to 6-31G Basis
[C epsilon]=sortEigs((diag(epsilon)),C); %Sort the Eigenvalues 
D=Build_Density(N,C); %Build initial Density
E(ncycle)=Fock_Energy(D,H0,F); %Get the Energy of the Core Hamiltonian
ECore=E(ncycle)+Nuc;

diff=1; %Dummy variable to break the SCF loop if Energy is converged

%Run the SCF Loop.  Should probably be its own function in order to extend
%to geometry optimization.
disp('Entering SCF Loop');
% [E D F G]=SCF(D,F,H0,X,S,gabcd,N,maxcycle);

while ncycle < maxcycle & diff ~=0
    ncycle=ncycle+1;
%     if ncycle > 2;
%     D_Old_Old=D_Old;
%     end
%     D_Old=D;
%     
    D=Build_Density(N,C);
%     if ncycle <3
%             D=0.5*D;% + 0.5*D_Old;
%     else
%             D= 0.5*D +0.25*D_Old + 0.25 * D_Old_Old;
%     end
    
    G=Build_Coulomb_Exchange(D,gabcd);
    %DelG=GOld-G;
    %F=F+0.3*DelG;
    F=H0+G;
    E(ncycle)=Fock_Energy(D, H0, F);
    fprintf(strcat('The Energy for cycle #', num2str(ncycle),' is ',num2str(E(ncycle)),'\n'));
    fprintf(strcat('The Energy difference is ',num2str(E(ncycle)-E(ncycle-1)),'\n'));
    if abs(E(ncycle) - E(ncycle-1)) < 1d-8
        diff=0;
    end
    Fprime=X'*F*X;
    [Cprime, epsilon]=eig(Fprime); %[V D]= eig(A);
    Y=X'*S*X;
    C=X'*Cprime;
    [C epsilon]=sortEigs(diag(epsilon),C);
end
disp(strcat('The total energy is ',num2str(E(end)+Nuc)));
 toc
% % EH2(NStep)=E(end)+Nuc;
% % end;