a=fopen('wateroverlap.txt','r');
N=10;
Nuc= 8.002367061810450;
while ~feof(a)
t=str2num(fgetl(a));
S(t(1),t(2))=t(3);
end
fclose(a);
nb=size(S,1);
for a= 1:nb
    for b=1:nb
        if a~=b
            if S(a,b)~=0
                S(b,a)=S(a,b);
            end
        end
    end
end
a=fopen('waterkinetic.txt','r');
while ~feof(a)
t=str2num(fgetl(a));
T(t(1),t(2))=t(3);
end
fclose(a);

for a= 1:nb
    for b=1:nb
        if a~=b
            if T(a,b)~=0
                T(b,a)=T(a,b);
            end
        end
    end
end
a=fopen('waternuclear.txt','r');
while ~feof(a)
t=str2num(fgetl(a));
vab(t(1),t(2))=t(3);
end
fclose(a);

for a= 1:nb
    for b=1:nb
        if a~=b
            if vab(a,b)~=0
                vab(b,a)=vab(a,b);
            end
        end
    end
end
H0=T+vab;
a=fopen('water2electron.txt','r');
while ~feof(a)
t=str2num(fgetl(a));
gabcd(t(1),t(2),t(3),t(4))=t(5);
end
fclose(a);

for a= 1:nb
    for b=1:nb
        for c=1:nb
            for d=1:nb
               if gabcd(a,b,c,d)~=0
                   gabcd(b,a,c,d)=gabcd(a,b,c,d);
                    gabcd(a,b,d,c)=gabcd(a,b,c,d);
                    gabcd(b,a,d,c)=gabcd(a,b,c,d);  
                    gabcd(c,d,a,b)=gabcd(a,b,c,d);
                    gabcd(d,c,a,b)=gabcd(a,b,c,d);
                    gabcd(c,d,b,a)=gabcd(a,b,c,d);
                    gabcd(d,c,b,a)=gabcd(a,b,c,d);
               end
            end
        end
    end
end
[V D]=eig(S);
X=V*D^(-.5)*V'; %X is used to orthogonalize the Fock matrix
ncycle=1;
maxcycle=100;
F=H0;  %Starting guess is the solution to core hamiltonian.  Not great
Fprime=X'*F*X;%Transform F into orthogonal basis functions
[Cprime, epsilon]=eig(Fprime);  %Get Eigenvalues of Fprime
C=X*Cprime; %Convert from orthogonal to 6-31G Basis
[C epsilon]=sortEigs(real(diag(epsilon)),C); %Sort the Eigenvalues 
D=Build_Density(N,C); %Build initial Density
E(ncycle)=Fock_Energy(D,H0,F); %Get the Energy of the Core Hamiltonian
ECore=Fock_Energy(D,H0,F);

diff=1; %Dummy variable to break the SCF loop if Energy is converged

Run the SCF Loop.  Should probably be its own function in order to extend
to geometry optimization.


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
%[E D F G]=SCF(D,F,H0,X,S,gabcd,N,maxcycle);
disp(strcat('The total energy is ',num2str(E(end)+Nuc)));