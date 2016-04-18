function out = DIIS(D,F,S, X, ncycle)
if ncycle > 4
    B=zeros(6,6);
    nb=size(D,1);
    for i = 1:5
        e{i}=F{ncycle-5+i}*D{ncycle -5+i}*S - S*D{ncycle-5+i}*F{ncycle-5+i};
        e{i}=X'*e{i}*X;
    end
    for i =1:5
        for  j=1:5
            B(i,j)=trace(e{i}*e{j}');
        end
    end
    B(6,:)=-1 * ones(1,6);
    B(:,6)=-1 * ones(6,1);
    B(6,6)=0;
    d=zeros(6,1);
    d(6,1)=-1;
    d=linsolve(B,d);
    Fout=zeros(nb,nb);
    for i = 1:5
        Fout=Fout + d(i)*F{ncycle-5+i};
    end
    out = Fout;
else 
    B=zeros(ncycle + 1, ncycle +1);
    nb=size(D,1);
    for i = 1:ncycle
        e{i}=F{i}*D{i}*S - S*D{i}*F{i};
        e{i}=X'*e{i}*X;
    end  
    for i =1:ncycle
       for  j=1:ncycle
           B(i,j)=trace(e{i}*e{j}');
       end
    end
    B(ncycle+1,:)=-1 * ones(1,ncycle+1);
    B(:,ncycle+1)=-1 * ones(ncycle+1,1);
    B(ncycle+1,ncycle+1)=0;
    d=zeros(ncycle+1,1);
    d(ncycle+1,1)=-1;
    d=linsolve(B,d);
    Fout=zeros(nb,nb);
    for i = 1:ncycle
        Fout=Fout + d(i)*F{ncycle};
    end
    if ncycle ==2
        Fout = 0.75 * Fout + 0.25 *F{1};
        disp('Damping initial DIIS cycle by 25%');
    end
    out = Fout;
end

    
    