function out = Energy_Nuclear(Z,AL);
natom=size(AL,1);
energy=0;
for i =1:natom
    for m = i+1:natom
        distance = sqrt(sum((AL(i,:)-AL(m,:)).^2));
        energy = energy+Z(i)*Z(m)/distance;
    end;
end
out=energy;