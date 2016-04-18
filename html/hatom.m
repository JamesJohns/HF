clear all
% Atomic units
m=1;hbar=1;e=1;
%length = bohr

grid =[-10:.3:10];
dx=grid(2)-grid(1);
k=2*pi/dx;
Gmax=k^2/2
gridsize=size(grid,2);
for x = 1:gridsize
    for y = 1:gridsize;
        for z=1:gridsize
            r=sqrt(grid(x).^2 + grid(y)^2 + grid(z)^2);
            VH(x,y,z)=-1/r;
        end;end;end
