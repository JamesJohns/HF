function output =primitive( x0, y0, z0, lx,ly, lz, alpha)

output.x0=x0;
output.y0=y0;
output.z0=z0;
output.lx=lx;
output.ly=ly;
output.lz=lz;
output.alpha=alpha;
% syms x y z
% f=((x-x0)^lx * exp(-alpha*(x-x0)^2))^2;
% g=((y-y0)^ly * exp(-alpha*(y-y0)^2))^2;
% h=((z-z0)^lz * exp(-alpha*(z-z0)^2))^2;
% %double(int(f,-inf,inf))
% %double(int(g,-inf, inf))
% %double(int(h,-inf, inf))
% N= double(int(f,-inf, inf))*double(int(g,-inf,inf))*double(int(h,-inf, inf));
% N=1/sqrt(N);
lxf2=1; lyf2=1; lzf2=1;
for i = 2*lx-1:-2:1
    lxf2=lxf2*i;
end
for i = 2*ly-1:-2:1
    lyf2=lyf2*i;
end
for i = 2*lz-1:-2:1
    lzf2=lzf2*i;
end
N=(2/pi)^.75 * (2^(lx+ly+lz)*alpha^((2*lx+2*ly+2*lz+3)/4))/((lxf2*lyf2*lzf2)^.5);
    
output.N = N;

