function S=goverlap(g1,g2);
%Determine overlap in x, y and z.  Expand each in terms of Hermite
%Gaussians using McMurchie Davidson.  Evaluate Overlap by taking each
%integral = Eij_0*sqrt(pi/p) where p = g1.a + g2.a;
%a=g1.alpha;b=g2.alpha
if    g1.lx < 0
    S=0;
elseif  g1.ly <0
    S=0;
elseif g1.lz < 0
    S=0;
elseif g2.lx <0
    S=0;
elseif g2.ly <0
    S=0;
elseif g2.lz <0
    S=0;
else 
    Ex=gprod_1D(g1.x0,g1.alpha,g1.lx,g2.x0,g2.alpha,g2.lx);
    Ey=gprod_1D(g1.y0,g1.alpha,g1.ly,g2.y0,g2.alpha,g2.ly);
    Ez=gprod_1D(g1.z0,g1.alpha,g1.lz,g2.z0,g2.alpha,g2.lz);

    S=Ex(1)*Ey(1)*Ez(1)*sqrt(pi/(g1.alpha + g2.alpha))^3    *g1.N*g2.N;
end
