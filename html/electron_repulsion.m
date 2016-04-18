function output=electron_repulsion(basis);
%output=electron_repulsion(basis,rho)Computing repulsion integrals via 
%mcmurchie davidson expansion.
%Step 1:  gaussian product rule to get two gaussians centered at p and p'
%expanded in Hermite gaussians
%Step 2:  Evaluate according to gabcd=2pi^2.5/(pp'sqrt(p+p') *
%(sum(t,u,v)E^ab_tuv*sum(t',u',v')E^cd_t'u'v' Rt+t'u+u'v+v'(alpha,Rpp')
%2Electron integrals tested against gaussian09 for S and P orbitals

nb=size(basis,2);

ee=zeros(nb,nb);
gabcd=zeros(nb,nb,nb,nb);
%Begin nested for loop hell to make gabcd, where each basis can be a
%contracted gaussian.
for a=1:nb  %Loop over basis a
for b=1:nb  %Loop over basis b
for c=1:nb
for d=1:nb
%     if gabcd(a,b,c,d)==0
     % Check if it's already been evaluated.    
     for nba=1:basis{a}.n; %For each gaussian in basis a
    for nbb=1 : basis{b}.n
         %%Use gaussian product theorem to define new gaussian on the line
         %%AB
            p=basis{a}.g(nba).alpha + basis{b}.g(nbb).alpha;
            Px=(basis{a}.g(nba).alpha * basis{a}.g(nba).x0 + basis{b}.g(nbb).alpha*basis{b}.g(nbb).x0)/p;
            Py=(basis{a}.g(nba).alpha * basis{a}.g(nba).y0 + basis{b}.g(nbb).alpha*basis{b}.g(nbb).y0)/p;
            Pz=(basis{a}.g(nba).alpha * basis{a}.g(nba).z0 + basis{b}.g(nbb).alpha*basis{b}.g(nbb).z0)/p;            
         %%Expand ga*gb in hermite gaussians centered at P with exponent p;
            Eabx=gprod_1D(basis{a}.g(nba).x0, basis{a}.g(nba).alpha,basis{a}.g(nba).lx,basis{b}.g(nbb).x0, basis{b}.g(nbb).alpha,basis{b}.g(nbb).lx);
            Eaby=gprod_1D(basis{a}.g(nba).y0, basis{a}.g(nba).alpha,basis{a}.g(nba).ly,basis{b}.g(nbb).y0, basis{b}.g(nbb).alpha,basis{b}.g(nbb).ly);
            Eabz=gprod_1D(basis{a}.g(nba).z0, basis{a}.g(nba).alpha,basis{a}.g(nba).lz,basis{b}.g(nbb).z0, basis{b}.g(nbb).alpha,basis{b}.g(nbb).lz);
            lxmax=size(Eabx,2)-1;
            lymax=size(Eaby,2)-1;
            lzmax=size(Eabz,2)-1;
      for nbc=1:basis{c}.n
      for nbd=1:basis{d}.n
            %Use gaussian product theorem to define new gaussian on the
            %line CD;
            pp=basis{c}.g(nbc).alpha + basis{d}.g(nbd).alpha;
            PPx=(basis{c}.g(nbc).alpha * basis{c}.g(nbc).x0 + basis{d}.g(nbd).alpha*basis{d}.g(nbd).x0)/pp;
            PPy=(basis{c}.g(nbc).alpha * basis{c}.g(nbc).y0 + basis{d}.g(nbd).alpha*basis{d}.g(nbd).y0)/pp;
            PPz=(basis{c}.g(nbc).alpha * basis{c}.g(nbc).z0 + basis{d}.g(nbd).alpha*basis{d}.g(nbd).z0)/pp;  
            %Expand gc*gd in hermite gaussians centered at PP with exponent
            %pp;
            Ecdx=gprod_1D(basis{c}.g(nbc).x0, basis{c}.g(nbc).alpha,basis{c}.g(nbc).lx,basis{d}.g(nbd).x0, basis{d}.g(nbd).alpha,basis{d}.g(nbd).lx);
            Ecdy=gprod_1D(basis{c}.g(nbc).y0, basis{c}.g(nbc).alpha,basis{c}.g(nbc).ly,basis{d}.g(nbd).y0, basis{d}.g(nbd).alpha,basis{d}.g(nbd).ly);
            Ecdz=gprod_1D(basis{c}.g(nbc).z0, basis{c}.g(nbc).alpha,basis{c}.g(nbc).lz,basis{d}.g(nbd).z0, basis{d}.g(nbd).alpha,basis{d}.g(nbd).lz);
            llxmax=size(Ecdx,2)-1;
            llymax=size(Ecdy,2)-1;
            llzmax=size(Ecdz,2)-1;
            alpha = p*pp/(p+pp);
            for lx=0:lxmax
                for ly=0:lymax
                    for lz=0:lzmax
                        for llx=0:llxmax
                            for lly=0:llymax
                                for llz=0:llzmax
                                    gabcd(a,b,c,d)=gabcd(a,b,c,d)+(2*pi^2.5/(p*pp*sqrt(p+pp)) *Eabx(lx+1)*Eaby(ly+1)*Eabz(lz+1)*(-1)^(llx + lly + llz) *Ecdx(llx+1)*Ecdy(lly+1)*Ecdz(llz+1)*Rntuv(0,lx+llx, ly+lly,lz+llz,alpha, [Px Py Pz], [PPx, PPy, PPz]))*basis{a}.g(nba).N*basis{b}.g(nbb).N*basis{c}.g(nbc).N*basis{d}.g(nbd).N;
                                end
                            end
                        end
                    end
                end
            end
      end
      end
    end
%      end
     %Fill in 8 fold symmetric permutations
%       if gabcd(a,b,c,d) ==0
%           gabcd(a,b,c,d) =1d-15;
%       end
%       gabcd(b,a,c,d)=gabcd(a,b,c,d);
%       gabcd(a,b,d,c)=gabcd(a,b,c,d);
%       gabcd(b,a,d,c)=gabcd(a,b,c,d);  
%       gabcd(c,d,a,b)=gabcd(a,b,c,d);
%       gabcd(d,c,a,b)=gabcd(a,b,c,d);
%       gabcd(c,d,b,a)=gabcd(a,b,c,d);
%       gabcd(d,c,b,a)=gabcd(a,b,c,d);
     end
end
end
end
a % Print a to monitor progress Precent done = a/nb *100;
end
output=gabcd;
            