%Testing Integrand methods
a=.05; b = 100;N=100
dx=.1;

grid=linspace(-20,20,N);
dx=grid(2)-grid(1);
ga=sqrt(2*a/pi)*exp(-2*a*grid.*grid);
gb=sqrt(2*b/pi)*exp(-2*b*grid.*grid);
figure(1);plot(grid,ga,grid,gb);
sga=sum(ga)*dx
sgb=sum(gb)*dx

mida=0;
midb=0;
for i=1:size(grid,2)-1;
    mida=mida + 0.5*(ga(i+1)+ga(i))*dx;
    midb=midb + 0.5*(gb(i+1)+gb(i))*dx;
end
mida
midb

simpa=0;
simpb=0;
simpa=dx/3 *(ga(1) +ga(N));
simpb=dx/3 * (gb(1) + gb(N));
for n = 2:2:N-1
    simpa=simpa+dx/3 * (4*ga(n) + 2*ga(n+1));
    simpb=simpb + dx/3 * (4*gb(n)+2*gb(n+1));
end
simpa
simpb

%map the grid from [-20 : 20] to [-1 1] to
%unitgrid(i)=1/20 * grid(i)
%figure(2);plot(r,grid);
lga=0;
[pt wt]=lgwt(N,-1,1);
for i = 1:n
    lga=lga+wt(i)*exp(-3*(20*pt(i))^2);
end

% lga=lga*20
% clear r;
% for N=4:500
% lgb=0;
% [pt wt]=lgwt(N,-20,20);
% for i = 1:N
%     lgb=lgb+wt(i)*sqrt(3/pi)*exp(-3*(pt(i))^2);
% end
% r(N)=1-lgb;
% end
% figure(4);plot(4:500,log10(r(4:500)))
%laga=0
% for N = 4:500
%  [pt wt]=Laguerre(N,a);
% % laga=0;lagb=0;
%  for i = 1:N
%      laga=laga+ sqrt(2*a/pi)*exp(-a*pt(i)^2)*wt(i)*2;
% 
%  end
%    r(N)=1-laga;
% end
[pt wt]=lgwt(200, -1 ,1);
clear alpha;
clear error;
alpha=logspace(-3, 12);

for m = 1:47
    a=alpha(m);
    r=linspace(-sqrt(20^2/a),sqrt(20^2/a),200);
    lga(m)=0;
    pttor= @(x) sqrt(400/a)*x;
  
    for i=1:200
          
        lga(m)=lga(m) + wt(i)*sqrt(a/pi)*exp(-a*( pttor(pt(i))^2)) * pttor(pt(i))^3;
    end
    lga=lga*sqrt(20^2/a);
    error(m)=1 -lga(m);
end
figure(4);plot(log10(alpha),log10(error));
    
%figure(4);plot(N(4:end),log10(r(4:end)));

%     lagb=lagb +sqrt(2*b/pi)*exp(-b*pt(i)^2)*wt(i)*2;
% end
% 
% [pt wt]=HermQuad(N);
% herma=0;hermb=0;
% for i =1:N
%     herma=herma + sqrt(2*a/pi)*exp(-a*pt(i)^2)*wt(i);
%     hermb=hermb + sqrt(2*b/pi)*exp(-b*pt(i)^2)*wt(i);
% end
% herma 
% hermb