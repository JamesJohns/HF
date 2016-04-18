function output = ChemGammainc(n, x)
%ChemGammainc(n,x) Returns the nth order boys function with argument x
if x==0
    output=1/(2*n+1);
else
    output=gamma(n+.5).*gammainc(x,n+.5)./(2*x.^(n+.5));
end
