function out=nuclear_attraction(basis, AList, Z);
  nA=size(AList,1);
  
  nb=size(basis,2);
  VAB=zeros(nb,nb);
  for N=1:nA
     for a = 1:nb
         for b = 1:nb
           for nba=1:basis{a}.n
              for nbb=1:basis{b}.n
                  VAB(a,b)=VAB(a,b)+coulombg(basis{a}.g(nba), basis{b}.g(nbb), AList(N,1), AList(N,2), AList(N,3),Z(N));
              end
           end
         end
     end
  end
  out=VAB;
            