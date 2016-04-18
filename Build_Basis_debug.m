function [out N]=Build_Basis(Z,AL)
%Out=Build_Basis(Z,AL) Building the set of basis functions to be used for
%the HF Calculation.  For each atom with atomic number Z add primitive
%gaussians corresponding to Pople's 6-31G basis set.  Starting with Be, C,
%H, and O.
natoms=size(Z,2);
nb=0;
N=0;
for i = 1:natoms
    Ax=AL(i,1);Ay=AL(i,2);Az=AL(i,3);
    switch Z(i)
        case 1
            N=N+1;
            %Add 4 s orbitals to basis set centered at [Ax, Ay, Az]
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,5);
            basis{nb}.n=1;
            SP2=2;
             for i = 1:1
%                 nb=nb+1;
%                 basis{nb}.n=1;
%                 basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,SP2(i));
%                 basis{nb}.c=1;
%                 nb=nb+1;
%                 basis{nb}.n=1;
%                 basis{nb}.g=primitive(Ax,Ay,Az,1,0,0,SP2(i));
%                 basis{nb}.c=1;
%                 nb=nb+1;
%                 basis{nb}.n=1;
%                 basis{nb}.g=primitive(Ax,Ay,Az,0,1,0,SP2(i));
%                 basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,1,SP2(i));
                basis{nb}.c=1;
            end
       
        case 4  %Beryllium
            N=N+4;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,5);
            basis{nb}.n=1;
                     SP2=[ 2];
            for i = 1:1
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,1,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,1,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,1,SP2(i));
                basis{nb}.c=1;
            end
        case 6  %Carbon 6-31G
            N=N+6;
            S=[3047.5249000              
                  457.3695100                   
                  103.9486900                   
                  29.2101550              
                   9.2866630                 
                   3.1639270 ];
             SP2=[  7.8682724                    
                          1.8812885                    
                          0.5442493                    
                          0.1687144];
              for i = 1:6
                  nb=nb+1;
                  basis{nb}.n=1;
                  basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,S(i));
                  basis{nb}.c=1;
              end
               for i = 1:4
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,1,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,1,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,1,SP2(i));
                basis{nb}.c=1;
               end
        case 7  %Nitrogen 6-31G
            N=N+7;
            S=[4173.5110000  
                  627.4579000  
                  142.9021000 
                  40.2343300
                  12.8202100 
                   4.3904370];
            SP2=[11.6263580               
                       2.7162800              
                       0.7722180               
                       0.2120313];
             for i = 1:6
                  nb=nb+1;
                  basis{nb}.n=1;
                  basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,S(i));
                  basis{nb}.c=1;
              end
               for i = 1:4
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,1,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,1,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,1,SP2(i));
                basis{nb}.c=1;
               end
        case 8  %Oxygen
            N=N+8;
            S=[5484.6717000               
                  825.2349500            
                  188.0469600                   
                  52.9645000               
                  16.8975700            
                  5.7996353];
              SP2=[15.5396160 
                         3.5999336   
                         1.0137618        
                         0.2700058];
                for i = 1:6
                  nb=nb+1;
                  basis{nb}.n=1;
                  basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,S(i));
                  basis{nb}.c=1;
              end
               for i = 1:4
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,1,0,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,1,0,SP2(i));
                basis{nb}.c=1;
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,1,SP2(i));
                basis{nb}.c=1;
               end
        otherwise
                disp('Atom Not Yet Coded');
    end
end
   out= basis;