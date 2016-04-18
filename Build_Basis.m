function [out N]=Build_Basis(Z,AL)
%Out=Build_Basis(Z,AL) Building the set of basis functions to be used for
%the HF Calculation.  For each atom with atomic number Z add primitive
%gaussians corresponding to Pople's 6-31G basis set taken from ESSL
% #  6-31G  EMSL  Basis Set Exchange Library  4/18/16 9:46 AM
% # Elements                             References
% # --------                             ----------
% # H - He: W.J. Hehre, R. Ditchfield and J.A. Pople, J. Chem. Phys. 56,
% # Li - Ne: 2257 (1972).  Note: Li and B come from J.D. Dill and J.A.
% # Pople, J. Chem. Phys. 62, 2921 (1975).
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
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,18.7311370);
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,2.8253937);
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0, 0.6401217 );
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0, 0.1612778 );
            basis{nb}.n=1;
            SP2=[.5];
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
        case 3 %Lithium 
            N=N+3;
            S=[642.4189200              
                  96.79851500                   
                  22.0911210                  
                  6.2010703              
                  1.9351177                
                   0.6367358 ];
             SP2=[ 2.3249184                   
                         0.6324306                     
                         0.0790534                   
                          0.0359620];
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
    
        case 4  %Beryllium
            N=N+4;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,1264.5857000);
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,189.9368100);
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,43.1590890);
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,12.0986630  );
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,3.8063232);
            basis{nb}.n=1;
            nb=nb+1;
            basis{nb}.c=1;
            basis{nb}.g=primitive(Ax, Ay, Az, 0,0,0,1.2728903 );
            basis{nb}.n=1;
            SP2=[  3.1964631        
            0.7478133                   
            0.2199663 
            0.0823099];
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
            for i = 1:4
                nb=nb+1;
                basis{nb}.n=1;
                basis{nb}.g=primitive(Ax,Ay,Az,0,0,0,SP2(i));
                basis{nb}.c=1;
            end
            for i = 1:4
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
        case 5 %Boron
            N=N+5;
            S=[2068.8823000              
                  310.6495700                   
                  70.6830330                   
                  19.8610800              
                  6.2993048                  
                   2.1270270 ];
             SP2=[  4.7279710                    
                   1.1903377                    
                   0.3594117                  
                   0.1267512];
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
        case 9 %Fluorine
            N=N+9;
            S=[7001.7130900               
                1051.3660900           
                  239.2856900                   
                  67.3974453               
                  21.5199573              
                  7.40310130 ];
              SP2=[20.8479528 
                    4.80830834 
                    1.34406986         
                    0.358151393 ];
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