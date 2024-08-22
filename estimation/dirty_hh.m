function hh = dirty_hh(hh_n)


eigen=eig(hh_n);
                check_eig = sum(eigen<0);
                if check_eig>0
                  [V,D] = eig(hh_n);
                  DD=abs(D);
                  hh=V*DD*inv(V);
                  
                  disp('Keep your fingers crossed!')
                else
                    
                    hh=hh_n;
                    disp('No adjustment necessary!')
                    
                end
                    
                
 sqrt(diag(inv(hh)))               