if F_m_new >= F_m
    F_m = F_m_new;
else
    F_m = F_m;
end

if F_z_new >= F_z
    F_z = F_z_new;
else
    F_z = F_z;
end



%  Damage initiation criteria %

if sig6_eff(1) >= 0

  F_f_new  =  sqrt((sig6_eff(1)/sig_11_f_t)**2 + (sig6_eff(4)/sig_12_f)**2 + (sig6_eff(5)/sig_13_f)**2);

elseif sig6_eff(1) < 0
   
  F_f_new =  sqrt((sig6_eff(1)/sig_11_f_c)**2);
  
endif




if sig6_eff(2)+sig6_eff(3) >= 0 
  
  F_m_new   =   sqrt(((sig6_eff(2)+sig6_eff(3))**2/(sig_22_f_t*sig_33_f_t)) -  (sig6_eff(2)*sig6_eff(3)/sig_23_f**2) +  (sig6_eff(4)/sig_12_f)**2 + (sig6_eff(5)/sig_13_f)**2 + (sig6_eff(6)/sig_23_f)**2);
 
elseif sig6_eff(2)+sig6_eff(3) < 0
  
  F_m_new =  sqrt(((sig6_eff(2)+sig6_eff(3))**2/(sig_22_f_c*sig_33_f_c)) + ((sig6_eff(2)+sig6_eff(3)/sig_22_f_c)*( (sig_22_f_c/2*sig_12_f)  - 1))   -  (sig6_eff(2)*sig6_eff(3)/sig_23_f**2) +  (sig6_eff(4)/sig_12_f)**2 + (sig6_eff(5)/sig_13_f)**2 + (sig6_eff(6)/sig_23_f)**2);

endif



if sig6_eff(3) >= 0

  F_z_new  =  sqrt((sig6_eff(3)/sig_33_f_t)**2 + (sig6_eff(5)/sig_23_f)**2 + (sig6_eff(6)/sig_13_f)**2);

elseif sig6_eff(3) < 0

  F_z_new  =  sqrt((sig6_eff(3)/sig_33_f_c)**2 + (sig6_eff(5)/sig_23_f)**2 + (sig6_eff(6)/sig_13_f)**2);
  
endif