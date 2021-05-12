%  Damage initiation criteria %

if sig3_eff(1) >= 0

  F_f_new  =  sqrt((sig3_eff(1)/sig_11_f_t)**2 + (sig3_eff(3)/sig_12_f)**2);

elseif sig3_eff(1) < 0
   
  F_f_new  =   sqrt((sig3_eff(1)/sig_11_f_c)**2);
  
endif


if sig3_eff(2) >= 0 
  
  F_m_new  =   sqrt( (sig3_eff(2)/sig_22_f_t)**2 +  (sig3_eff(3)/sig_12_f)**2);

elseif sig3_eff(2) < 0
  
  F_m_new  =  sqrt( (sig3_eff(2)/(2*sig_12_f))**2 + ((sig3_eff(2)/sig_22_f_c)*( (sig_22_f_c/2*sig_12_f)**2  - 1)) +  (sig3_eff(3)/sig_12_f)**2);

endif