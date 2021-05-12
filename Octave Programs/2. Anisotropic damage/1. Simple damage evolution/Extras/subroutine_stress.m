function [sig6,A66,sdvl]=subroutine_stress(eps6,sdvl,ttype)



%%%% material parameters %%%%

matp       = inputmat();
young_x    = matp(1);
young_y    = matp(2);
young_z    = matp(3);
pr_xy      = matp(4);
pr_yz      = matp(5);
pr_xz      = matp(6);
g_xy       = matp(7);
g_yz       = matp(8);
g_xz       = matp(9);
sig_11_f_t = matp(10);
sig_11_f_c = matp(11);
sig_22_f_t = matp(12);
sig_22_f_c = matp(13);
sig_33_f_t = matp(12);
sig_33_f_c = matp(13);
sig_12_f   =  sig_13_f = sig_23_f = matp(14);
G_c_1      = matp(15);
G_c_2      = matp(17);
G_c_3      = matp(17);
L_c        = matp(19);



pr_yx = (young_y * pr_xy) / young_x;
pr_zy = (young_z * pr_yz) / young_y;
pr_zx = (young_z* pr_xz) / young_x;


xy_yx = pr_xy*pr_yx;
yz_zy = pr_yz*pr_zy;
zx_xz = pr_zx*pr_xz;
xyz   = 2*pr_xy*pr_yz*pr_zx;
E_xyz = young_x*young_y*young_z;

delta = (1 - (xy_yx) - (yz_zy) - (zx_xz) - (xyz)) / E_xyz;



% restore the strain tensor in voigt notation
eps = [eps6(1); eps6(2); eps6(3); eps6(4); eps6(5); eps6(6);];

% restore the strain tensor
eps_tensor = [eps6(1) eps6(4)/2 eps6(6)/2;
              eps6(4)/2 eps6(2) eps6(5)/2;
              eps6(6)/2 eps6(5)/2 eps6(3)];


              
% restore the internal variables at tn
d1  = sdvl(1);
d2  = sdvl(2);
d3  = sdvl(3);
F_f = sdvl(4);
F_m = sdvl(5);
F_z = sdvl(6);




C = zeros(6,6);
% Elastic stiffness matix (6*6)
C(1,1) = ((1 -yz_zy) / (young_y*young_z*delta))*(1 - d1)**2;
C(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d1)*(1 - d2);
C(1,3) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d1)*(1 - d3);
C(1,4) = 0;
C(1,5) = 0;
C(1,6) = 0;
C(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d1)*(1 - d2);
C(2,2) = ((1 -zx_xz) / (young_x*young_z*delta))*(1 - d2)**2;
C(2,3) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d3)*(1 - d2);
C(2,4) = 0;
C(2,5) = 0;
C(2,6) = 0;
C(3,1) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d3)*(1 - d1);
C(3,2) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d3)*(1 - d2);
C(3,3) = ((1 -xy_yx) / (young_x*young_y*delta))*(1 - d3)**2;
C(3,4) = 0;
C(3,5) = 0;
C(3,6) = 0;
C(4,1) = 0;
C(4,2) = 0;
C(4,3) = 0;
C(4,4) = g_xy*(1 - d1)*(1 - d2);
C(4,5) = 0;
C(4,6) = 0;
C(5,1) = 0;
C(5,2) = 0;
C(5,3) = 0;
C(5,4) = 0;
C(5,5) = g_yz*(1 - d1)*(1 - d3);
C(5,6) = 0;
C(6,1) = 0;
C(6,2) = 0;
C(6,3) = 0;
C(6,4) = 0;
C(6,5) = 0;
C(6,6) = g_xz*(1 - d3)*(1 - d2);
C;

eps_11_f_t = sig_11_f_t / ((1 -yz_zy) / (young_y*young_z*delta));
eps_11_f_c = sig_11_f_c / ((1 -yz_zy) / (young_y*young_z*delta));
eps_22_f_t = sig_22_f_t / ((1 -zx_xz) / (young_x*young_z*delta));
eps_22_f_c = sig_22_f_c / ((1 -zx_xz) / (young_x*young_z*delta));
eps_33_f_t = sig_33_f_t / ((1 -xy_yx) / (young_x*young_y*delta));
eps_33_f_c = sig_33_f_c / ((1 -xy_yx) / (young_x*young_y*delta));
eps_12_f   = sig_12_f / g_xy;
eps_13_f   = sig_13_f / g_yz;
eps_23_f   = sig_23_f / g_xz;


%%%%%%%%%%%%  Compute stress   %%%%%%%%%%%%%

% Create an empty stress vector
sig6 = zeros(6,1);
for i = 1:6
   for j = 1:6
      sig6(i) = sig6(i) + C(i,j)*eps(j); 
   end
end
    

M = zeros(6,6);
M(1,1) = 1/(1-d1);
M(2,2) = 1/(1-d2);
M(3,3) = 1/(1-d3);
M(4,4) = 1/sqrt((1-d1)*(1-d2));
M(5,5) = 1/sqrt((1-d2)*(1-d3));
M(6,6) = 1/sqrt((1-d1)*(1-d3));
M;
% Create an empty effective stress vector
sig6_eff = zeros(6,1);
for i = 1:6
   for j = 1:6
      sig6_eff(i) = sig6_eff(i) + M(i,j)*sig6(j); 
   end
end

sig6_eff


%  Damage initiation criteria %

if sig6_eff(1) >= 0

  F_f_new   =  sig6_eff(1)/sig_11_f_t;

elseif sig6_eff(1) < 0
   
  F_f_new   =  sig6_eff(1)/sig_11_f_c;
  
endif




if sig6_eff(2) >= 0 
  
  F_m_new   =  sig6_eff(2)/sig_22_f_t; 
 
elseif sig6_eff(2) < 0
  
  F_m_new   =  sig6_eff(2)/sig_22_f_c;

endif



if sig6_eff(3) >= 0

  F_z_new =  sig6_eff(3)/sig_33_f_t;

elseif sig6_eff(3) < 0

  F_z_new  =  sig6_eff(3)/sig_33_f_c;
  
endif

%%%%%%%% To make sure damage initiation criteria is greater than or equal to previous step  %%%%%%%%%%

if F_f_new >= F_f
    F_f = F_f_new;
else
    F_f = F_f;
end

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

F_f
F_m
F_z

%%%%%%%%  Check whether damage has initiated or not  %%%%%%%%%

if F_f<=1 && F_m<=1 && F_z<=1
  

    sig6 = sig6;
    C_T  =  C;
    

    
else
    %fprintf('yes\n');

%%%%%%  Terms in damage evolution equations  %%%%%%%%
    if sig6_eff(1) >= 0

       epsilon_f = eps_11_f_t;
       P1 = 25;
    elseif sig6_eff(1) < 0
   
       epsilon_f = eps_11_f_c; 
       P1 = -25;
    endif
    
    
    
    if sig6_eff(2) >= 0 
      %fprintf('Tension\n')
      epsilon_m = eps_22_f_t;
      P2 = 25;
    elseif sig6_eff(2) < 0
      %fprintf('Compression\n')
      epsilon_m = eps_22_f_c;
      P2 = -25;
    endif

    
    
    if sig6_eff(3) >= 0

      epsilon_z =  eps_33_f_t;
      P3 = 25;
    elseif sig6_eff(3) < 0

      epsilon_z =  eps_33_f_c; 
      P3 = -25;
    endif   

   

    %%%%%%%% Damage evolution equations  %%%%%%%%%%
    if F_f > 1
      
      d1_new =  1  -  (exp(-P1*(eps(1) - epsilon_f)))    %d1

       
      if d1_new >= d1
          d1 = d1_new;
      else
          d1 = d1;
      end

    endif
    
    if F_m > 1

      d2_new = 1  -  (exp(-P2*(eps(2) - epsilon_m)))   %d2
      
      if d2_new >= d2
          d2 = d2_new;
      else
          d2 = d2;
      end
      

    endif

    
    if F_z > 1
      
      d3_new = 1  -  (exp(-P3*(eps(3) - epsilon_z)))   %d3
      
      if d3_new >= d3
          d3 = d3_new;
      else
          d3 = d3;
      end
      
    endif
  
    
    
    C_d = zeros(6,6);
    C_d(1,1) = ((1 -yz_zy) / (young_y*young_z*delta))*(1 - d1)**2;
    C_d(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d1)*(1 - d2);
    C_d(1,3) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d1)*(1 - d3);
    C_d(1,4) = 0;
    C_d(1,5) = 0;
    C_d(1,6) = 0;
    C_d(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d1)*(1 - d2);
    C_d(2,2) = ((1 -zx_xz) / (young_x*young_z*delta))*(1 - d2)**2;
    C_d(2,3) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d3)*(1 - d2);
    C_d(2,4) = 0;
    C_d(2,5) = 0;
    C_d(2,6) = 0;
    C_d(3,1) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d3)*(1 - d1);
    C_d(3,2) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d3)*(1 - d2);
    C_d(3,3) = ((1 -xy_yx) / (young_x*young_y*delta))*(1 - d3)**2;
    C_d(3,4) = 0;
    C_d(3,5) = 0;
    C_d(3,6) = 0;
    C_d(4,1) = 0;
    C_d(4,2) = 0;
    C_d(4,3) = 0;
    C_d(4,4) = g_xy*(1 - d1)*(1 - d2);
    C_d(4,5) = 0;
    C_d(4,6) = 0;
    C_d(5,1) = 0;
    C_d(5,2) = 0;
    C_d(5,3) = 0;
    C_d(5,4) = 0;
    C_d(5,5) = g_yz*(1 - d1)*(1 - d3);
    C_d(5,6) = 0;
    C_d(6,1) = 0;
    C_d(6,2) = 0;
    C_d(6,3) = 0;
    C_d(6,4) = 0;
    C_d(6,5) = 0;
    C_d(6,6) = g_xz*(1 - d3)*(1 - d2);
    C_d;

    
    
    if d1 == 0
      
      C_T_1 = zeros(6,6);   
      
    else
    
      %%%%%%%%% First term C_T_1 ((d_C_d/d1 : eps) outerProduct (d_d1/d_epsilon))   %%%%%%%%%%

      d_C_d_d1  =  zeros(6,6);

      d_C_d_d1(1,1) = -2*((1 -yz_zy) / (young_y*young_z*delta))*(1 - d1);
      d_C_d_d1(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(d2 - 1);
      d_C_d_d1(1,3) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(d3 - 1);
      d_C_d_d1(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(d2 - 1);
      d_C_d_d1(3,1) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(d3 - 1);
      d_C_d_d1(4,4) = g_xy*(d2 -1);
      d_C_d_d1(5,5) = g_yz*(d3 -1);

    %%%%  (d_C_d/d1 : eps)  %%%%%
      C_T_1_a = zeros(6,1);
      for i = 1:6
         for j = 1:6
            C_T_1_a(i) = C_T_1_a(i) + d_C_d_d1(i,j)*eps(j); 
         end
      end

      
      %%%%%%%%%%%%%%%%   Derivative of d1 with respect to strain (d_d1/d_epsilon)  %%%%%%%%%%%%%
      
      %%%%%%   For Tension   %%%%%%
      if sig6_eff(1) > 0

        C_T_1_b  = [ P1*exp(-P1*(eps(1) - eps_11_f_t)); 0; 0; 0; 0; 0; ];
        
        
      %%%%%   For Compression  %%%%%%
      elseif sig6_eff(1) < 0

        C_T_1_b  =  [  P1*exp(-P1*(eps(1) - eps_11_f_c)); 0; 0; 0; 0; 0;];
       
      endif
        
      C_T_1  =  C_T_1_a*C_T_1_b';
    
    endif
    
    

    if d2 == 0
      
      C_T_2  =  zeros(6,6);
     
    else
      

      %%%%%%%%% Second term C_T_2 ((d_C_d/d2 : eps) outerProduct (d_d2/d_epsilon))   %%%%%%%%%%

      d_C_d_d2  =  zeros(6,6);

      d_C_d_d2(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(d1 - 1);
      d_C_d_d2(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(d1 - 1); 
      d_C_d_d2(2,2) = -2*((1 -zx_xz) / (young_x*young_z*delta))*(1 - d2);
      d_C_d_d2(2,3) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(d3 - 1);   
      d_C_d_d2(3,2) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(d3 - 1);
      d_C_d_d2(4,4) = g_xy*(d1 - 1);
      d_C_d_d2(6,6) = g_xz*(d3 - 1);

   %%%%%  (d_C_d/d2 : eps)  %%%%%
      C_T_2_a = zeros(6,1);
      for i = 1:6
         for j = 1:6
            C_T_2_a(i) = C_T_2_a(i) + d_C_d_d2(i,j)*eps(j); 
         end
      end

    %%%%%%%%%%%%%%%%   Derivative of d2 with respect to strain (d_d2/d_epsilon)  %%%%%%%%%%%%%

      %%%%%%   For Tension   %%%%%%   
      if sig6_eff(2)  >= 0
        
        
        C_T_2_b  = [0; P2*exp(-P2*(eps(2) - eps_22_f_t)); 0; 0; 0; 0;];
        

      %%%%%   For Compression  %%%%%%  
      elseif sig6_eff(2) < 0
        
  
        C_T_2_b  =  [0; P2*exp(-P2*(eps(2) - eps_22_f_c)); 0; 0; 0; 0;];
          
      endif

      C_T_2  =  C_T_2_a*C_T_2_b';
      
    endif
    
    
    
   if d3  == 0
     
     C_T_3  =  zeros(6,6);
   
   else
         
      
      %%%%%%%%% Third term C_T_3 ((d_C_d/d3 : eps) outerProduct (d_d3/d_epsilon))  %%%%%%%%%%
      
      d_C_d_d3  =  zeros(6,6);
      
      d_C_d_d3(1,3)  =  ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(d1 - 1);
      d_C_d_d3(2,3)  =  ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(d2 - 1);
      d_C_d_d3(3,1)  =  ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(d1 - 1);
      d_C_d_d3(3,2)  =  ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(d2 - 1);
      d_C_d_d3(3,3)  =  -2*((1 -xy_yx) / (young_x*young_y*delta))*(1 - d3);
      d_C_d_d3(5,5)  =  g_yz*(d1 - 1);
      d_C_d_d3(6,6)  =  g_xz*(d2 - 1);
      
      %%%%% (d_C_d/d3 : eps) %%%%
      C_T_3_a = zeros(6,1);
      for i = 1:6
         for j = 1:6
            C_T_3_a(i) = C_T_3_a(i) + d_C_d_d3(i,j)*eps(j); 
         end
      end

    %%%%%%%%%%%%%%%%   Derivative of d3 with respect to strain (d_d3/d_epsilon)  %%%%%%%%%%%%%
     
      %%%%%%   For Tension   %%%%%%       
      if sig6_eff(3) >= 0 
        
        C_T_3_b = [0; 0; P3*exp(-P3*(eps(3) - eps_33_f_t)); 0; 0; 0;];
        
      %%%%%   For Compression  %%%%%%
      elseif sig6_eff(3) < 0
        
        C_T_3_b = [0; 0; P3*exp(-P3*(eps(3) - eps_33_f_c)); 0; 0; 0;];

      endif
      
      
      C_T_3  =  C_T_3_a*C_T_3_b';
     
    endif
    
    
    

      
    sig6 = zeros(6,1);
    % Compute stress using Hookes law  %
    for i = 1:6
       for j = 1:6
          sig6(i) = sig6(i) + C_d(i,j)*eps(j);
       end
    end
    
    
    
    
    %%%%%%%%%  Tangent stiffness %%%%%%%%%

    C_T  =  C_d + C_T_1 + C_T_2 + C_T_3;
    

  
endif
  
  
  
  

if ttype==0
    for i=1:6
        for j=1:6
        A66(i,j) = C_T(i,j);
        end
    end
    
    
%%%%%%%%%%%   Numerical tangent   %%%%%%%%%%%%%%

elseif ttype == 1
    hper=1e-8;
    %perturbation of the strain entries (here total strain, maybe this has to be modified)
    for ieps=1:1:length(eps6)
        epsper=eps6;
        epsper(ieps)=epsper(ieps)+hper;
        %recursiv call of your material routine with ttype=0 to avoid
        %endless loop
        %Calculate perturbed stress, sdv are not overwritten
        [sig6per,Adummy,sdvldummy] = subroutine_stress(epsper,sdvl,0);
        %Simple differential quotient
        A66_num(:,ieps)=(sig6per-sig6)/hper;
        
    end
    A66=A66_num;

end




% store history variables
sdvl(1) = d1;
sdvl(2) = d2;
sdvl(3) = d3;
sdvl(4) = F_f;
sdvl(5) = F_m;
sdvl(6) = F_z;


end
