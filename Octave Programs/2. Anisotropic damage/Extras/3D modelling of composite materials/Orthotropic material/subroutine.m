function [sig6,A66,sdvl]=subroutine(eps6,sdvl,ttype)



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
d_f  = sdvl(1);
d_m  = sdvl(2);
d_z  = sdvl(3);
g_f = sdvl(4);
g_m = sdvl(5);
g_z = sdvl(6);
R1  = sdvl(7);
R2  = sdvl(8);
R3  = sdvl(9);
R4  = sdvl(10);
R5  = sdvl(11);
R6  = sdvl(12);
gamma_f_0 = 1;
gamma_m_0 = 1;
gamma_z_0 = 1;


C = zeros(6,6);
% Elastic stiffness matix (6*6)
C(1,1) = ((1 -yz_zy) / (young_y*young_z*delta))*(1 - d_f)**2;
C(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d_f)*(1 - d_m);
C(1,3) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d_f)*(1 - d_z);
C(1,4) = 0;
C(1,5) = 0;
C(1,6) = 0;
C(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d_f)*(1 - d_m);
C(2,2) = ((1 -zx_xz) / (young_x*young_z*delta))*(1 - d_m)**2;
C(2,3) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d_z)*(1 - d_m);
C(2,4) = 0;
C(2,5) = 0;
C(2,6) = 0;
C(3,1) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d_z)*(1 - d_f);
C(3,2) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d_z)*(1 - d_m);
C(3,3) = ((1 -xy_yx) / (young_x*young_y*delta))*(1 - d_z)**2;
C(3,4) = 0;
C(3,5) = 0;
C(3,6) = 0;
C(4,1) = 0;
C(4,2) = 0;
C(4,3) = 0;
C(4,4) = g_xy*(1 - d_f)*(1 - d_m);
C(4,5) = 0;
C(4,6) = 0;
C(5,1) = 0;
C(5,2) = 0;
C(5,3) = 0;
C(5,4) = 0;
C(5,5) = g_yz*(1 - d_f)*(1 - d_z);
C(5,6) = 0;
C(6,1) = 0;
C(6,2) = 0;
C(6,3) = 0;
C(6,4) = 0;
C(6,5) = 0;
C(6,6) = g_xz*(1 - d_z)*(1 - d_m);
C;
vec = [1,2,3,4,5,6];

PD_C = vec*C*vec';

eps_11_f_t = sig_11_f_t / ((1 -yz_zy) / (young_y*young_z*delta))
eps_11_f_c = sig_11_f_c / ((1 -yz_zy) / (young_y*young_z*delta));
eps_22_f_t = sig_22_f_t / ((1 -zx_xz) / (young_x*young_z*delta));
eps_22_f_c = sig_22_f_c / ((1 -zx_xz) / (young_x*young_z*delta));
eps_33_f_t = sig_33_f_t / ((1 -xy_yx) / (young_x*young_y*delta));
eps_33_f_c = sig_33_f_c / ((1 -xy_yx) / (young_x*young_y*delta));
eps_12_f   = sig_12_f / g_xy;
eps_13_f   = sig_13_f / g_yz;
eps_23_f   = sig_23_f / g_xz;



############      Damage initiation criteria      ##############

###   Fibre direction
if eps(1) > 0
  
  g_f_cap = eps(1)/eps_11_f_t;
  
  g_f_new =  g_f_cap  - (gamma_f_0 + R1);
 
  if g_f_new >= g_f
    g_f = g_f_new;
  else
    g_f = g_f;
  endif

else 

  g_f_cap = -eps(1)/eps_11_f_c;
  
  g_f_new =  g_f_cap  - (gamma_f_0 + R4);

  g_f = g_f_new;
 

endif




###   Matrix 1 direction
if eps(2) > 0
  
  g_m_cap = eps(2)/eps_22_f_t;
  
  g_m_new =  g_m_cap  - (gamma_m_0 + R2);
  
  if g_m_new >= g_m
    g_m = g_m_new;
  else
    g_m = g_m;
  endif

else 

  g_m_cap = -eps(2)/eps_22_f_c;
  
  g_m_new =  g_m_cap  - (gamma_m_0 + R5);

  g_m = g_m_new;
endif




###   Matrix 2 direction
if eps(3) > 0
  
  g_z_cap = eps(3)/eps_33_f_t;
  
  g_z_new =  g_z_cap  - (gamma_z_0 + R3);
  
  if g_z_new > g_z
    g_z = g_z_new;
  else
    g_z = g_z;
  endif

else 

  g_z_cap = -eps(3)/eps_33_f_c;
  
  g_z_new =  g_z_cap  - (gamma_z_0 + R6);
  
  g_z = g_z_new;
endif



g_f;
g_m;
g_z;


% Create an empty stress vector
sig6 = zeros(6,1);

if g_f<=0 && g_m<=0 && g_z<=0
  
    
    % Compute stress using Hookes law  %
    for i = 1:6
       for j = 1:6
          sig6(i) = sig6(i) + C(i,j)*eps(j);
       end
    end
    
    C_T  = C;
    
    
else
    
    if eps(1) > 0
      R1 = g_f_cap  -  1;
      gamma_f = gamma_f_0 + R1;
    else
      g_f_cap
      R4 = g_f_cap  -  1

      R1 = R1 + R4;
      gamma_f = gamma_f_0 + R4;
    endif
    
    
    if eps(2) > 0
      R2 = g_m_cap  -  1
    
      gamma_m = gamma_m_0 + R2;
    else
      R5 = g_m_cap  -  1;
      R2 = R2 + R5;
      gamma_m = gamma_m_0 + R5;
    endif

    
    if eps(3) > 0
      R3 = g_z_cap  -  1;
      gamma_z = gamma_z_0 + R3;
    else
      R6 = g_z_cap  -  1;
      R3 = R3 + R6;
      gamma_z = gamma_z_0 + R6;
    endif
    g_f_cap;
    g_m_cap;
    g_z_cap;
    gamma_f;
    gamma_m;
    gamma_z;
    
    alpha_f = 0; 
    alpha_m = 0; 
    alpha_z = 0;
    ####    Damage  evolution   #####
    if g_f > 0
      
      
      d_f_new =  (1 + alpha_f/2)*(1 - 1/gamma_f);
      
      
      if d_f_new >= d_f
          d_f = d_f_new;
      else
          d_f = d_f;
      end  
      
    endif
    
    
    if g_m > 0
      
      
      d_m_new =  (1 + alpha_m/2)*(1 - 1/gamma_m);
      
      if d_m_new > d_m
          d_m = d_m_new;
      else
          d_m = d_m;
      end  
      
    endif
    
    
     
    if g_z > 0
      
      
      d_z_new =  (1 + alpha_z/2)*(1 - 1/gamma_z);
      
      if d_z_new > d_m
          d_z = d_z_new;
      else
          d_z = d_z;
      end  
      
    endif
    d_f;
    d_m;
    d_z;
    
    
    C_d = zeros(6,6);
    % Elastic stiffness matix (6*6)
    C_d(1,1) = ((1 -yz_zy) / (young_y*young_z*delta))*(1 - d_f)**2;
    C_d(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d_f)*(1 - d_m);
    C_d(1,3) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d_f)*(1 - d_z);
    C_d(1,4) = 0;
    C_d(1,5) = 0;
    C_d(1,6) = 0;
    C_d(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta))*(1 - d_f)*(1 - d_m);
    C_d(2,2) = ((1 -zx_xz) / (young_x*young_z*delta))*(1 - d_m)**2;
    C_d(2,3) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d_z)*(1 - d_m);
    C_d(2,4) = 0;
    C_d(2,5) = 0;
    C_d(2,6) = 0;
    C_d(3,1) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta))*(1 - d_z)*(1 - d_f);
    C_d(3,2) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta))*(1 - d_z)*(1 - d_m);
    C_d(3,3) = ((1 -xy_yx) / (young_x*young_y*delta))*(1 - d_z)**2;
    C_d(3,4) = 0;
    C_d(3,5) = 0;
    C_d(3,6) = 0;
    C_d(4,1) = 0;
    C_d(4,2) = 0;
    C_d(4,3) = 0;
    C_d(4,4) = g_xy*(1 - d_f)*(1 - d_m);
    C_d(4,5) = 0;
    C_d(4,6) = 0;
    C_d(5,1) = 0;
    C_d(5,2) = 0;
    C_d(5,3) = 0;
    C_d(5,4) = 0;
    C_d(5,5) = g_yz*(1 - d_f)*(1 - d_z);
    C_d(5,6) = 0;
    C_d(6,1) = 0;
    C_d(6,2) = 0;
    C_d(6,3) = 0;
    C_d(6,4) = 0;
    C_d(6,5) = 0;
    C_d(6,6) = g_xz*(1 - d_z)*(1 - d_m);
 
 
    % Compute stress using Hookes law  %
    for i = 1:6
       for j = 1:6
          sig6(i) = sig6(i) + C_d(i,j)*eps(j);
       end
    end
  
    C_T = C_d; 
    
    
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
        [sig6per,Adummy,sdvldummy] = subroutine(epsper,sdvl,0);
        %Simple differential quotient
        A66_num(:,ieps)=(sig6per-sig6)/hper;
        
    end
    A66=A66_num;
    %A66(1,4)=A66(1,5)=A66(1,6)=0;
    %A66(2,4)=A66(2,5)=A66(2,6)=0;
    %A66(3,4)=A66(3,5)=A66(3,6)=0;
    A66
    vec = [1,2,3,4,5,6];
    
    PD_NT = vec*A66*vec'
end




% store history variables
sdvl(1) = d_f;
sdvl(2) = d_m;
sdvl(3) = d_z;
sdvl(4) = g_f;
sdvl(5) = g_m;
sdvl(6) = g_z;
sdvl(7) =  R1;
sdvl(8) =  R2;
sdvl(9) =  R3;
sdvl(10) = R4;
sdvl(11) = R5;
sdvl(12) = R6;

end
