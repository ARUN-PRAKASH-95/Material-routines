function [sig3,A33,sdvl]=subroutine_stress(eps3,sdvl,ttype)




% material parameters
matp       = inputmat_planestress();
young_x    = matp(1);
young_y    = matp(2);
pr_xy      = matp(3);
g_xy       = matp(4);
sig_11_f_t = matp(5);
sig_11_f_c = matp(6);
sig_22_f_t = matp(7);
sig_22_f_c = matp(8);
pr_yx      = (young_y * pr_xy) / young_x;




% restore the strain tensor in voigt notation
eps = [eps3(1); eps3(2); eps3(3);];

% restore the internal variables at tn
damage  = sdvl(1);
d1      = sdvl(2);
d2      = sdvl(3);
F_f     = sdvl(4);
F_m     = sdvl(5);



C = zeros(3,3);
% Elastic stiffness matix (6*6)
C(1,1) = young_x/(1 - pr_xy*pr_yx);
C(1,2) = pr_xy * young_y /(1 - pr_xy*pr_yx);
C(1,3) = 0;
C(2,1) = pr_xy * young_y /(1 - pr_xy*pr_yx);
C(2,2) = young_y/(1 - pr_xy*pr_yx);
C(2,3) = 0;
C(3,1) = 0;
C(3,2) = 0;
C(3,3) = g_xy;


%%%%%  Failure strains  %%%%%%%%
eps_11_f_t = sig_11_f_t / C(1,1);
eps_11_f_c = sig_11_f_c /  C(1,1);
eps_22_f_t = sig_22_f_t /  C(2,2);
eps_22_f_c = sig_22_f_c /  C(2,2);



eta = 0;

%%%%%%%%%%%   Compute stress   %%%%%%%%

sig3 = zeros(3,1);
for i = 1:3
  for j = 1:3
    sig3(i) = sig3(i) +  ((1+eta-damage)*C(i,j)*eps(j));  % At the beginning of loading damage will be zero, while unloading damage value is the value recorded during the end of tensile loading
  end
end


M = zeros(3,3);
M(1,1) = 1/(1-d1);
M(2,2) = 1/(1-d2);
M(3,3) = 1/sqrt((1-d1)*(1-d2));


% Create an empty effective stress vector
sig3_eff = zeros(3,1);
for i = 1:3
   for j = 1:3
      sig3_eff(i) = sig3_eff(i) + M(i,j)*sig3(j); 
   end
end

sig3_eff;



%  Damage initiation criteria %

if sig3_eff(1) >= 0

  F_f_new   =  sig3_eff(1)/sig_11_f_t;

elseif sig3_eff(1) < 0
   
  F_f_new   =  sig3_eff(1)/sig_11_f_c;
  
endif




if sig3_eff(2) >= 0 
  
  F_m_new   =  sig3_eff(2)/sig_22_f_t; 
 
elseif sig3_eff(2) < 0
  
  F_m_new   =  sig3_eff(2)/sig_22_f_c;

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



%%%%%%%%  Check whether damage has initiated or not  %%%%%%%%%

if F_f<1 && F_m<1 
  
  sig3 = sig3;
  C_T =  (1 + eta - damage)*C;
  

else
    fprintf('yes\n');

%%%%%%  Terms in damage evolution equations  %%%%%%%%
    if sig3_eff(1) >= 0

       epsilon_f = eps_11_f_t
       P1 = 25;
    elseif sig3_eff(1) < 0
   
       epsilon_f = eps_11_f_c; 
       P1 = -25;
    endif
    
    
    
    if sig3_eff(2) >= 0 
      
      epsilon_m = eps_22_f_t;
      P2 = 25;
    elseif sig3_eff(2) < 0
  
      epsilon_m = eps_22_f_c;
      P2 = -25;
    endif


    
    
  
    %%%%%%%% Damage evolution equations  %%%%%%%%%%
    
    if F_f > 1
      
      d1_new =  1  -  (exp(-P1*(eps(1) - epsilon_f)));     %d1
      
     if d1_new >= d1
          d1 = d1_new;
      else
          d1 = d1;
      end
    else
      d1 = 0;      
    endif
    
    
    if F_m > 1
       
      d2_new = 1  -  (exp(-P2*(eps(2) - epsilon_m)));    %d2
      
      if d2_new >= d2
          d2 = d2_new;
      else
          d2 = d2;
      end
    else
      d2 = 0;
    endif

    damage_new = 1  - ((1 - d1)*(1 - d2))

  
% To make sure the damage evolution is greater than or equal to zero
    if damage_new >= damage
      damage = damage_new;
    else
      damage = damage;
    end
  
    sig3 = zeros(3,1);
    for i = 1:3
      for j = 1:3
        sig3(i) = sig3(i) +  ((1 + eta - damage)*C(i,j)*eps(j));
      end
    end
  
  
%%%%%%%%%   Tangent stiffness   %%%%%%%%%
    a = (1 + eta - damage)*C;

%%%%%% Second term of tangent stiffness  %%%%%    
    b = [ P1*exp(-P1*(eps(1)-epsilon_f))*(1 - d2); P2*exp(-P2*(eps(2)-epsilon_m))*(1 - d1); 0];
    
    c = C*eps; 
    
    d = c*b';
    
%%%%  Tangent stiffness   %%%%
    C_T = a - d;
    
    

end



%sig6
%fprintf('s11 %f\n', sig6(1));
%fprintf('C11 %f\n', C(1,1));
%fprintf('C12 %f\n', C(1,2));
%fprintf('C13 %f\n', C(1,3));
%fprintf('eps11 %f\n', eps_tensor(1,1));
%fprintf('eps22 %f\n', eps_tensor(2,2));
%fprintf('eps33 %f\n', eps_tensor(3,3));
%fprintf('damage %f\n', damage);
%fprintf('ef %f\n', epsilon_f);
%fprintf('************************\n');



A33=zeros(3,3);




if ttype==0
    for i=1:3
        for j=1:3
        A33(i,j) = C_T(i,j);
        end
    end
elseif ttype == 1
    hper=1e-8;
    %perturbation of the strain entries (here total strain, maybe this has to be modified)
    for ieps=1:1:length(eps3)
        epsper=eps3;
        epsper(ieps)=epsper(ieps)+hper;
        %recursiv call of your material routine with ttype=0 to avoid
        %endless loop
        %Calculate perturbed stress, sdv are not overwritten
        [sig3per,Adummy,sdvldummy]=subroutine_stress(epsper,sdvl,0);
        %Simple differential quotient
        A33_num(:,ieps)=(sig3per-sig3)/hper;
        
    end
    A33=A33_num;
end


%end

% store history variables
sdvl(1) = damage;
sdvl(2) = d1;
sdvl(3) = d2;
sdvl(4) = F_f;
sdvl(5) = F_m;

end
