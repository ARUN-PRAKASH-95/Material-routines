function [sig6,A66,sdvl]=vmises(eps6,sdvl,ttype)


tol=1e-8;

% material parameters
matp      = inputmat();
xE        = matp(1);       % Young's modulus
xnu       = matp(2);       % Poisson's ratio
xsigy0    = matp(3);       % initial yield stress
xk        = xE/(3*(1-2*xnu));   % bulk modulus

% Lames constants
lambda  =  xE * xnu / ((1 + xnu) * (1 - 2*xnu));
xmu  = xE/(2*(1+xnu));     % shear modulus


% general 
ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
xid = eye(3);

% restore the strain tensor
eps = [eps6(1); eps6(2); eps6(3); eps6(4); eps6(5); eps6(6);];

% restore the internal variables at tn
damage  =  sdvl(1);    %Scalar damage variable


C = zeros(6,6);

C(1,1) = lambda + 2*xmu;
C(1,2) = lambda;
C(1,3) = lambda;
C(1,4) = 0;
C(1,5) = 0;
C(1,6) = 0;
C(2,1) = lambda;
C(2,2) = lambda + 2*xmu;
C(2,3) = lambda;
C(2,4) = 0;
C(2,5) = 0;
C(2,6) = 0;
C(3,1) = lambda;
C(3,2) = lambda;
C(3,3) = lambda + 2*xmu;
C(3,4) = 0;
C(3,5) = 0;
C(3,6) = 0;
C(4,1) = 0;
C(4,2) = 0;
C(4,3) = 0;
C(4,4) = xmu;
C(4,5) = 0;
C(4,6) = 0;
C(5,1) = 0;
C(5,2) = 0;
C(5,3) = 0;
C(5,4) = 0;
C(5,5) = xmu;
C(5,6) = 0;
C(6,1) = 0;
C(6,2) = 0;
C(6,3) = 0;
C(6,4) = 0;
C(6,5) = 0;
C(6,6) = xmu;


% Maximum failure strain
epsilon_f = xsigy0 / C(1,1)
eps(1)

% Create an empty stress vector
sig6 = zeros(6,1);


%Check if the strain in 11 directioin is smaller than failure strain

if eps(1) < epsilon_f
    
  % Compute stress using Hookes law
  for i = 1:6
    for j = 1:6
      sig6(i) = sig6(i) +  C(i,j)*eps(j);
    end
  end

else
  
  damage = 1 - (exp(-eps(1)*100))
  
  for i = 1:6
    for j = 1:6
      if i==j
        sig6(i) = sig6(i) +  (1 - damage)*C(i,j)*eps(j);
      else
        sig6(i) = sig6(i) +  C(i,j)*eps(j);
      end
    end
  end

end





% restore stiffness tensor as matrix
ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];

A66=zeros(6,6);
%if ttype==0
    for i=1:6
        for j=1:6
        A66(i,j) = C(i,j);
        end
    end
%end

% store history variables
sdvl(1) = damage;
end
