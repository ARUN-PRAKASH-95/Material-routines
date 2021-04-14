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


% Maximum failure strain in 11 direction
epsilon_f = xsigy0 / C(1,1);

fprintf('epsilon_f %f\n',epsilon_f);

% Create an empty stress vector
sig6 = zeros(6,1);
eta = 0.001;
P = 300;

%Check if the strain in 11 direction is smaller than failure strain

if eps(1) < epsilon_f
    
  % Compute stress using Hookes law
  for i = 1:6
    for j = 1:6
      sig6(i) = sig6(i) +  C(i,j)*eps(j);
    end
  end

elseif eps(1) >= epsilon_f
  
  damage = 1 - (exp(-eps(1)*P))
  
  for i = 1:6
    for j = 1:6
      sig6(i) = sig6(i) +  (1 + eta - damage)*C(i,j)*eps(j);
    end
  end

end
fprintf('C11 %f\n', C(1,1));
fprintf('s11 %f\n', sig6(1));
fprintf('********************\n')



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
