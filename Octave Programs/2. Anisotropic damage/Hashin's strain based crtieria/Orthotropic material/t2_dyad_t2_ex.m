a = eye(3);
b = [1,0,0;0,1,0;0,0,2];
matp      = inputmat();
xE        = matp(1);            % Young's modulus
xnu       = matp(2);            % Poisson's ratio
xsigy0    = matp(3);            % initial yield stress
xk        = xE/(3*(1-2*xnu));   % bulk modulus

c = xk*t2_otimes_t2(a,b)