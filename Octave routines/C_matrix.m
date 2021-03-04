% material parameters
matp      = inputmat();
xE        = matp(1);       % Young's modulus
xnu       = matp(2);       % Poisson's ratio
xsigy0    = matp(3);       % initial yield stress
xk        = xE/(3*(1-2*xnu))   % bulk modulus
xid = eye(3);

% Lames constants
lambda  =  xE * xnu / ((1 + xnu) * (1 - 2*xnu))
xmu  = xE/(2*(1+xnu));     % shear modulus

C = 2*xmu*getP4sym()  + xk*t2_otimes_t2(xid,xid) 

ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
A66=zeros(6,6);
%if ttype==0
    for i=1:6
        for j=1:6
        A66(i,j) = C(ii(i),jj(i),ii(j),jj(j));
        end
    end
A66