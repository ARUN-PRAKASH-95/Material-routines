% material parameters
matp      = inputmat();
xE        = 200e6;       % Young's modulus
xnu       = 0.33;       % Poisson's ratio
xk        = xE/(3*(1-2*xnu))   % bulk modulus
xid = eye(3);

% Lames constants
lambda  =  xE * xnu / ((1 + xnu) * (1 - 2*xnu))
mu  = xE/(2*(1+xnu))    % shear modulus
lambda+2*mu

C = 2*mu*getP4sym() +xk*t2_otimes_t2(xid,xid);

ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
A66=zeros(6,6);
%if ttype==0
    for i=1:6
        for j=1:6
        A66(i,j) = C(ii(i),jj(i),ii(j),jj(j));
        end
    end

    
 C2 = zeros(3,3,3,3);
    %if ttype==0
    for i=1:6
        for j=1:6
         C2(ii(i),jj(i),ii(j),jj(j)) = A66(i,j);
        end
    end
    
 C2