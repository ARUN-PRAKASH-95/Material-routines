function [sig6,A66,sdvl]=vmises(eps6,sdvl,ttype)



% material parameters
matp     = inputmat();
young_x  = matp(1);
young_y  = matp(2);
young_z  = matp(3);
pr_xy    = matp(4);
pr_yz    = matp(5);
pr_xz    = matp(6);
g_xy     = matp(7);
g_yz     = matp(8);
g_xz     = matp(9);

disp(young_x)
% general 
ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
xid = eye(3);

% restore the strain tensor
eps = [eps6(1) eps6(4)/2 eps6(6)/2;
       eps6(4)/2 eps6(2) eps6(5)/2;
       eps6(6)/2 eps6(5)/2 eps6(3)];


pr_yx = (young_y * pr_xy) / young_x;
pr_zy = (young_z * pr_yz) / young_y;
pr_zx = (young_z * pr_xz) / young_x;
	
	
xy_yx = pr_xy*pr_yx;
yz_zy = pr_yz*pr_zy; 
zx_xz = pr_zx*pr_xz;
xyz   = TWO*pr_xy*pr_yz*pr_zx;
E_xyz = young_x*young_y*young_z;
	
delta = (ONE - (xy_yx) - (yz_zy) - (zx_xz) - (xyz)) / E_xyz;
 

dsdeEl(1,1) = (ONE -yz_zy) / (young_y*young_z*delta);
dsdeEl(1,2) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta);
dsdeEl(1,3) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta);
dsdeEl(1,4) = 0;
dsdeEl(1,5) = 0;
dsdeEl(1,6) = 0;
dsdeEl(2,1) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta);
dsdeEl(2,2) = (ONE -zx_xz) / (young_x*young_z*delta);
dsdeEl(2,3) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta);
dsdeEl(2,4) = 0;
dsdeEl(2,5) = 0;
dsdeEl(2,6) = 0;
dsdeEl(3,1) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta);
dsdeEl(3,2) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta);
dsdeEl(3,3) = (ONE -xy_yx) / (young_x*young_y*delta);
dsdeEl(3,4) = 0;
dsdeEl(3,5) = 0;
dsdeEl(3,6) = 0;
dsdeEl(4,1) = 0;
dsdeEl(4,2) = 0;
dsdeEl(4,3) = 0;
dsdeEl(4,4) = TWO * g_yz;
dsdeEl(4,5) = 0;
dsdeEl(4,6) = 0;
dsdeEl(5,1) = 0;
dsdeEl(5,2) = 0;
dsdeEl(5,3) = 0;
dsdeEl(5,4) = 0;
dsdeEl(5,5) = TWO * g_xz;
dsdeEl(5,6) = 0;
dsdeEl(6,1) = 0;
dsdeEl(6,2) = 0;
dsdeEl(6,3) = 0;
dsdeEl(6,4) = 0;
dsdeEl(6,5) = 0;
dsdeEl(6,6) = TWO * g_xy; 

C = zeros(6,6);       
  for i=1:6
      for j=1:6
        C(i,j) = dsdeEl(i,j);
      end
    end
        
    
sig = C * eps;


% restore stress tensor as vector
sig6=zeros(6,1);
for i=1:6
    sig6(i) = sig(ii(i),jj(i));
end

% restore stiffness tensor as matrix
ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
A66=zeros(6,6);
%if ttype==0
    for i=1:6
        for j=1:6
        A66(i,j) = C(ii(i),jj(i),ii(j),jj(j));
        end
    end
%end

% store history variables
sdvl(1:6,1)=[epsp(1,1) epsp(2,2) epsp(3,3)...
             2*epsp(1,2) 2*epsp(2,3) 2*epsp(1,3)]';
sdvl(7:12,1)=[Balpha(1,1) Balpha(2,2) Balpha(3,3)...
             2*Balpha(1,2) 2*Balpha(2,3) 2*Balpha(1,3)]';
sdvl(13) = alpha;
end
