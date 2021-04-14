% material parameters
matp      = inputmat();
young_x   = matp(1);
young_y   = matp(2);
young_z   = matp(3);
pr_xy     = matp(4);
pr_yz     = matp(5);
pr_xz     = matp(6);
g_xy      = matp(7);
g_yz      = matp(8);
g_xz      = matp(9);
xsigy0    = matp(10);


pr_yx = (young_y * pr_xy) / young_x;
pr_zy = (young_z * pr_yz) / young_y;
pr_zx = (young_z * pr_xz) / young_x;


xy_yx = pr_xy*pr_yx;
yz_zy = pr_yz*pr_zy;
zx_xz = pr_zx*pr_xz;
xyz   = 2*pr_xy*pr_yz*pr_zx;
E_xyz = young_x*young_y*young_z;

delta = (1 - (xy_yx) - (yz_zy) - (zx_xz) - (xyz)) / E_xyz;


C = zeros(6,6);
% ELastic stiffness matix (6*6)
C(1,1) = (1 -yz_zy) / (young_y*young_z*delta);
C(1,2) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta);
C(1,3) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta);
C(1,4) = 0;
C(1,5) = 0;
C(1,6) = 0;
C(2,1) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta);
C(2,2) = (1 -zx_xz) / (young_x*young_z*delta);
C(2,3) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta);
C(2,4) = 0;
C(2,5) = 0;
C(2,6) = 0;
C(3,1) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta);
C(3,2) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta);
C(3,3) = (1 -xy_yx) / (young_x*young_y*delta);
C(3,4) = 0;
C(3,5) = 0;
C(3,6) = 0;
C(4,1) = 0;
C(4,2) = 0;
C(4,3) = 0;
C(4,4) = g_xy;
C(4,5) = 0;
C(4,6) = 0;
C(5,1) = 0;
C(5,2) = 0;
C(5,3) = 0;
C(5,4) = 0;
C(5,5) = g_yz;
C(5,6) = 0;
C(6,1) = 0;
C(6,2) = 0;
C(6,3) = 0;
C(6,4) = 0;
C(6,5) = 0;
C(6,6) = g_xz;
C

a = [-1,0,0,0,0,-1000];

B = a*C*a'

B