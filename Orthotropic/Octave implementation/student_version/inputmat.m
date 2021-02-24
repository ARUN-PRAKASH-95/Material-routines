function [matp] = inputmat()
%inputmat.m define material parameters

matp(1)  = 5e10;                % young_x
matp(2)  = 8e9;                 % young_y 
matp(3)  = 8e9;                 % young_z
matp(4)  = 0.3;                 % pr_xy 
matp(5)  = 0.4;                 % pr_yz
matp(6)  = 0.3;                 % pr_xz          
matp(7)  = 5e9;                 % g_xy
matp(8)  = 3.8462e9;            % g_yz
matp(9)  = 5e9;                 % g_xz

end
