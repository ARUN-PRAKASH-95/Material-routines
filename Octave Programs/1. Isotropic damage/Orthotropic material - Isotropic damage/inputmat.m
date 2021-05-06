function [matp] = inputmat()
%inputmat.m define material parameters

matp(1)  = 55000e6;           % young_x
matp(2)  = 9500e6;           % young_y
matp(3)  = 9500e6;           % young_z
matp(4)  = 0.33;              % Poisson_xy
matp(5)  = 0.27;              % Poisson_yz
matp(6)  = 0.33;              % Poisson_xz
matp(7)  = 5500e6;            % Shear_mod_xy
matp(8)  = 3000e6;            % Shear_mod_yz
matp(9)  = 5500e6;            % Shear_mod_xz
matp(10) = 2500e6;            % Longitudinal tensile strength     


end
