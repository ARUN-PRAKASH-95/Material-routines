function [matp] = inputmat()
%inputmat.m define material parameters

matp(1)  = 38500e6;           % young_x
matp(2)  = 16500e6;           % young_y
matp(3)  = 16500e6;           % young_z
matp(4)  = 0.27;              % Poisson_xy
matp(5)  = 0.40;              % Poisson_yz
matp(6)  = 0.27;              % Poisson_xz
matp(7)  = 4700e6;            % Shear_mod_xy
matp(8)  = 4700e6;            % Shear_mod_yz
matp(9)  = 4700e6;            % Shear_mod_xz
matp(10) = 1250e6;            % Longitudinal tensile strength     


end
