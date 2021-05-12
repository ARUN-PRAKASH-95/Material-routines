function [matp] = inputmat_planestress()
%inputmat.m define material parameters

matp(1)  = 55000e6;          % young_x
matp(2)  = 9500e6;           % young_y
matp(3)  = 0.33;             % Poisson_xy
matp(4)  = 5500e6;           % Shear_mod_xy
matp(5)  = 2500e6;           % Longitudinal tensile strength     
matp(6)  = -2000e6;            % Longitudinal compressive strength
matp(7)  = 50e6;              % Transverse tensile strength
matp(8)  = -150e6;             % Transverse compressive strength
matp(9)  = 50e6;              % Shear strength
matp(10) = 12.5e3;            % Longitudinal Tensile fracture energy
matp(11) = 12.5e3;            % Longitudinal Compressive fracture energy
matp(12) = 1e3;               % Transverse tensile fracture energy
matp(13) = 1e3;              % Transverse compressive fracture energy
matp(14) = 1**(1/3);          % Characteristic length

end
