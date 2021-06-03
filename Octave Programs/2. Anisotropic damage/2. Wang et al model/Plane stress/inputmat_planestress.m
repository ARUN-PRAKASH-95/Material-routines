function [matp] = inputmat_planestress()
%inputmat.m define material parameters

matp(1)  = 55.8e9;          % young_x
matp(2)  = 54.9e9;           % young_y
matp(3)  = 0.043;             % Poisson_xy
matp(4)  = 4.2e9;           % Shear_mod_xy
matp(5)  = 910.1e6;           % Longitudinal tensile strength     
matp(6)  = -710.2e6;            % Longitudinal compressive strength
matp(7)  = 772.2e6;              % Transverse tensile strength
matp(8)  = -703.3e6;             % Transverse compressive strength
matp(9)  = 131e6;              % Shear strength
matp(10) = 12.5e6;            % Longitudinal Tensile fracture energy
matp(11) = 12.5e6;            % Longitudinal Compressive fracture energy
matp(12) = 1e6;               % Transverse tensile fracture energy
matp(13) = 1e6;              % Transverse compressive fracture energy
matp(14) = 1e-3;          % Characteristic length

end
