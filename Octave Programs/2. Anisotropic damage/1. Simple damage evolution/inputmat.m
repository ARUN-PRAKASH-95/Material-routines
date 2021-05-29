function [matp] = inputmat();



matp(1)  = 55000e6;           % young_x
matp(2)  = 9500e6;            % young_y
matp(3)  = 9500e6;            % young_z
matp(4)  = 0.33;              % Poisson_xy
matp(5)  = 0.27;              % Poisson_yz
matp(6)  = 0.33;              % Poisson_xz
matp(7)  = 5500e6;            % Shear_mod_xy
matp(8)  = 3000e6;            % Shear_mod_yz
matp(9)  = 5500e6;            % Shear_mod_xz

matp(10) = 2500e6;            % Longitudinal tensile strength
matp(11) = -2000e6;            % Longitudinal compressive strength
matp(12) = 50e6;              % Transverse tensile strength
matp(13) = -200e6;             % Transverse compressive strength
matp(14) = 50e6;              % Shear strength
matp(15) = 12.5e3;            % Longitudinal Tensile fracture energy
matp(16) = 12.5e3;            % Longitudinal Compressive fracture energy
matp(17) = 1e3;               % Transverse tensile fracture energy
matp(18) = 1e3;              % Transverse compressive fracture energy
matp(19) = 1;          % Characteristic length

end
