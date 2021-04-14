function [matp] = inputmat();



matp(1)  = 146.9e9;           % young_x
matp(2)  = 11.4e9;            % young_y
matp(3)  = 11.4e9;            % young_z
matp(4)  = 0.3;              % Poisson_xy
matp(5)  = 0.4;              % Poisson_yz
matp(6)  = 0.3;              % Poisson_xz
matp(7)  = 6.18e9;            % Shear_mod_xy
matp(8)  = 6.18e9;            % Shear_mod_yz
matp(9)  = 6.18e9;            % Shear_mod_xz

matp(10) = 1370.6e6;            % Longitudinal tensile strength
matp(11) = -1379e6;            % Longitudinal compressive strength
matp(12) = 66.5e6;              % Transverse tensile strength
matp(13) = -268.2e6;             % Transverse compressive strength
matp(14) = 133.8e6;              % Shear strength
matp(15) = 12.5e6;   %1000000000;            % Longitudinal Tensile fracture energy
matp(16) = 12.5e6;  %1000000000;            % Longitudinal Compressive fracture energy
matp(17) = 1e6;    %1000000000;               % Transverse tensile fracture energy
matp(18) = 1e6;    %1000000000;               % Transverse compressive fracture energy
matp(19) = 1**(1/3);          % Characteristic length

end
