function [matp] = inputmat()
%inputmat.m define material parameters

matp(1)  = 210.0e3;         % xE
matp(2)  = 0.33;            % xnu 
matp(3)  = 200.0;           % xsigy0
matp(4)  = 10;            % G_c
matp(5)  = 1**(1/3);         % Characteristic length

end
