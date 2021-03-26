function [out]=kelvintm(in)

sq2 = sqrt(2);

if numel(in)==9
    
    % second order tensor -> 6x1-matrix
    out=[in(1,1); in(2,2); in(3,3); sq2*in(1,2); sq2*in(2,3); sq2*in(1,3)];
    
elseif numel(in)==81

    % Kelvin-notation
    % forth order tensor -> 6x6-matrix
    out=[    in(1,1,1,1)     in(1,1,2,2)     in(1,1,3,3) sq2*in(1,1,1,2) sq2*in(1,1,2,3) sq2*in(1,1,1,3);
             in(2,2,1,1)     in(2,2,2,2)     in(2,2,3,3) sq2*in(2,2,1,2) sq2*in(2,2,2,3) sq2*in(2,2,1,3);
             in(3,3,1,1)     in(3,3,2,2)     in(3,3,3,3) sq2*in(3,3,1,2) sq2*in(3,3,2,3) sq2*in(3,3,1,3);
         sq2*in(1,2,1,1) sq2*in(1,2,2,2) sq2*in(1,2,3,3)   1*in(1,2,1,2)   1*in(1,2,2,3)   1*in(1,2,1,3);
         sq2*in(2,3,1,1) sq2*in(2,3,2,2) sq2*in(2,3,3,3)   1*in(2,3,1,2)   1*in(2,3,2,3)   1*in(2,3,1,3);
         sq2*in(1,3,1,1) sq2*in(1,3,2,2) sq2*in(1,3,3,3)   1*in(1,3,1,2)   1*in(1,3,2,3)   1*in(1,3,1,3)];
 
end