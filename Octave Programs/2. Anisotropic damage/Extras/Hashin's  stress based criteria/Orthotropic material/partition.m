function [out]=partition(in)

if numel(in)==3
    
    % 3x1-matrix -> 2x1-matrix
    out = [in(2); 
           in(3);];
    
elseif numel(in)==9

    % 3x3-matrix -> 2x2-matrix
    out=[in(2,2) in(2,3);
         in(3,2) in(3,3);];
 
end