function [out]=partitiontri(in)

if numel(in)==6
    
    % 6x1-matrix -> 4x1-matrix
    out=[in(4); 
         in(5); 
         in(6)];
    
    
elseif numel(in)==36

    % 6x6-matrix -> 4x4-matrix
    out=[  in(4,4) in(4,5) in(4,6);
           in(5,4) in(5,5) in(5,6);
           in(6,4) in(6,5) in(6,6)];
 
end