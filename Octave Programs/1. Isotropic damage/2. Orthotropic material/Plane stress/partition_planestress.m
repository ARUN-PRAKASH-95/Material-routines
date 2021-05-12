function [out]=partition_planestress(in)

if numel(in)==3
    
    % 6x1-matrix -> 5x1-matrix
    out=[in(2); 
         in(3);]; 
         
    
    
elseif numel(in)==9

    % 6x6-matrix -> 5x5-matrix
    out=[in(2,2) in(2,3) 
         in(3,2) in(3,3)];

 
end