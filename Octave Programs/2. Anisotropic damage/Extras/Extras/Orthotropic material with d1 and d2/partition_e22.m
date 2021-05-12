function [out]=partition_e22(in)

if numel(in)==6
    
    % 6x1-matrix -> 5x1-matrix
    out=[in(1); 
         in(3); 
         in(4); 
         in(5); 
         in(6)];
    
elseif numel(in)==9
    
    % 3x3-matrix -> 5x1-matrix
    out=[in(2,2); 
         in(3,3); 
         in(1,2); 
         in(2,3); 
         in(1,3)];
    
elseif numel(in)==36

    % 6x6-matrix -> 5x5-matrix
    out=[in(1,1) in(1,3) in(1,4) in(1,5) in(1,6);
         in(3,1) in(3,3) in(3,4) in(3,5) in(3,6);
         in(4,1) in(4,3) in(4,4) in(4,5) in(4,6);
         in(5,1) in(5,3) in(5,4) in(5,5) in(5,6);
         in(6,1) in(6,3) in(6,4) in(6,5) in(6,6)];
 
end