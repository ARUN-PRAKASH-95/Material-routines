function [out]=partition_planestress_Bi(in)

if numel(in)==3
    
    % 6x1-matrix -> 5x1-matrix
    out=[in(3);]; 
         
    
    
elseif numel(in)==9

    % 6x6-matrix -> 5x5-matrix
    out=[in(3,3)];
          

 
end