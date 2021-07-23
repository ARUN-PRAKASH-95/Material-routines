function [ Reducedmatrix ] = DeleteColumn( Matrix,VariationU )
%gives a column reduced Matrix

%read size of matrix 
Size    =size(Matrix);

%initiate new matrix
K1      =[];

%assembly new matrix
for append = 1:Size(2)
    
    Assembly = 1;
    
    %look, wheter column should be deleted
    for delete = 1:length(VariationU)
        
        if VariationU(delete) == append
            Assembly    =0;
           
        end
    end
    
    %assembly the column to the rest
    if Assembly == 1
        
        K1  =horzcat(K1,Matrix(:,append));
            
    end
end
   Reducedmatrix    =K1;
   
end


