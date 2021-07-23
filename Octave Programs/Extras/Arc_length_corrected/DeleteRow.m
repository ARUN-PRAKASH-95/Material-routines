function [ Reducedmatrix ] = DeleteRow( Matrix,VariationU )
%gives a row reduced Matrix

%read size of matrix 
Size    =size(Matrix);

%initiate new matrix
Reducedmatrix  =[];

%assembly new matrix
for append = 1:Size(1)
    
    Assembly    =1;
    
    %look, wheter row should be deleted
    for delete = 1:length(VariationU)
        
        if VariationU(delete) == append
            
            Assembly    =0;
            
        end
        
    end
        
    %assembly the row to the rest
    if Assembly == 1
        
        Reducedmatrix  =vertcat(Reducedmatrix,Matrix(append,:));
        
    end
    
end

end

