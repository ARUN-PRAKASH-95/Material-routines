function [ Extendedvector ] = CreateRow( Vector,VariationU )
%extends the displacementvector because of calculating the displacement by
%deleting rows and columns 

%read length of the whole displacement vector
Len     =length(Vector)+length(VariationU);

%set place 1
e   =1;

%set empty vector
V   =[];

%assembly the vector
for i = 1:Len
    
    control=1;
    
    %look wheter the row was deleted
    for j = 1:length(VariationU)
        
        if i == VariationU(j)
            
            control     =0;
            
        end
        
    end
    
    %if rows was deleted, set value to zero   
    if control == 0
         
        V(i)    =0;
         
    else
         
        V(i)    =Vector(e);
        
        %next place
        e       =e+1;
        
    end
        
end

%transpose the vector
Extendedvector=V';

end

