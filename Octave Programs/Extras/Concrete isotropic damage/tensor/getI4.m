
function [I4] = getI4()

  I4  = zeros(3,3,3,3);
    for i = 1:3
      for j = 1:3
        for k = 1:3
          for l = 1:3
            if i == j
              I4(i,j,k,l) = 1;
            else
              I4(i,j,k,l) = 0;
            end
          end
        end
      end 
    end 

  
end