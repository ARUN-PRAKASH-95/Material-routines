function [out]=kelvinmt(in)
%kelvintm.m maps 4th order tensor to 6x6 matrix in Kelvin notation
% coded by: B. Kiefer 9 Feb 2012

sq2 = sqrt(2);

fak4=[  1   1   1 sq2 sq2 sq2;
        1   1   1 sq2 sq2 sq2;
        1   1   1 sq2 sq2 sq2;
      sq2 sq2 sq2   2   2   2;
      sq2 sq2 sq2   2   2   2;
      sq2 sq2 sq2   2   2   2;];
  
 nn=[ 1 4 6;
      4 2 5;
      6 5 3;]; 

temp=zeros(6);  
for i=1:6
    for j=1:6
        temp(i,j) = in(i,j)/fak4(i,j);
    end
end

out = zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
               out(i,j,k,l)=temp(nn(i,j),nn(k,l));
            end
        end
    end
end
