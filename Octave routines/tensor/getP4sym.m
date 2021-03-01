% returns 4th order symmetric deviatoric projection tensor
% [P4sym] = getP4sym()
% b. kiefer, iofm

function [P4sym] = getP4sym()

	% 2nd order identity
	I = eye(3);
	
	% initialize and assemble 4th order symmetric identity tensor
	I4sym = zeros(3,3,3,3);
	for i = 1:3
		for j = 1:3
			for k = 1:3
				for l = 1:3
					I4sym(i,j,k,l) = 1/2*(I(i,k)*I(j,l)+I(i,l)*I(j,k)); 
				end
			end
		end 
	end % for
    
    % initialize and assemble 4th order symmetric identity tensor
	P4sym = zeros(3,3,3,3);
	for i = 1:3
		for j = 1:3
			for k = 1:3
				for l = 1:3
                    P4sym(i,j,k,l) = I4sym(i,j,k,l)-1/3*I(i,j)*I(k,l); 
				end
			end
		end 
	end % for
	
end % function