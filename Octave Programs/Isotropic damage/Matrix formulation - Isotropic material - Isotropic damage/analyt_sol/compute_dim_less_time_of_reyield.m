function Phi = compute_dim_less_time_of_reyield(t_dst_d_T,q_el,strain_ampl,strain_T,epsilon_pl,alpha,mu,sigma_y0,H,h)
%
% compute dimensionless time of re-yielding from the evaluation of the
% von Mises yield function
%
tr_strain_T=sum(strain_T);
dev_strain_T=strain_T-tr_strain_T/3;
%
dev_strain_dst=dev_strain_T-strain_ampl*[2/3*(1-q_el);-1/3*(1-q_el);-1/3*(1-q_el)]*(t_dst_d_T-1);
xi=2*mu*(dev_strain_dst-epsilon_pl)-H*epsilon_pl;
norm_xi=norm(xi);
Phi=norm_xi-sqrt(2/3)*(sigma_y0+h*alpha);
end

