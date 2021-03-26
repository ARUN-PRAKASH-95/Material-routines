function [stress,strain,epsilon_plast] = analyt_sol_vM_plast_strain_contr_uniaxial_stress(T,T_sim,strain_ampl,n_mult)
% *************************************************************************
% analytical solution for von Mises plasticity with linear isotropic and
% kinematic hardening subject to uniaxial stress state (under
% strain-control)
% *************************************************************************
%
% material parameters
mat_param = inputmat();
xE = mat_param(1); xnu = mat_param(2);
sigma_y0 = mat_param(3);
H = mat_param(4);
h = mat_param(5);
mu = xE/(2*(1+xnu));
kappa = xE/(3*(1-2*xnu));
%
%% computations
%
% lateral contraction for elastic and elastic-plastic loading
q_el=-0.5*(kappa-2/3*mu)/(kappa+1/3*mu);
hard_mod=2*mu/(2*mu+H+2/3*h);
q_ep=-0.5*(kappa-2/3*mu+2/3*mu*hard_mod)/(kappa+1/3*mu-1/3*mu*hard_mod);
%
% strain amplitude for triangle wave function
time_prop = loading_history(T,T_sim,strain_ampl);
%
% initialize stress and internal variables
stress=zeros(3,5);
strain=zeros(3,5);
epsilon_plast=zeros(3,5);
alpha_iso_hard=zeros(5,1);
back_stress_beta=zeros(3,5);
lagr_mult=zeros(5,1);
%
% -------------------------------------------------------------------------
% initial yield (starting from alpha=0, epsilon_pl=0) t^{*}
% -------------------------------------------------------------------------
t_st=1/n_mult*T; % /(1-q_el)
epsilon_11=comput_strain_prscr_comp(t_st,time_prop);
epsilon_strain=[epsilon_11; q_el*epsilon_11; q_el*epsilon_11];
stress_temp = linear_isotr_elasticity_vol_isochoric_split(kappa,mu,epsilon_strain,epsilon_plast(:,2));
%
stress(:,2)=stress_temp;
strain(:,2)=epsilon_strain;
%
% -------------------------------------------------------------------------
% integration of evolution eq. for elastic plastic loading (t_st < t < T)
% -------------------------------------------------------------------------
lagr_mult_t=2*mu*(1-q_ep)/(2*mu+H+2/3*h)*sqrt(2/3)*(1-1/n_mult)*strain_ampl;
epsilon_pl_t=lagr_mult_t*sign(stress(1,2))*[2/3;-1/3;-1/3]/sqrt(2/3);
%
epsilon_11=comput_strain_prscr_comp(T,time_prop);
epsilon_strain=[epsilon_11; strain(2,2)+q_ep*(1-1/n_mult)*strain_ampl; strain(3,2)+q_ep*(1-1/n_mult)*strain_ampl];
%
stress_temp = linear_isotr_elasticity_vol_isochoric_split(kappa,mu,epsilon_strain,epsilon_pl_t);
%
stress(:,3)=stress_temp;
strain(:,3)=epsilon_strain;
epsilon_plast(:,3)=epsilon_pl_t;
alpha_iso_hard(3)=sqrt(2/3)*lagr_mult_t;
back_stress_beta(:,3)=H*epsilon_pl_t;
lagr_mult(3)=lagr_mult_t;
%
% -------------------------------------------------------------------------
% re-yielding under compressive loading t^{**}
% -------------------------------------------------------------------------
x0=[1+sqrt(eps),3];
t_dst_d_T=fzero(@(t_dst_d_T) compute_dim_less_time_of_reyield(t_dst_d_T,q_el,strain_ampl,strain(:,3),epsilon_plast(:,3),alpha_iso_hard(3),mu,sigma_y0,H,h),x0);
% t_temp=1:(3-1)/(99):3;
% residual=zeros(100,1);
% for i=1:100
%     residual(i)=compute_dim_less_time_of_reyield(t_temp(i),q_el,strain_ampl,strain(:,3),epsilon_plast(:,3),alpha_iso_hard(3),mu,sigma_y0,H,h);
% end
epsilon_strain=strain(:,3)-strain_ampl*[1;q_el;q_el]*(t_dst_d_T-1);
stress_temp = linear_isotr_elasticity_vol_isochoric_split(kappa,mu,epsilon_strain,epsilon_plast(:,3));
%
stress(:,4)=stress_temp;
strain(:,4)=epsilon_strain;
epsilon_plast(:,4)=epsilon_plast(:,3);
alpha_iso_hard(4)=alpha_iso_hard(3);
back_stress_beta(:,4)=back_stress_beta(:,3);
lagr_mult(4)=lagr_mult(3);
%
% -------------------------------------------------------------------------
% integration of evolution eq. for elastic plastic loading (t_dst < t < 3T)
% -------------------------------------------------------------------------
%
delta_lambda=2*mu*(1-q_ep)/(2*mu+H+2/3*h)*sqrt(2/3)*(T_sim/T-t_dst_d_T)*strain_ampl;
lagr_mult_t=lagr_mult(4)+delta_lambda;
epsilon_pl_t=epsilon_plast(:,4)+delta_lambda*sign(stress(1,4)-back_stress_beta(1,4))/sqrt(2/3)*[2/3;-1/3;-1/3];
epsilon_strain=strain(:,4)-strain_ampl*[1;q_ep;q_ep]*(T_sim/T-t_dst_d_T);
stress_temp = linear_isotr_elasticity_vol_isochoric_split(kappa,mu,epsilon_strain,epsilon_pl_t);
%
stress(:,5)=stress_temp;
strain(:,5)=epsilon_strain;
epsilon_plast(:,5)=epsilon_pl_t;
alpha_iso_hard(5)=alpha_iso_hard(4)+sqrt(2/3)*delta_lambda;
back_stress_beta(:,5)=H*epsilon_pl_t;
lagr_mult(5)=lagr_mult_t;
%
end