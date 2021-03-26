function [eps_t,eps_pl,sigma,alpha] = isochoric_tension_ref_sol(n_ampl,mat_param)
%
%
% material parameter
xE = mat_param(1);
xnu = mat_param(2);
sigY0 = mat_param(3);
H = mat_param(4);
h = mat_param(5);
mu = xE/(2*(1+xnu));
kappa = xE/(3*(1-2*xnu));
%
% strain amplitude for monotonically increasing strain
epsilon_bar_11=sigY0/(2*mu)*n_ampl;
%
eps_t=zeros(6,3);
eps_pl=zeros(6,3);
sigma=zeros(6,3);
alpha=zeros(6,3);
%
% strain history
eps_t(1:3,2)=[1;-0.5;-0.5]*2/3/n_ampl*epsilon_bar_11;
eps_t(1:3,3)=[1;-0.5;-0.5]*epsilon_bar_11;
%
% initial yield
sigma(:,2)=2*mu*eps_t(:,2);
%
% elastic-plastic loading
lambda=2*mu/(2*mu+H+2/3*h)*sqrt(3/2)*epsilon_bar_11*(1-2/3/n_ampl);
eps_pl(1:3,3)=[1;-0.5;-0.5]/sqrt(3/2)*lambda;
sigma(:,3)=2*mu*(eps_t(:,3)-eps_pl(:,3));
%
end