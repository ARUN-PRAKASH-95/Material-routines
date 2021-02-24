clear all;
% *************************************************************************
% test suite for three-dimensional von Mises plasticity 
% isochoric tension (fully-deformation controlled test, ramp function)
% *************************************************************************
%
% strain amplitude in terms of multiple of normalized yield stress
% strain_ampl=sigma_y0/2/mu*n_ampl
n_ampl=2;      % n_ampl>1
%
mat_param = inputmat();
xE = mat_param(1); xnu = mat_param(2); sigma_y0 = mat_param(3);
mu = xE/(2*(1+xnu));
%
% define loading
% 1: linear ramping of load
ltype=1;
dt=0.1;
if ltype==1
    t=[0 10];
    lam=[0 n_ampl*sigma_y0/2/mu];
end
%
% computation of tangent moduli
ttype = 0; % 0: analytical
%
% path to auxiliary functions
addpath('tensor/');
%% computation
%
% prescribed load/time step
% start and end-time of loading, time-scale, no. of steps
ta=t(1);
te=t(end);
time=ta:dt:te;
steps=size(time,2)-1;
e11=loading(ltype,dt,t,lam);
%
% initialize internal variables
sdv=zeros(13,steps);
%
% initialise quantities for post-processing
s11=zeros(1,steps); s22=zeros(1,steps); s33=zeros(1,steps);
%
% initialize waitbar
wb=waitbar(0,'computation in progress...');

for n=1:steps
%
% display waitbar
    waitbar(n/steps);
% display current time step
    disp(['n = ', num2str(n)]);
%
    epsilon=[1;-0.5;-0.5;0;0;0]*e11(n+1);
%
% constitutive law: algorithmic stresses and moduli 
    [s,A,sdvup]=vmises(epsilon,sdv(:,n),ttype);
%
% update of internal variables after obtaining convergence
    sdv(:,n+1) = sdvup;
%    
% store quantities for post-processing
    s11(n+1)=s(1); s22(n+1)=s(2); s33(n+1)=s(3);
%
end
close(wb);
%
% reference solution
[eps_t,eps_pl,sigma,alpha] = isochoric_tension_ref_sol(n_ampl,mat_param);
%
%% visualization
figure(1)
clf;
hold on
plot(e11,s11,'r-')
plot(e11,s22,'b-')
plot(e11,s33,'g-')
%
plot(eps_t(1,:),sigma(1,:),'ko-.')
plot(eps_t(1,:),sigma(2,:),'ko-.')
plot(eps_t(1,:),sigma(3,:),'ko-.')
%
xlabel('\epsilon_{11}','FontSize',12)
ylabel('\sigma_{11}, \sigma_{22}, \sigma_{33} in MPa','FontSize',12)
legend('\sigma_{11}','\sigma_{22}','\sigma_{33}','Ref.','Location','East')
%
figure(2)
clf;
hold on
plot(e11,sdv(1,:),'r-')
plot(e11,sdv(2,:),'b-')
plot(e11,sdv(3,:),'g-')
%plot(e11,sdv(13,:),'m-')
%
plot(eps_t(1,:),eps_pl(1,:),'ko-.')
plot(eps_t(1,:),eps_pl(2,:),'ko-.')
plot(eps_t(1,:),eps_pl(3,:),'ko-.')
%
xlabel('\epsilon_{11}','FontSize',12)
ylabel('\epsilon_{11}^{pl}, \epsilon_{22}^{pl}, \epsilon_{33}^{pl}','FontSize',12)
legend('\epsilon_{11}^{pl}','\epsilon_{22}^{pl}','\epsilon_{33}^{pl}','Ref.','Location','NorthWest')