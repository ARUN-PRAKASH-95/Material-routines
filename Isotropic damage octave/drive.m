%==========================================================================
%   drive.m
%
%   driver for uni-axial isothermal stress tests 
% 
%==========================================================================
% 
%   coded by: B. Kiefer 15 Nov 2011
%   analytical solution added: S. Prueger 10 Dec 2018
%
%   comments: --> numerical tangent moduli computation for ttype = 1
%
%==========================================================================
clear all
close all
clc
addpath('tensor/'); 
addpath('analyt_sol/'); 
%format long;
%
% strain amplitude in terms of multiple of normalized yield stress
% strain_ampl=sigma_y0/2/mu*n_ampl/(1-q_el)
n_ampl=3;      % n_ampl>1
%
mat_param = inputmat();
xE = mat_param(1); xnu = mat_param(2); sigma_y0 = mat_param(3);
mu = xE/(2*(1+xnu)); kappa = xE/(3*(1-2*xnu));
q_el=-0.5*(kappa-2/3*mu)/(kappa+1/3*mu);
%
%--------------------------------------------------------------------------
%   define loading
%--------------------------------------------------------------------------
% 1: linear ramping of load
% 2: half cycle with linear load change
% 3: linear loading and unloading, start and end point different
% 4: full cycle with linear load change
% 5: two cycles
ltype=3;
if ltype==1
    t=[0 10];
    lam=[0 0.001];
elseif ltype==2
    t=[0 5 10];
    lam=[0 1.59155e-3];
elseif ltype==3
    t=[0 5 10];
    lam=[0 0.0015 0.0000];
elseif ltype==4
    t=[0 2.5 7.5 10];
    lam=[0 n_ampl*sigma_y0/2/mu/(1-q_el) -n_ampl*sigma_y0/2/mu/(1-q_el)];
elseif ltype==5
    t=[0 2.5 7.5 12.5 17.5 22.5 25.0];
    lam=[0 0.005 -0.005];
end

%--------------------------------------------------------------------------

% prescribed load/time step
dt=0.1;
% start and end-time of loading, time-scale, no. of steps
ta=t(1);
te=t(end);
time=ta:dt:te;
steps=size(time,2)-1;
e11=loading(ltype,dt,t,lam);

%-------------------------------------------------------------------------

% initialize strains, temperature and internal variables
epsbar=zeros(5,1);
sdv = zeros(1,steps);
sdv(1,1)=0; 

% initialise quantities for post-processing
s11=zeros(1,steps);
eps22=zeros(1,steps); eps33=zeros(1,steps);

% tolerance and maximum no. of iterations for Newton iteration
tol=1e-10;
maxit=100;
ttype = 0; % 0: analytical, 1: numerical tangent moduli computation

% initialize waitbar
wb=waitbar(0,'computation in progress...');

for n=1:steps
    
    % display waitbar
    waitbar(n/steps)
    
    % display current time step
    disp(['n = ', num2str(n)])
    
    % initialize rnorm in order to enter the while-loop and set no. of
    % iterations to zero
    sbar=ones(5,1);
    iter=0;
    
    while norm(sbar) > tol % check convergence
        iter = iter+1;
        if iter > maxit
            error(['No convergence after ', num2str(maxit), ... 
                ' global iterations'])
        end
        
        % 1.) total deformation
         epsilon(1,1) = e11(n+1);
         epsilon(1,1)
         epsilon(2:6,1) = epsbar;
         
        % 2.) constitutive law: algorithmic stresses and moduli 
        [s,A,sdvup]=vmises(epsilon,sdv(:,n),ttype);
        
        % 3.) partitioning
        sbar=partition(s);
        Abar=partition(A);

        % 4.) update of lateral deformations
        epsbar=epsbar-Abar\sbar;

        % display convergence
        disp(['|sbar| = ', num2str(norm(sbar))])

    end % while
    
    
    % update of internal variables after obtaining convergence
    sdv(:,n+1) = sdvup;
    
    % store quantities for post-processing
    s11(n+1)=s(1);
    eps22(n+1)=epsilon(2); eps33(n+1)=epsilon(3); 

end % for
fprintf('********************\n')
close(wb)
%

figure(1)
%subplot(2,1,1)
plot(time,e11)
xlabel('time')
ylabel('e11')


% plot stress-strain response
figure(2);
 %subplot(2,1,2)
 %plot(e11,s11, paperX,paperY, '-.rx')
plot(e11,s11,'r-')
legend('\sigma_{11}','Ref.','Location','NorthWest')
xlabel('eps11')
ylabel('sig11')
 

 
 % plot lateral strains
figure(3)
plot(e11,eps22,'b-')
hold on
plot(e11,eps33,'g-')
legend('\epsilon_{22}','\epsilon_{33}','Ref.','Location','NorthEast')
xlabel('eps11')
ylabel('eps22, eps33')


% plot damage
figure(4)
plot(e11,sdv(1,:),'r-')
legend('damage','Ref.','Location','West')
xlabel('eps11')
ylabel('damage')