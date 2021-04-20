function stress = linear_isotr_elasticity_vol_isochoric_split(kappa,mu,epsilon_strain,epsilon_pl)
%
% isotropic, linear elasticity with volumetric / isochoric split for small
% strains (Hooke's law)
%
tr_epsilon=sum(epsilon_strain);
dev_epsilon_strain=epsilon_strain-tr_epsilon/3;
stress=kappa*tr_epsilon*[1;1;1]+2*mu*(dev_epsilon_strain-epsilon_pl);
end

