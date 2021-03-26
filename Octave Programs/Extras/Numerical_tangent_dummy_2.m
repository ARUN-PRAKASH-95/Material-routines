if ttype = 0
    for i=1:6
        for j=1:6
            if eps(1) < epsilon_f
              A66(i,j) = (1+eta-damage)*C(i,j);
            else
              A66(i,j) = C_T(ii(i),jj(i),ii(j),jj(j));
            end  
        end
    end
elseif ttype = 1
    %perturbation, here very naiv
    hper=1e-8;
    %perturbation of the strain entries (here total strain, maybe this has to be modified)
    for ieps=1:1:length(eps6)
        epsper=eps6;
        epsper(ieps)=epsper(ieps)+hper;
        %recursiv call of your material routine with ttype=0 to avoid
        %endless loop
        %Calculate perturbed stress, sdv are not overwritten
        [sig6per,Adummy,sdvldummy]=vmises(epsper,sdvl,0);
        %Simple differential quotient
        A66_num(:,ieps)=(sig6per-sig6)/hper;
        
    end
    A66=A66_num;
end