function epsilon = comput_strain_prscr_comp(time,time_prop)
%
% computes the value of the strain
%
ind=find(time_prop(:,1)>time);
if isempty(ind)
    if abs(time-time_prop(end,1))<sqrt(eps)
        epsilon=time_prop(end,2);
    else
        error('Current time is outside time range specified in "loading_conditions_deformation_control.m"')
    end;
else
    t_para_temp=(time-time_prop(ind(1)-1,1))/(time_prop(ind(1),1)-time_prop(ind(1)-1,1));
    epsilon=(time_prop(ind(1),2)-time_prop(ind(1)-1,2))*t_para_temp+time_prop(ind(1)-1,2);
end
%
end