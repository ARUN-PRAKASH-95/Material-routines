function time_prop = loading_history(T,T_sim,ampl)
%
% load history in terms of applied strain
%
% time_prop.........matrix containing in each row (time, prop-factor) value
%
time_prop=[[0.0,0.0];
           [T,ampl];
           [T_sim,-ampl]];
%
end

