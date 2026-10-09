function power = benchmark_mean_metabolic_power(R)
% Integrate muscle metabolic power over the complete periodic gait cycle.
t = R.time.mesh_GC(:);
rates = R.metabolics.Bhargava2004.Edot_gait;
if size(rates,1) == numel(t)-1
    % PredSim omits the repeated endpoint in gait-cycle result arrays.
    rates(end+1,:) = rates(1,:);
end
assert(size(rates,1)==numel(t) && numel(t)>1 && all(diff(t)>0),...
    'PredSim:InvalidMetabolicTime','Metabolic samples must match the increasing gait-cycle mesh.');
power = sum(trapz(t,rates,1))/(t(end)-t(1));
end
