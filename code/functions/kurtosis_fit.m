function mfs = kurtosis_fit(s_simulated, xps)
%
unit_to_SI = [1e-9   1]; % MD, MK

t_ub      = [3   5];
t_lb      = [0   -5];

fun = @(x)get_target(x); %objective

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');

% vary the initial conditions
Ntrials = 30;
RSS = inf;
for trial = 1:Ntrials
t_0 = t_lb + rand*(t_ub-t_lb);
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
tmp_signal = signal_kurtosis(t.*unit_to_SI, xps);
tmp_RSS = sum((tmp_signal - s_simulated).^2);
if tmp_RSS < RSS; best_t = t; RSS = tmp_RSS; end
end
% 
t = best_t;

%% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = signal_kurtosis(m, xps);
        target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.MD = t(1);
mfs.KT = t(2);

end