function mfs = cti_fit(s_pa, cti_xps)
%
if ~isequal(size(s_pa), size(cti_xps.b1)); error('Stop'); end
unit_to_SI = [1e-9   1  1   1]; % MD, KT, Kaniso, Kiso

t_ub      = [5   5   5    5];
t_lb      = [0   -5   0   -5];

% t_0 = [1   1    1   1];

fun = @(x)get_target(x); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');
% t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

%now trying to vary the initial conditions
Ntrials = 30;
RSS = inf;
for trial = 1:Ntrials
t_0 = t_lb + rand*(t_ub-t_lb);
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
tmp_signal = signal_cti(t.*unit_to_SI, cti_xps);
tmp_RSS = sum((tmp_signal - s_pa).^2);
if tmp_RSS < RSS; best_t = t; RSS = tmp_RSS; end
end
t = best_t;

% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = signal_cti(m, cti_xps);
        target = (signal)-(s_pa);
    end
%Evaluate final results
mfs.signal = get_target(t)+ s_pa;
mfs.MD = t(1);
mfs.KT = t(2);
mfs.Kaniso = t(3);
mfs.Kiso = t(4);

end