function mfs = micro_mge_fit(s_pa, xps)
if ~isequal(size(s_pa), size(xps.b1)); error('Stop'); end

unit_to_SI = [1e-9   1  1  1 1  1  1]; % MD, Kiso, Kaniso, Kinf_iso, Kinf_aniso, muK, k

t_ub      = [3  5  5   5   5 5  500];
t_lb      = [0  -5  -5  -5  -5   -5    0];

t0_lb = 0*t_ub+1e-9;
t0_ub = [1  3  3  3  3 3  100 ];

fun = @(x)get_target(x); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');
% t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
% 
% %now trying to vary the initial conditions
Ntrials = 10;
RSS = inf;
for trial = 1:Ntrials
t_0 = t0_lb + rand*(t0_ub-t0_lb);
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
tmp_signal = signal_micro_mge(t.*unit_to_SI, xps);
tmp_RSS = sum((tmp_signal - s_pa).^2);
if tmp_RSS < RSS; best_t = t; RSS = tmp_RSS; end
end
t = best_t;

% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = signal_micro_mge(m, xps);
        target = (signal)-(s_pa);
    end
%Evaluate final results
mfs.signal = get_target(t)+ s_pa;
mfs.MD = t(1);
mfs.Kiso = t(2);
mfs.Kaniso = t(3);
mfs.Kinf_iso = t(4);
mfs.Kinf_aniso = t(5);
mfs.muK = t(6);
mfs.k = t(7);
end
