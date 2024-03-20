function mfs = mge_fit(s_pa, xps)
%Fits tensorial MGE to signals

if ~isequal(size(s_pa), size(xps.b1)); error('Stop'); end

unit_to_SI = [1e-9   1  1   1  1 1]; % MD, Kiso, Kaniso, k, Kinf_iso, Kinf_aniso

t_ub      = [3  5  5   300  5  5];
t_lb      = [0  -5  -5   0    0  0];

t0_ub = [1  3  3  20   3  3];
t0_lb = [0  0  0  0    0  0];

fun = @(x)get_target(x); %objective

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');

%vary the initial conditions
Ntrials = 20;
RSS = inf;
for trial = 1:Ntrials
t_0 = t0_lb + rand*(t0_ub-t0_lb);
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
tmp_signal = signal_mge(t.*unit_to_SI, xps);
tmp_RSS = sum((tmp_signal - s_pa).^2);
if tmp_RSS < RSS; best_t = t; RSS = tmp_RSS; end
end
t = best_t;

% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = signal_mge(m, xps);
        target = (signal)-(s_pa);
    end

%Evaluate final results
mfs.signal = get_target(t)+ s_pa;
mfs.MD = t(1);
mfs.Kiso = t(2);
mfs.Kaniso = t(3);
mfs.k = t(4);
mfs.Kinf_iso = t(5);
mfs.Kinf_aniso = t(6);
mfs.KT = mfs.Kiso + mfs.Kaniso + mfs.Kinf_iso + mfs.Kinf_aniso;
end
