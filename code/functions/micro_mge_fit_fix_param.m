function mfs = micro_mge_fit_fix_param(s_pa, xps, par_id, val)
if ~isequal(size(s_pa), size(xps.b1)); error('Stop'); end

unit_to_SI = [1e-9   1  1  1  1  1  1]; % MD, Kiso, Kaniso, Kinf_iso, Kinf_aniso, muK, k

t_ub      = [3  5  5   5   5  500];
t_lb      = [0  -5  -5  -5  -5   0];
t0_ub = [1  3  3  3  3  100 ];


switch par_id
    case "Kiso"
        id = 2;
    case "Kaniso"
        id = 3;
    case "Kiso_p"
        id = 4;
    case "Kaniso_p"
        id = 5;
    case "Kmu"
        id = 6;
    case "k"
        id = 7;
        t_ub      = [3  5  5   5   5 5 ];
        t_lb      = [0  -5  -5  -5  -5   -5];
        t0_ub = [1  3  3  3  3 3 ];
end

t0_lb = 0*t_ub+1e-9;

fun = @(x)get_target(x, id, val); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');
% t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
%
% %now trying to vary the initial conditions
Ntrials = 5;
RSS = inf;
for trial = 1:Ntrials
    t_0 = t0_lb + rand*(t0_ub-t0_lb);
    t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
    t = [t(1:(id-1)) val t(id+1:end)];
    tmp_signal = signal_micro_mge(t.*unit_to_SI, xps);
    tmp_RSS = sum((tmp_signal - s_pa).^2);
    if tmp_RSS < RSS; best_t = t; RSS = tmp_RSS; end
end
t = best_t;

% Objective
    function target = get_target(t_0, id, val)
        tmp_t0 = [t_0(1:(id-1)) val t_0(id+1:end)];
        m = tmp_t0.*unit_to_SI;
        
        signal = signal_micro_mge(m, xps);
        target = (signal)-(s_pa);
    end
%Evaluate final results
t = [t(1:(id-1)) val t(id+1:end)];
mfs.signal = get_target(t, id, val)+ s_pa;
mfs.MD = t(1);
mfs.Kiso = t(2);
mfs.Kaniso = t(3);
mfs.Kinf_iso = t(4);
mfs.Kinf_aniso = t(5);
mfs.muK = t(6);
mfs.k = t(7);
end
