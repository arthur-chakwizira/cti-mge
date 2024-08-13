function mfs = mge_1D_fit(s_pa, mge_xps)
%
unit_to_SI = [1e-9   1   1  ]; % MD, KT, k

t_ub      = [3.5       3        500  ];
t_lb      = [0         0         0   ];

% t_0 = [1         2         0.5           10         ];

fun = @(x)get_target(x); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');
% t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

Ntrials = 30;
RSS = inf;
for trial = 1:Ntrials
t_0 = t_lb + rand*(t_ub-t_lb);
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);
tmp_signal = signal_mge_1D(t.*unit_to_SI, mge_xps);
tmp_RSS = sum((tmp_signal - s_pa).^2);
if tmp_RSS < RSS; best_t = t; RSS = tmp_RSS; end
end
t = best_t;

%Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = signal_mge_1D(m, mge_xps);
               target = (signal)-(s_pa);
    end
%Evaluate final results
mfs.signal = get_target(t)+ s_pa;
mfs.MD = t(1);
mfs.KT = t(2);
mfs.k = t(3);
end
