function signal = signal_mge_1D(m, mge_xps)

b = mge_xps.b;
%Model parameters
MD = m(1);
KT = m(2);
k = m(3);

dt = mge_xps.t(2)-mge_xps.t(1);
h = 2*sum(mge_xps.Q4.*exp(-k*mge_xps.t), 2)*dt;

signal = exp(-b.*MD + (1/6)*b.^2*MD^2.*(KT.*h));
end