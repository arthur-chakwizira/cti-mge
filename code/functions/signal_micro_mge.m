
function signal = signal_micro_mge(m, xps)
MD = m(1);
Kiso = m(2);
Kaniso = m(3);
Kinf_iso = m(4);
Kinf_aniso = m(5);
muK = m(6);
k = m(7);

b1 = xps.b1;
b2 = xps.b2;

bdelta = xps.bdelta;
bmu = xps.bmu;

dt = xps.t(2)-xps.t(1);
c2 = -(b1+b2)*MD;

b2_k = 2*sum(xps.Q4_I.*exp(-k*xps.t), 2)*dt;
bdelta2_k = 2*sum(xps.Q4_A.*exp(-k*xps.t), 2)*dt./(b2_k+eps);

KT =   (Kiso.*b2_k + Kaniso.*b2_k.*bdelta2_k) + (b1+b2).^2.*(Kinf_iso + bdelta.^2*Kinf_aniso + bmu.*muK);
c4 = (1/6)*MD^2.*(KT);
signal = exp(c2 + c4);
end
