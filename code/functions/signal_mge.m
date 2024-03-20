function signal = signal_mge(m, xps)
%Returns signal according to tensorial MGE theory
    MD = m(1);
    Kiso = m(2);
    Kaniso = m(3);
    k = m(4);
    Kinf_iso = m(5);
    Kinf_aniso = m(6);
    
    b1 = xps.b1;
    b2 = xps.b2;
    bdelta = xps.bdelta;
    
    dt = xps.t(2)-xps.t(1);    
    c2 = -(b1+b2)*MD;

    b2_k = 2*sum(xps.Q4_I.*exp(-k*xps.t), 2)*dt;
    
    bdelta2_k = 2*sum(xps.Q4_A.*exp(-k*xps.t), 2)*dt./(b2_k+eps);
    
    KT =   (Kiso.*b2_k + Kaniso.*b2_k.*bdelta2_k) + (b1+b2).^2*Kinf_iso + bdelta.^2.*(b1+b2).^2*Kinf_aniso;

    c4 = (1/6)*MD^2.*(KT);
    signal = exp(c2 + c4);
end
