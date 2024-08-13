function signal = signal_cti(m, cti_xps)
    MD = m(1);
    KT = m(2);
    Kaniso = m(3);
    Kiso = m(4);
    b1 = cti_xps.b1;
    b2 = cti_xps.b2;
    theta = cti_xps.theta;
    signal = exp(-(b1+b2)*MD + (1/6)*(b1.^2+b2.^2)*MD^2*KT + (1/2)*b1.*b2.*( (cos(theta)).^2 )*MD^2*Kaniso  +  (1/6)*b1.*b2.*MD^2*(2*Kiso - Kaniso)  );
end
