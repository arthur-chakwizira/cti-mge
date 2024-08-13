function signal = signal_kurtosis(m, xps)

b = xps.b;
%Model parameters
MD = m(1);
MK = m(2);

signal = exp(-b.*MD + (1/6)*MK*(MD^2)*b.^2);

end