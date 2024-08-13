function s_noise = add_noise_to_signal(signal, snr, N_samp)
%adds Rice distributed noise to signal
%at snr = snr
%generates N_samp realisations of noise
n_sig = numel(signal);
s_noise = zeros(N_samp, n_sig);
for n = 1:N_samp
    noise_real = randn(size(signal))*(1/snr);
    noise_imag = randn(size(signal))*(1/snr);
    tmp_s_noise = sqrt( (signal + noise_real).^2 + noise_imag.^2);
    s_noise(n,:) = tmp_s_noise;
end

end