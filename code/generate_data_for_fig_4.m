
%Generate the data needed to make Fig.4

base_fn = fileparts(fileparts(mfilename('fullpath')));  %main folder containing all subfolders


prot_fn = base_fn + "/protocols/protocol_cti_mge.mat";
load(prot_fn, "xps")

xps.mge_s_ind = xps.mge_s_ind & (xps.b <= 0.5e9);


substrate_names = ["mge_iso", "mge_iso_aniso", "spheres_d6_regular_st"];

for c_sn = 1:numel(substrate_names)
    
    substrate_name = substrate_names(c_sn);
    
    true_k = [10 20 30 40 50];% 60 80 100];
    k_string = ""; for i = 1:numel(true_k); k_string(i) = num2str(true_k(i)); end
    pa_sig_fn = base_fn + "/signals/" + "signal_"+ substrate_name + "_k" + k_string + "_pa_v11";
    
    if strcmp(substrate_name, "mge_iso")
        true_Kiso = true_k*0 + 1;
        true_Kinf_iso = true_k*0;
        true_Kaniso = true_k*0;
        true_Kinf_aniso = true_k*0;
    end
    if strcmp(substrate_name, "mge_iso_aniso")
        true_Kiso = true_k*0 + 0.33;
        true_Kinf_iso = true_k*0;
        true_Kaniso = true_k*0 + 0.33;
        true_Kinf_aniso = true_k*0 + 0.33;
    end
    
    if strcmp(substrate_name, "spheres_d6_regular_st")
        true_Kiso = true_k*0 + 1.8;
        true_Kinf_iso = true_k*0 + 0.15;
        true_Kaniso = true_k*0;
        true_Kinf_aniso = true_k*0;
    end
    
    
    %now fit tensorial MGE to powdered signals
    
    kfit_tensorial = true_k*0;
    N_samples = 100;
    sn_r = 200;
    MD = zeros(numel(true_k), N_samples);
    Kiso = zeros(numel(true_k), N_samples);
    Kaniso = zeros(numel(true_k), N_samples);
    k = zeros(numel(true_k), N_samples);
    Kinf_iso = zeros(numel(true_k), N_samples);
    Kinf_aniso = zeros(numel(true_k), N_samples);
    
    for c_k = 1:numel(true_k)
        load(pa_sig_fn(c_k), 's_pa');
        %add noise
        s_pa_noisy = add_noise_to_signal(s_pa, sn_r, N_samples);
        parfor c_s = 1:N_samples
            tmp_s_pa = s_pa_noisy(c_s, :)';
            mfs_tensor = mge_fit(tmp_s_pa, xps); %Tensorial MGE
            MD(c_k, c_s) = mfs_tensor.MD;
            Kiso(c_k, c_s) = mfs_tensor.Kiso;
            Kaniso(c_k, c_s) = mfs_tensor.Kaniso;
            k(c_k, c_s) = mfs_tensor.k;
            Kinf_iso(c_k, c_s) = mfs_tensor.Kinf_iso;
            Kinf_aniso(c_k, c_s) = mfs_tensor.Kinf_aniso;
        end
        disp("Done fitting " + num2str(c_k) + " of " + num2str(numel(true_k)))
    end
    
    save(base_fn + "/fit_results/tensorial_MGE_fit_" + substrate_name, "MD", "Kiso", "Kaniso", "k", "Kinf_iso", "Kinf_aniso", "N_samples", "sn_r", "true_k",...
        "true_Kiso", "true_Kinf_iso", "true_Kaniso", "true_Kinf_aniso")
    
end

