function mfs = cti_fit_cascade(s_pa, cti_xps)
%CTI fit according to Henriques 2021
if ~isequal(size(s_pa), size(cti_xps.b1)); error('Stop'); end

sde_ind = (cti_xps.b2 == 0);
xps.b = cti_xps.b1(sde_ind) + cti_xps.b2(sde_ind);
mfs = kurtosis_fit(s_pa(sde_ind), xps);
MD = mfs.MD*1e-9;

%get uK
sde_ind = (cti_xps.b1 == max(cti_xps.b1 + cti_xps.b2))&(cti_xps.b2 ==0);
dde_par_ind = (cti_xps.b2 == 0.5*max(cti_xps.b1 + cti_xps.b2))&(cti_xps.b1 == 0.5*max(cti_xps.b1 + cti_xps.b2))&(cti_xps.theta == 0);
sde_dde_diff = log(s_pa(sde_ind)) - log(s_pa(dde_par_ind));
mfs.uK = 12*sde_dde_diff/((max(cti_xps.b1 + cti_xps.b2)).^2 * MD^2);

%get Kaniso
dde_orth_ind = (cti_xps.b2 == 0.5*max(cti_xps.b1 + cti_xps.b2))&(cti_xps.b1 == 0.5*max(cti_xps.b1 + cti_xps.b2))&(cti_xps.theta == pi/2);
dde_par_orth_diff = log(s_pa(dde_par_ind)) - log(s_pa(dde_orth_ind));
mfs.Kaniso = 2*dde_par_orth_diff/((max(cti_xps.b1 + cti_xps.b2))^2 * MD^2);

mfs.Kiso = mfs.KT - mfs.Kaniso - mfs.uK;

m = [MD, mfs.KT, mfs.Kaniso, mfs.Kiso];
mfs.signal = signal_cti(m, cti_xps);
end