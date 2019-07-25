function [cur_n_gam,gam_mtx,gam_val,gam_idx] = ...
   add_new_wavenum(new_wavenum,cur_n,cur_n_gam,gam_mtx,gam_val,gam_idx, ...
      wh_w_est,wh_w_inv,Gam_vec,w_est,s_power)
% add new candidate wavenumbers after one has been used in the
% approximation
d = size(new_wavenum,2);
cur_n_gam = cur_n_gam + 1; %one addition to the wavenumber pool
add_wavenum = new_wavenum; %start with one we just had
imp_k = max(wh_w_inv(new_wavenum > 0)); %find the least important coordinate that is nonzero
if isempty(imp_k)
   imp_k = 1;
end
wh_least_k = wh_w_est(imp_k);
add_wavenum(wh_least_k) = add_wavenum(wh_least_k) + 1; %add one to that value and add to pool
gam_mtx(cur_n_gam,:) = add_wavenum;
gam_val(cur_n_gam) = ... %compute gamma value for this new wavenumber
comp_wts_iter(Gam_vec,w_est,s_power,gam_mtx(cur_n_gam,:));
%gam_idx points to where the gammas are in order of largest to
%smallest starting with those after cur_n
wh_insert = find(gam_val(cur_n_gam) > gam_val(gam_idx(cur_n+1:cur_n_gam-1)),1,'first'); %find where this fits in the order
if ~isempty(wh_insert)
   gam_idx(cur_n+wh_insert:cur_n_gam) = [cur_n_gam; gam_idx(cur_n+wh_insert:cur_n_gam-1)];
else
   gam_idx(cur_n_gam) = cur_n_gam;
end
for k = imp_k+1:d
   %add additional interactions for all other 
   cur_n_gam = cur_n_gam + 1; %one addition to the wavenumber pool
   add_wavenum = new_wavenum; %start with one we just had
   add_wavenum(wh_w_est(k)) = 1; %add one where there was a zero
   gam_mtx(cur_n_gam,:) = add_wavenum;
   gam_val(cur_n_gam) = ...
    comp_wts_iter(Gam_vec,w_est,s_power,gam_mtx(cur_n_gam,:));
   wh_insert = find(gam_val(cur_n_gam) > gam_val(gam_idx(cur_n+1:cur_n_gam-1)),1,'first'); %find where this fits in the order
   if ~isempty(wh_insert)
      gam_idx(cur_n+wh_insert:cur_n_gam) = [cur_n_gam; gam_idx(cur_n+wh_insert:cur_n_gam-1)];
   else
      gam_idx(cur_n_gam) = cur_n_gam;
   end
end
