function [gam_val,w_est,gam_idx,f_hat_nm,errBd] = ...
   samp_sz_iter(four_coef,Gam_vec,w_vec,s_power,s_sum,gam_mtx,C,n0, ...
   samp_idx,nm_flg,w_flg,gam_val)
    
   % Computes sample size given desired error tolerance eps
   % - four_coef : true or estimated fourier coefficients
   % - Gam_vec : order weights
   % - w_vec : coordinate weights
   % - s_vec : smoothness weights
   % - gam_mtx : matrix of possible wavenumbers
   % - epsTol : desired error tolerance (don't use eps = machine epsilon)
   % - C : inflation factor
   % - n0 : initial samples per coordinate
   % - samp_idx : index of bases used for computing weights
   % - nm_flg : true if we know ||\hat{f}||_\gamma
   % - w_flg : true if we know product wts
   % - w_ini : initial w for optimization
   % - gam_val : gamma values computed from Gam_vec, w_vec, and s_vec
   
    
        
   % Split by flags
   if ~w_flg %Need to estimate the coordinate weights
      w_max = 1;
      [nBasis,d] = size(gam_mtx);

      %if samp_idx not provided, sample exact coefficients for 1-d weights (up to n0 smoothness)
      if isempty(samp_idx)
         samp_idx = find(sum(gam_mtx ~= 0,2) == 1 ... %just one nonzero element in row
            & max(gam_mtx,[],2) <= n0); %largest wavenumber is small enough
      else
         samp_idx(samp_idx==1) = []; %remove intercept (if in samp_idx)
      end
      nSamp_idx = size(samp_idx,1);
      [whCoord,~] = ind2sub([d,nSamp_idx],find(gam_mtx(samp_idx,:)' ~= 0));

      %Estimate product weights
      fpOvers = abs(four_coef(samp_idx)) ...
         ./ ((gam_mtx(sub2ind([nBasis d],samp_idx, whCoord))).^s_power)';
      uell = zeros(1,d);
      for ell = 1:d
         uell(ell) = max(fpOvers(whCoord == ell));
      end
      f_hat_pre = max(max(uell)/w_max); %smallest the norm can be
      w_est = uell/f_hat_pre; %smallest the w_l can be
      % samp_idx = [1; samp_idx]; %put back in the intercept
   else
      w_est = w_vec; %use inputted coordinate weights
   end
    
   %Compute gamma weights
   if ~exist('gam_val','var')
      gam_val = comp_wts_iter(Gam_vec,w_est,s_power,gam_mtx);
      %[~,gam_idx] = sort(gam_val,'descend');
      gam_idx = 1:length(gam_val);
   end
   
    
   if nm_flg %if we know norm...       
      % Compute (exact) f_hat from all fourier coefficients
      f_hat_nm = max(abs(four_coef) ./ gam_val); 
   else %we do not know the norm but bound it using Fourier coefficients indexed by samp_idx   
      % Approx. norm of f_hat
      f_hat_nm = C * max( abs(four_coef(samp_idx)) ./ gam_val(samp_idx)); 
   end
    
    % Set sample size
    tail_sum = prod(1 + w_est*s_sum) - sum(gam_val);
    errBd = f_hat_nm * tail_sum;

end
