%% Compute and plot the sample size required for different error criteria

%% Simulation settings
close all
clearvars

d_vec = [1 2 4];
nD = length(d_vec);
s_max = 30; %maximum smoothness for approximated function
s_power_vec = [3 4 8];
nS = length(s_power_vec);
num_eps = 25; %no. of errors
min_log10_eps = -3;
max_log10_eps = -1;
C = 1.5; % inflation factor
w_vec_all = 0.5.^[0 0 1:10];

% Error tolerances
eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance

%% Set initial design and estimate Fourier coefficients
% Set desired sampling indices J for initial design
nBLength = 1e3;


%% Run algorithm for different error tolerances
err_vec(length(eps_vec),1) = 0; %container for errors
n_vec(length(eps_vec),1) = 0; %container for sample sizes
gam_val(nBLength,1) = 0;
gam_idx(nBLength,1) = 0;
tail_sum = inf(nBLength,nD,nS);
n_req(nD,nS,num_eps) = 0;
gam_idx(1) = 1;

min_tol = min(eps_vec);

for ii = 1:nD
   d = d_vec(ii)
   Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
   w_vec = w_vec_all(1:d);
   gam_mtx = zeros(nBLength,d);
   gam_mtx(1,:) = zeros(1,d);
   wh_w_inv = zeros(1,d);
   for jj = 1:nS
      s_power = s_power_vec(jj)
      s_sum = zeta(s_power);
      cur_n = 1;
      cur_n_gam = 1;
      [~,wh_w_vec] = sort(w_vec,'descend'); %sort the coordinate weights
      wh_w_inv(wh_w_vec) = 1:d;
      gamma_sum = prod(1 + w_vec*s_sum);
      gam_val(cur_n_gam) = ... %compute gamma value for this new wavenumber
         comp_wts_iter(Gam_vec,w_vec,s_power,gam_mtx(cur_n_gam,:));
      tail_sum(ii,jj,cur_n) = gamma_sum - sum(gam_val(gam_idx(1:cur_n)));
      while (tail_sum(ii,jj,cur_n) > min_tol) %... if sample size not enough, add largest unobserved gamma
         new_idx = gam_idx(cur_n);
         new_wavenum = gam_mtx(new_idx,:);
         [cur_n_gam,gam_mtx,gam_val,gam_idx] = ...
            add_new_wavenum(new_wavenum,cur_n,cur_n_gam,gam_mtx,gam_val,gam_idx, ...
               wh_w_vec,wh_w_inv,Gam_vec,w_vec,s_power);
         cur_n = cur_n+1;
         tail_sum(ii,jj,cur_n) = gamma_sum - sum(gam_val(gam_idx(1:cur_n)));
      end
      for m = 1:num_eps
         n_req(ii,jj,m) = find(tail_sum(ii,jj,:) < eps_vec(m),1,'first');
      end
   end
end

InitializeDisplay
close all
figure
log10epsVec = log10(eps_vec);
hLeg{nD*nS,1} = 0;
hds(nD*nS,1) = 0;
for ii = 1:nD
   for jj = 1:nS
      h = loglog(eps_vec,reshape(n_req(ii,jj,:),[num_eps 1]),'.', ...
         'marker',markerSequence{ii},'color',colorSequence{jj},'markerSize',10); 
      hold on
      hds((ii-1)*nS+jj) = h;
      hLeg{(ii-1)*nS+jj} = ['\(d = ' int2str(d_vec(ii)) ', \ s_j = j^{-' int2str(s_power_vec(jj)) '}\)'];
   end
end
xlim([min(eps_vec)*0.8 max(eps_vec)*1.2])
ylim([1 10000])
[~,leg_icons] = legend(hds,hLeg, ...
   'box','off','location','eastoutside','orientation','vertical');
xlabel({'\(\varepsilon/||\hat{f}||_{\infty,}\)\boldmath\({}_{\gamma}\)'})
ylabel({'Sample size \(n\)'})
set(gcf,'Position',[200,200,1000,500]) %make figure big enough and the right aspect ratio
title('\boldmath\(w \, \)\unboldmath\( = (1,1,0.5, 0.25, \ldots ) \)')
for ii = 1:nD*nS
   leg_icons(ii).FontSize = 30;
end
print -depsc SampleSizeForDifferGamma.eps

