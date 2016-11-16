% alignment_quality_IK calculates the spectral data alignment quality
% measure based on the correlations between the samples
%
% INPUT
% sp            : spectral data matrix
% ppm           : x-axis label corresponding to NMR chemical shift in ppm
% bin_size      : size of each bin in ppm for initial auto-scaling step of the bins
%
% OUTPUT
% aq_bin        : alignment quality with respect to the predefined bin size
% cc_bin_vec    : vector of correlation coefficients for each sample pair
%
% Author, I. Karaman, Imperial College London, 2016
%
function [aq_bin,cc_bin_vec] = alignment_quality_IK(sp,ppm,bin_size)

[m,n] = size(sp);
min_ppm = min(ppm);
no_bins = 1;
bin_gr = zeros(n,1);

for i = 1:n
    if ppm(i) >= min_ppm && ppm(i) < min_ppm + bin_size
        bin_gr(i,1) = no_bins;
    else
        bin_gr(i,1) = no_bins+1;
        no_bins = no_bins+1;
        min_ppm = ppm(i);
    end
end

sp_uv = [];
for k = 1:no_bins
    sp_mean = mean(sp(:,bin_gr == k),2);
    sp_std = std(sp(:,bin_gr == k),0,2); sp_std(sp_std == 0) = 1;
    sp_uv = [sp_uv (sp(:,bin_gr == k) - (sp_mean*ones(1,sum(bin_gr == k))))./(sp_std*ones(1,sum(bin_gr == k)))];
end

cc_bin_mat = corr(sp_uv');

cc_bin_vec = [];
for i = 1:m
    cc_bin_vec = [cc_bin_vec; cc_bin_mat(i+1:end,i)];
end
    
aq_bin = mean(cc_bin_vec);

end




