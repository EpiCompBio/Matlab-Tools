function [data_norm,factors,refsample]= normalize_samples(data_ini,method,ind_refsample)
% normalizes samples (rows) of a matrix using a normalization method
%
% INPUT:
%	data_ini        : data matrix
% 	method          : 'TA': total area normalization
%                     'PQ': probabilistic quotient normalization
%                     'QMEDIAN': quantile equating/normalization
%	ind_refsample   : a vector containing the index of the reference sample(s). 
%                     Default is the column-wise median of data_ini, i.e., the median spectrum
%
% OUTPUT:
%   data_norm       : normalized data matrix
%   factors         : row vector of normalization factors estimated for each sample
%   refsample       : row vector of reference sample

factors = nan(size(data_ini,1),1);
data_norm = nan(size(data_ini));

switch method
    case 'TA' % total area
        for i=1:size(data_ini,1)
            factors(i)=sum(data_ini(i,:));
            data_norm(i,:) = 10^6*data_ini(i,:)./factors(i);
        end
    case 'PQ' % probabilistic quotient
        for i=1:size(data_ini,1)
            factors(i)=sum(data_ini(i,:));
            data_norm(i,:) = 10^6*data_ini(i,:)./factors(i);
        end
        data_norm(data_norm==0)= 10^-15;
        if nargin < 3
            refsample = nanmedian(data_norm);
        else
            refsample = nanmedian(data_norm(ind_refsample,:),1);
        end
        quotients = data_norm./(refsample(ones(1,size(data_ini,1)),:));
        for i=1:size(data_ini,1)
            data_norm(i,:) = data_norm(i,:)./nanmedian(quotients(i,:));
            factors(i)=(factors(i)*nanmedian(quotients(i,:)))/10^6;
        end
    case 'QMEDIAN' % quantile equating
        data_ini = data_ini';
        [sorted_data_ini,ind1] = sort(data_ini);
        [~,ind2] = sort(ind1);
        median_sorted_data_ini = nanmedian(sorted_data_ini,2);
        sorted_data_ini_new = median_sorted_data_ini*ones(1,size(data_ini,2));
        data_norm = [];
        for j = 1:size(sorted_data_ini,2)
            data_norm(:,j) = sorted_data_ini_new(ind2(:,j),j);
        end
        data_norm = data_norm';
end


end
