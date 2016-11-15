% pca_nipals_IK performs PCA using NIPALS algorithm
%
% INPUT
% A         : data matrix
% K         : number of principal components
% scaling   : scaling method
%             'mean' - mean-centering
%             'auto' - auto-scaling
%             'pareto' - Pareto scaling
%             0 - no scaling
%
% OUTPUT
% res       : variable for the results containing
%             T - scores
%             P - loadings
%             ExpVar - explained variance
%             PEV - percent explained variance
%             DModX - distance to the model values for each sample
%             Dcrit - critical distance (1st row: 95%, 2nd row: 99%)
%             HotellingT2 - Hotelling's T2 values for each sample
%             CorrLoad - correlation loadings
%             P_transformed - loadings scaled to original variable scale (useful if data matrix is auto-scaled)
%
% Author, I. Karaman, Imperial College London, 2016
%
function res = pca_nipals_IK (A,K,scaling)

%% Calculating the components

[m,n] = size(A);

switch scaling
    case 'mean'
        B = A - ones(m,1)*mean(A);
    case 'auto'
        B = (A - ones(m,1)*mean(A))./(ones(m,1)*std(A));
        B(isnan(B)) = 0;
    case 'pareto'
        B = (A - ones(m,1)*mean(A))./sqrt(ones(m,1)*std(A));
        B(isnan(B)) = 0;
    case 0
        B = A;
end

E = B;
T = zeros(m,K);
P = zeros(n,K);
ExpVar = zeros(K,1);

for i = 1:K
    
    [~,idx]=max(sum(E.*E));
    t_old=E(:,idx);
    it = 0;
    converged = false;
    while ~converged
        p = E'*t_old;
        p = p/norm(p);
        t_new = E*p;
        converged = (t_new - t_old)'*(t_new - t_old) < 1e-12;
        t_old = t_new;
        it = it + 1;
        if it == 5000
            disp(['PC' num2str(i) ' not converged!'])
            break
        end
    end
    
    if mod(i,10) == 0
        disp(['PC' num2str(i) ' with ' num2str(it) ' iterations' ])
    end
    
    T(:,i) = t_new;
    P(:,i) = p;
    
    E = E - t_new*p';
    
    ExpVar(i) = t_new'*t_new;
    
end

perExpVar = 100*ExpVar ./ sum(sum(B.^2));

res.T = T;
res.P = P;
res.ExpVar = ExpVar;
res.PEV = perExpVar;

%% DModX
for k = 1:K
    EE = B - T(:,1:k)*P(:,1:k)';
    s_i(:,k) = (m/(m-k-1))*sqrt(sum(EE.^2,2)/(n-k));
    s_0(1,k) = sqrt(sum(sum(EE.^2))/((m-k-1)*(n-k)));
    
    %     DFmod(1,k) = sqrt((m-k-1)*(n-k));
    %     if n >DFmod(1,k)
    %         DFobs(1,k) = (min([n,100,DFmod(1,k)])+sqrt(n-DFmod(1,k))-k)/(m/(m-k-1));
    %     else
    %         DFobs(1,k) = (min([n,100,DFmod(1,k)])-k)/(m/(m-k-1));
    %     end
    
    DFmod(1,k) = (m-k-1)*(n-k);
    DFobs(1,k) = (n-k);
    
    Dcrit(1,k) = sqrt(finv(0.95,DFobs(1,k),DFmod(1,k)));
    Dcrit(2,k) = sqrt(finv(0.99,DFobs(1,k),DFmod(1,k)));
    
end
res.DModX = s_i./(ones(m,1)*s_0);
res.Dcrit = Dcrit;

%% Hotelling T2

for k = 1:K
    S2(k,k) = (T(:,k)'*T(:,k))/(m-1);
    T2crit(1,k) = ((k*(m-1))/(m-k))*finv(0.95,k,m-k);
    T2crit(2,k) = ((k*(m-1))/(m-k))*finv(0.99,k,m-k);
end

for i = 1:m
    for k = 1:K
        res.HotellingT2(i,k) = T(i,1:k)*(inv(S2(1:k,1:k)))*T(i,1:k)';
    end
end
res.T2crit = T2crit;

%% Correlation Loadings and Transformed Loadings

stdevA = std(A);
stdevT = std(res.T);
for i = 1:n
    for j = 1:K
        temp = corrcoef(res.T(:,j),A(:,i));
        res.Corr_Load(i,j) = temp(2,1);
    end
end
res.P_transformed = res.Corr_Load.*((stdevA'*ones(1,K))./(ones(n,1)*stdevT));

end

