function Re=SolarRecharge(K,F,N)
% Minimum and maximum values for the maximum energy consumption matrix.
% Variations are U ~ [B_min, B_max]:
Re_min = 0; Re_max = 0.01;

% Correlation coefficient between the energy consumptions in one frame to the next.
% 0 <= rho < 1. Higher correlation means lower variations due to mobility speed,
% unused activity level, etc.:
rho = 0.98;

%% Main Loop
% Correlation matrix of the actual maximum energy consumptions:
M_corr=zeros(K*F,K*F);
for i=1:K*F; for j=1:K*F; M_corr(i,j) = 1-(abs(i-j)*(1-rho)); end; end
% Correlation matrix for Gaussian-distributed consumptions:
M_corr = M_corr.*(M_corr>=0);
% Adjusting correlations for uniformly-distributed consumptions:
for i = 1:K*F
    for j = max(K*F-1,i):K*F
        if i ~= j
            M_corr(i, j) = 2 * sin(pi * M_corr(i, j) / 6);
            M_corr(j, i) = 2 * sin(pi * M_corr(j, i) / 6);
        end
    end
end
% Inducing correlation:
Ch = chol(M_corr); G = randn(N,K*F); Re = G * Ch;
% Actual maximum consumption matrix with imposed limits:
Re = normcdf(Re)*(Re_max-Re_min) + Re_min;