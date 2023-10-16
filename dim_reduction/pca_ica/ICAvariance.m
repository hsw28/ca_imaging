function var_accounted_for = ICAvariance(Zica, W, T, mu, OG_data)

% Assuming Zica, W, T, and mu are already defined
[nComponents, nSamples] = size(Zica);

% Original data reconstruction
Z_reconstructed = T \ (W' * Zica) + repmat(mu,1,nSamples);

% Compute the total variance of the original data
total_variance = sum(var(OG_data,0,2));

%Compute the total variance of the reconstructed data:
reconstructed_variance = sum(var(Z_reconstructed,0,2));

%compute residual variances
residual_variance = total_variance - reconstructed_variance

total_variance
residual_variance

var_accounted_for = reconstructed_variance./total_variance;

%To determine the variance captured by each IC, one can project the data onto each IC and compute the variance.
%However, keep in mind that in ICA, components are typically scaled to have unit variance,
%so this might not give the intuitive result of explaining variance like in PCA.
%Still, if desired, it can be done as:
ic_variances = var(Zica,0,2);

%With the above, reconstructed_variance provides the variance captured by all the ICs combined.
%The residual_variance gives you the variance that's not captured by the ICs.
%If you need a percentage representation, you can divide the reconstructed_variance by total_variance.


%{
% Initialize the variance accounted for by each IC
variancesIC = zeros(1, nComponents);

for i = 1:nComponents
    % Set all ICs to zero except the i-th
    tempZica = zeros(size(Zica));
    tempZica(i, :) = Zica(i, :);

    % Reconstruct the data using only the i-th IC
    tempReconstruction = T \ W' * tempZica + repmat(mu, 1, nSamples);

    % Compute the variance of the reconstruction
    variancesIC(i) = sum(var(tempReconstruction, 0, 2));
end

% Normalize the variances by the total variance to get the proportion
proportionVariances = variancesIC / totalVariance;
%}
