function [dBn]=dB_n_adjust(m, s, n, alpha)

% Function to adjust the site scatter of a paleointensity estimate for
% sample size following [1]. Please cite [1] if used.
%
% Note: This requires the Statistics Toolbox
%
% Written by Greig A. Paterson
%
% Input:
%       m - vector of mean values
%       s - vector of standard deviations
%       n - vector of the number of specimens
%       alpha - confidence level for adjustment
%
% Output:
%       dBn - the adjusted site scatter as a % of the mean
%
%
% References:
%       [1] Paterson, G. A., D. Heslop, and A. R. Muxworthy (2010), 
%           Deriving confidence in paleointensity estimates, 
%           Geochem. Geophys. Geosyst., 11, Q07Z18, doi:10.1029/2010GC003071.
%

%% Input checking
if nargin < 3
    error('dB_n_adjust:Input', 'At least 3 input arguments are required.');
end

if ~isequal(size(m), size(s), size(n))
    error('dB_n_adjust:Input', 'm, s, and n must the same size.');
end

if nargin < 4
    alpha=0.05;
end

%% Main function

tnc = nctinv(alpha, (n-1), m.*sqrt(n)./s);

dBn = 100.*abs(sqrt(n)./tnc);

end