function vq = multi_interp1(x, v, xq, interp_dim, varargin)
% Interpolate along a specific dimension of an N-D array.
% 
%% Syntax
% vq = multi_interp1(x, v, xq, interp_dim)
% vq = multi_interp1(x, v, xq, interp_dim, varargin)
%
%% Description 
% vq = multi_interp1(x, v, xq, interp_dim) interpolates along a specific
%       dimension of an N-D array using interp1.
% vq = multi_interp1(x, v, xq, interp_dim, varargin) specifies the
%       interpolate options of interp1.
%
%% Examples
%
%
%% Input Arguments
%   x          - 1D coordinate array along `interp_dim` (length m)
%   v          - N-D array of values (size(interp_dim) must equal length(x))
%   xq         - Interpolation points (vector)
%   interp_dim - The dimension of `v` along which to interpolate
%   varargin   - Optional parameters passed to interp1 (method, 'extrap', etc.)
%
%% Output Arguments
%   vq         - Interpolated array with size(vq, interp_dim) = length(xq)
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 17 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
% Move interpolation dimension to the first dimension
perm_order = 1:ndims(v);
perm_order([1, interp_dim]) = perm_order([interp_dim, 1]);  % swap interp_dim <-> 1
v_perm = permute(v, perm_order);

% Reshape to 2D for interpolation
sz = size(v_perm);
v_2d = reshape(v_perm, sz(1), []);

% Interpolation (now every column is a curve along x)
vq_2d = interp1(x(:), v_2d, xq(:), varargin{:});  % [length(xq) x N]

% Reshape back to permuted shape
sz(1) = length(xq);
vq_perm = reshape(vq_2d, sz);

% Restore original dimension order
vq = ipermute(vq_perm, perm_order);

end