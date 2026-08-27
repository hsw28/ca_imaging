function c = nanconv(a, kernel, varargin)
%NANCONV Convolve while ignoring NaN samples.
%   C = NANCONV(A, KERNEL) performs a two-dimensional convolution with
%   output size 'same'. NaN values in A do not contribute to the result.
%
%   Supported options (kept compatible with the calls in this project):
%     'same', 'full', 'valid' - convolution output shape
%     'edge'                  - correct for a truncated kernel at edges
%     'noedge'                - retain ordinary convolution edge scaling
%     'nanout'                - restore NaNs at NaN locations in A
%     'nonanout'              - return interpolated values at those locations
%
%   This local implementation removes the dependency on the third-party
%   File Exchange nanconv utility formerly required by the analysis code.

shape = 'same';
correctEdges = false;
restoreNaNs = false;

for i = 1:numel(varargin)
    option = lower(char(varargin{i}));
    switch option
        case {'same', 'full', 'valid'}
            shape = option;
        case 'edge'
            correctEdges = true;
        case 'noedge'
            correctEdges = false;
        case 'nanout'
            restoreNaNs = true;
        case 'nonanout'
            restoreNaNs = false;
        otherwise
            error('nanconv:UnknownOption', 'Unknown option ''%s''.', option);
    end
end

if ~isnumeric(a) || ~isnumeric(kernel) || ndims(a) > 2 || ndims(kernel) > 2
    error('nanconv:InvalidInput', 'A and KERNEL must be numeric 2-D arrays.');
end
if isempty(a) || isempty(kernel)
    c = conv2(a, kernel, shape);
    return
end

nanMask = isnan(a);
valid = ~nanMask;
a(nanMask) = 0;

% The numerator is the ordinary convolution with missing samples set to
% zero. availableWeight measures how much of the kernel had valid data.
numerator = conv2(a, kernel, shape);
availableWeight = conv2(double(valid), kernel, shape);

if correctEdges
    % Renormalize both missing samples and portions of the kernel outside A.
    referenceWeight = sum(kernel(:)) .* ones(size(numerator));
else
    % Renormalize missing samples only, preserving conv2's edge behavior.
    referenceWeight = conv2(ones(size(a)), kernel, shape);
end

c = numerator .* referenceWeight ./ availableWeight;
c(availableWeight == 0) = NaN;

if restoreNaNs
    if ~strcmp(shape, 'same')
        error('nanconv:NanOutShape', ...
            '''nanout'' is supported only with the ''same'' output shape.');
    end
    c(nanMask) = NaN;
end
end
