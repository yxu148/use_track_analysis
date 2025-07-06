function cmap = makeColormap(colors, weights, n)
% makePiecewiseColormapWeighted generates a custom colormap with adjustable segment lengths.
%
% Inputs:
%   - colors: Mx3 matrix of RGB keypoints (e.g., [0 0 1; 0.5 0.5 0.5; 1 0 0])
%   - weights: (M-1)x1 vector of relative lengths of each color segment
%   - n: total number of entries in the colormap (default 256, i.e. continous)
%
% Output:
%   - cmap: nx3 colormap matrix

    if nargin < 3
        n = 256;
    end

    if size(colors, 1) - 1 ~= length(weights)
        error('Number of weights must be one less than number of colors.');
    end

    % Normalize weights to sum to 1
    weights = weights / sum(weights);

    % Compute segment lengths in colormap
    segmentLengths = round(weights * n);

    % Adjust rounding errors to ensure total length = n
    diff = n - sum(segmentLengths);
    if diff ~= 0
        segmentLengths(end) = segmentLengths(end) + diff;
    end

    cmap = [];

    for i = 1:length(segmentLengths)
        ns = segmentLengths(i);
        r = linspace(colors(i,1), colors(i+1,1), ns)';
        g = linspace(colors(i,2), colors(i+1,2), ns)';
        b = linspace(colors(i,3), colors(i+1,3), ns)';
        cmap = [cmap; [r g b]];
    end
end
