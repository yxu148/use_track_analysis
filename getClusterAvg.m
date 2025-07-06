function avgSpeeds = getClusterAvg(speed, logicArray)
%GETCLUSTERAVG Calculates average quantity (eg. speed) for each cluster of 1s indicated in logicArray.
%
%   avgSpeeds = GETCLUSTERAVG(speed, logicArray)
%
%   Inputs:
%     - speed:      Numeric array of speeds (e.g., 1x1000).
%     - logicArray: Logical array of the same size, with clusters of 1s.
%
%   Output:
%     - avgSpeeds:  Array containing average speed for each cluster of 1s.

    % Ensure inputs are the same size
    if length(speed) ~= length(logicArray)
        error('Inputs "speed" and "logicArray" must be the same length.');
    end

    % Find start and end indices of clusters of 1s
    diffLogic = diff([0, logicArray, 0]);
    startIdx = find(diffLogic == 1);
    endIdx   = find(diffLogic == -1) - 1;

    % Number of clusters
    nClusters = length(startIdx);
    avgSpeeds = zeros(1, nClusters);

    % Compute average speed for each cluster
    for i = 1:nClusters
        avgSpeeds(i) = mean(speed(startIdx(i):endIdx(i)));
    end
end
