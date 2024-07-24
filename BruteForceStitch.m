function BruteForceStitch(expt, track_index_array)

% track_index_array, an array containing tracks to be stitched into one track in sequence,
% e.g. [1, 2, 3], which means track1, track2,  and track3 exist in order
% and not overlapped.
% After brute-force stitching, only the stitched long track is left.
% This function stitches the tracks without judging if they are right to
% be stitched. You need to judge.

for j = length(track_index_array):-1:2
    expt.track(track_index_array(j-1)).merge(expt.track(track_index_array(j)));
end

delete(expt.track(track_index_array(2:end)));
expt.track = expt.track(track_index_array(1));