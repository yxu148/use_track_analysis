% funtion tonoff_from_squarewave(sw, period)
%     sw, squarewave is a 1-d array, close to square wave when
%     plot(squarewave), with completely flat steps but not ideally vertical
%     edges.
%     period is a double, the time you want to map the one period of square
%     wave into time [0, period].
function [ton, toff] = tontoff_from_squarewave(sw, period)


slop_sw = zeros(1, length(sw));
slop_sw(1 : end-1) = diff(sw);

decreasing_sw = slop_sw < 0;  % 1 for decreasing, 0 for anything else, length(sw)
start_decreasing_sw  = [0, diff(decreasing_sw)] > 0.5;  % 1 for start decreasing, 0 for anything else, length(sw)
end_decreasing_sw = [0, diff(decreasing_sw)] < -0.5;  % 1 for start decreasing, 0 for anything else, length(sw)
decreasing_index_middle = round(0.5 * (find(end_decreasing_sw) + find(start_decreasing_sw)));  % frame number

increasing_sw = slop_sw > 0;  % 1 for decreasing, 0 for anything else
start_increasing_sw  = [0, diff(increasing_sw)] > 0.5;  % 1 for start decreasing, 0 for anything else
end_increasing_sw = [0, diff(increasing_sw)] < -0.5;  % 1 for start decreasing, 0 for anything else
increasing_index_middle = round(0.5 * (find(end_increasing_sw) + find(start_increasing_sw)));


% Discard the decreasing index if there is a very close increasing index
% nearby (within 0.5 s), which will be 10 frames apart from each other given 20 Hz
to_discard = [];  % indexes of elements to discard in decreasing_index_middle
for i = 1 : length(increasing_index_middle)
    distance = abs(decreasing_index_middle - increasing_index_middle(i));
    to_discard = [to_discard, find(distance < 10)];  %  10 frames
end
decreasing_index_middle(to_discard) = -10;
decreasing_index_middle_truly = decreasing_index_middle(decreasing_index_middle > 0);


% Create the zig-zag curve toff
toff = -100 * ones(1, length(sw));
toff(decreasing_index_middle_truly(1)) = period;
for i = 1 : length(decreasing_index_middle_truly) - 1
    left_index = decreasing_index_middle_truly(i) + 1;
    right_index = decreasing_index_middle_truly(i + 1);
    num_frames = right_index - left_index;
    toff(left_index : right_index) = 0 : period/num_frames : period;
end


% Create the zig-zag curve ton
ton = -100 * ones(1, length(sw));
ton(increasing_index_middle(1)) = period;
for i = 1 : length(increasing_index_middle) - 1
    left_index = increasing_index_middle(i) + 1;
    right_index = increasing_index_middle(i + 1);
    num_frames = right_index - left_index;
    ton(left_index : right_index) = 0 : period/num_frames : period;
end


% figure; 
% plot(sw, 'r'); hold on; 
% plot(decreasing_index_middle, 70, 'g.');
% plot(increasing_index_middle, 68, 'b.');
% plot(decreasing_index_middle(to_discard), 70, 'ro');
% plot(toff + 71, 'k');
% plot(ton + 50, 'm');
% hold off;