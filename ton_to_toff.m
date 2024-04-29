function toff = ton_to_toff(ton, duration_high, duration_low)
% function toff = ton_to_toff(ton, duration_high, duration_low)
%
% ton is 1-by-length(ton) array, with the time (s) as elements. The time is
% in period starting with intensity high (on).
% duration_high is a double, means the duration of high intensity in one
% period in seconds.
% duration_low is a double, means the duration of low intensity in one
% period in seconds.

toff = zeros(1, length(ton));
toff(ton > duration_high) = ton(ton > duration_high) - duration_high;
toff(ton <= duration_high) = ton(ton <= duration_high) + duration_low;