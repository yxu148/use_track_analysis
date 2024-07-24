% Transferring from the time series in period of some incidents happen, t, to the
% rate of it happening r, and the corresponding standard deviation
% rate = count / time; rate_std = sqrt(count) / time;
%[0, T] is the possible period the incident happens
% s stepsize, b binsize

function [r, r_std] = rate_from_time(t, T, s, b)

m = fix(T/s);
r = zeros(1, m + 1);
r_std = zeros(1, m + 1);

for j = 0 : m
    tleft = mod(j*s-b, T);  % use mod to have periodic boundary condition
    tright = mod(j*s, T);  % r(t) is the calculated between [t-b/2, t+b/2), or [t, t+b)
    if tleft>tright
        r(j+1) = nnz(t >= tleft | t < tright ) / b;
        r_std(j+1) = sqrt( nnz(t >= tleft | t < tright) )/ b;
    else
        r(j+1) = nnz(t >= tleft & t < tright ) / b;
        r_std(j+1) = sqrt( nnz(t >= tleft & t < tright) )/ b;
    end
end

end