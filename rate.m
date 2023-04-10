% Transferring from the time series of some incidents happen, t, to the
% rate of it happening r.
%[0, T] is the possible period the incident happens

function r = rate_from_time(t, T, s, b)

m = fix(T/s);
r = NaN(1, m + 1);

for j = 0 : m
    r[j+1] = nnz(t >= mod(j*s-b/2, T) & t < mod(j*s +b/2, T)) / b;
end

end