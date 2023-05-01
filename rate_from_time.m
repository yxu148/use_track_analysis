% Transferring from the time series in period of some incidents happen, t, to the
% rate of it happening r.
%[0, T] is the possible period the incident happens
% s stepsize, b binsize

function r = rate_from_time(t, T, s, b)

m = fix(T/s);
r = zeros(1, m + 1);

for j = 0 : m
    tleft = mod(j*s-b, T);  % periodic boundary condition
    tright = mod(j*s, T);
    if tleft>tright
        r(j+1) = nnz(t >= tleft | t < tright ) / b;
    else
        r(j+1) = nnz(t >= tleft & t < tright ) / b;
    end
end

end