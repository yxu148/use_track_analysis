function v_prop = proplize(v)
% Proportion-lize a non-negative vector so that the sum of v is 1
% input v should be a 1-dimension vector, row or column

if sum(v) ~= 0
    v_prop = v / sum(v);
else
    v_prop = v;
end
