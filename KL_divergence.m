% function kld = KL_divergence(P, Q)
% The inputs, P and Q, are 1-dimensional arrays with the same length
% Only keep the positive terms
% Return a double, D_KL(P||Q), called Kullback–Leibler divergence (also called relative entropy and I-divergence) 

function [kld, kld_terms] = KL_divergence(P, Q)

index = (P>0)&(Q>0);  % select positive values
kld = sum(P(index) .* log(P(index) ./ Q(index)));
kld_terms = zeros(1, length(P));
kld_terms(index) = P(index) .* log(P(index) ./ Q(index));  