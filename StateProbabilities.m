function [gamma,xi,logL] = StateProbabilities(a,b,alpha,beta,scale)

[K,T] = size(b);

logalpha = log(alpha+eps);
logbeta = log(beta+eps);

loga = log(a+eps);
logb = log(b+eps);
logscale = log(scale+eps);

% state occupation probability to be in state j at time t 
% gamma = alpha.*beta;
gamma = exp(logalpha + logbeta);

%normalize gamma

% state transition probability to be in state i at time t and in state j at time t+1 
xi = zeros(K,K,T-1);
for t = 1:T-1
    for i = 1:K
        for j = 1:K
%             xi(i,j,t) = alpha(i,t)*a(i,j)*b(j,t+1)*beta(j,t+1)/scale(t+1);
            xi(i,j,t) = logalpha(i,t) + loga(i,j) + logb(j,t+1) + logbeta(j,t+1) - logscale(t+1);
        end
    end
end
xi = exp(xi);

% log likelihood
logL = sum(logscale);





