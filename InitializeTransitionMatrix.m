function a = InitializeTransitionMatrix_v4(K,A)

% initialize transmission matrix
if K ~= 1
a = diag(ones(1,K)*A) + (ones(K,K)-eye(K))*(1-A)/(K-1);
a = a./(sum(a,2)*ones(1,K));
end 

if K == 1 %FH
    a = ones(1,1); %FH
end %FH

