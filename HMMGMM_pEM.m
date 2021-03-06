function hmmmodel = HMMGMM_pEM(deltaX, trackInfo, params)

hmmmodel = InitializeHMM(deltaX, params.pik, params.vacfk, params.transProb);

maxiter = params.maxiter;
tolerance = params.tolerance;
verbose = params.verbose;
splitIndex = trackInfo.splitIndex;
numFeatures = trackInfo.numFeatures;

p = hmmmodel.p;
a = hmmmodel.a;
b = hmmmodel.b;
sigma = hmmmodel.sigma;

K = size(a,1);

pbest = p;
abest = a;
bbest = b;
sigmabest = sigma;
gammankbest = 0;

% HMM model parameters
Lold = -Inf;
for i = 1:maxiter
    
    % expectation step
    [est,Lnew,gammank] = ExpectationHMM(p,a,b,splitIndex);
    disp([num2str(i) ': ' num2str(Lnew)]);
    
    % maximization step
    [p,a,b,sigma] = MaximizationHMM(deltaX,gammank,est,K,trackInfo);
    
    % check percentage similar
    [MAX,index] = max(gammank,[],2);

    if abs((Lnew - Lold)/(Lold)) < tolerance
        if verbose == 1
            disp('Converged...');
            disp(['iteration ' num2str(i) ': logL = ' num2str(Lnew)]);
        end
        break;
    elseif (Lold - Lnew) > tolerance
        if verbose == 1
            disp('Not Stable...');
            disp(['iteration ' num2str(i) ': logL = ' num2str(Lnew) 'oldlogL = ' num2str(Lold)]);
            
        end
        break;
    else
        Lold = Lnew;
        pbest = p;
        abest = a;
        bbest = b;
        sigmabest = sigma;
        gammankbest = gammank;
    end
    
end

if i == maxiter
    if verbose == 1
        disp('Did not converge...');
    end
end

vacfbest = zeros(K,numFeatures);
for i = 1:K
vacfbest(i,:) = sigmabest(1,1:numFeatures,i);
end
hmmmodel.vacf = vacfbest;
hmmmodel.p = pbest;
hmmmodel.a = abest;
hmmmodel.b = bbest;
hmmmodel.sigma = sigmabest;
hmmmodel.gammank = gammankbest;
hmmmodel.logL = Lold;

