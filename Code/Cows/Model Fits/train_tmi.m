function [Z,thetB,thetA,thetF,C] = train_tmi(xs_b,xs_a,xs_f,metz,nTot,alph,ctf_grad)

[nSamp nOb] = size(xs_b);
nOa = size(xs_a,2);nOf = size(xs_f,2);nMet = size(metz,2);

Z     = 0.1*randn(nSamp,nTot);
C     = 0.1*randn(nTot,nMet);
thetA = 0.1*randn(nTot,nOa);
thetB = 0.1*randn(nTot,nOb);
thetF = 0.1*randn(nTot,nOf);


%
% % hyperparameters
etaZ = 0.001;etaT = 0.001;etaC = 0.001;
QA = exp(-Z*thetA);QA = normalize(QA,2,'norm',1);
QB = exp(-Z*thetB);QB = normalize(QB,2,'norm',1);
QF = exp(-Z*thetF);QF = normalize(QF,2,'norm',1);


grdnorm = 1;
iter = 1;
while grdnorm > ctf_grad
    QA = exp(-Z*thetA);QA = normalize(QA,2,'norm',1);
    QB = exp(-Z*thetB);QB = normalize(QB,2,'norm',1);
    QF = exp(-Z*thetF);QF = normalize(QF,2,'norm',1);
    deltA = xs_a-QA;
    deltB = xs_b-QB;
    deltF = xs_f-QF;
    
    % gradients
    grthetA = (1-alph)*Z'*deltA;
    grthetB = (1-alph)*Z'*deltB;
    grthetF = (1-alph)*Z'*deltF;
    
    grzQA   = (1-alph)*deltA*thetA';
    grzQB   = (1-alph)*deltB*thetB';
    grzQF   = (1-alph)*deltF*thetF';
    grzM    = -2*alph*(metz-Z*C)*C';
    grz     = grzQA + grzQB + grzQF + grzM;
    
    grc     = -2*alph*Z'*(metz - Z*C);
    
    % update the variables
    Z      = Z - etaZ*grz;
    thetA  = thetA - etaT*grthetA;
    thetB  = thetB - etaT*grthetB;
    thetF  = thetF - etaT*grthetF;
    
    C      = C - etaC*grc;
    
    % errors    
    grdnorm = norm(grz)/norm(Z) + norm(grc)/norm(C) + norm(grthetA)/norm(thetA) + norm(grthetB)/norm(thetB) + norm(grthetF)/norm(thetF);
    
    % Output
    if mod(iter,500) == 0
        grdnorm
    end
    iter = iter + 1;
    
end

end

