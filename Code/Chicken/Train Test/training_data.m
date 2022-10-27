clc
clear
load ../OTUTables
load testing_data
swabOTU(testing_data,:) = [];
cecumOTU(testing_data,:) = [];
ileumOTU(testing_data,:) = [];


[nSamp nOs] = size(swabOTU);
[nSamp nOi] = size(ileumOTU);
[nSamp nOc] = size(cecumOTU);
% %

nC = 20;
filen = strcat('combined_',num2str(nC),'.mat')

% %
eta = 0.005;
thetS = 0.1*randn(nC,nOs);
thetI = 0.1*randn(nC,nOi);
thetC = 0.1*randn(nC,nOc);
Zcom  = 0.1*randn(nSamp,nC);
% 
flg   = 1;iter = 1;

while flg > 0    
    % predictions
    QS    = exp(-[Zcom]*thetS);QS = normalize(QS,2,'norm',1);
    deltS = swabOTU-QS;
    QI    = exp(-[Zcom]*thetI);QI = normalize(QI,2,'norm',1);
    deltI = ileumOTU-QI;
    QC    = exp(-[Zcom]*thetC);QC = normalize(QC,2,'norm',1);
    deltC = cecumOTU-QC;
    
    % Gradients
    grthetS     = [Zcom]'*deltS;
    grthetI     = [Zcom]'*deltI;
    grthetC     = [Zcom]'*deltC;
    grZS_tot    = deltS*thetS';
    grZI_tot    = deltI*thetI';
    grZC_tot    = deltC*thetC';
    
    grZcom      = grZS_tot + grZI_tot + grZC_tot;

    nrg = norm(grthetS)/norm(thetS) + norm(grthetI)/norm(thetI) + norm(grthetC)/norm(thetC) + norm(grZcom)/norm(Zcom);
    
    if nrg < 0.005
        flg = 0;
    end
    
    % update
    thetS    = thetS - eta*grthetS;
    thetI    = thetI - eta*grthetI;
    thetC    = thetC - eta*grthetC;
    Zcom     = Zcom  - eta*grZcom;
    %
    

    if mod(iter,500) == 0
        nrg
    end
    iter = iter + 1;
end
save(filen,'Zcom','thetI','thetS','thetC')
