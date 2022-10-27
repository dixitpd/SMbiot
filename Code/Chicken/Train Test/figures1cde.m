clc
clear
load ../OTUTables
load testing_data
nTest = length(testing_data);
swab_test  = swabOTU(testing_data,:);
cecum_test = cecumOTU(testing_data,:);
ileum_test = ileumOTU(testing_data,:);


nC = 20;
filen = strcat('combined_',num2str(nC),'.mat')
load(filen)

eta = 0.005;
Ztest  = 0.1*randn(nTest,nC);
%
flg   = 1;iter = 1;

while flg > 0
    % predictions
    QS    = exp(-[Ztest]*thetS);QS = normalize(QS,2,'norm',1);
    deltS = swab_test-QS;
    
    % Gradients
    grZcom    = deltS*thetS';
    nrg       = norm(grZcom)/norm(Ztest);
    
    if nrg < 0.005
        flg = 0;
    end
    
    % update
    Ztest     = Ztest  - eta*grZcom;
    if mod(iter,500) == 0
        nrg
    end
    iter = iter + 1;
end

QI    = exp(-[Ztest]*thetI);QI = normalize(QI,2,'norm',1);
QC    = exp(-[Ztest]*thetC);QC = normalize(QC,2,'norm',1);

cc = 1;
for i=1:nTest-1
    for j=i+1:nTest
        cic_0(cc) = corr(ileum_test(i,:)',ileum_test(j,:)','type','spearman');
        ccc_0(cc) = corr(cecum_test(i,:)',cecum_test(j,:)','type','spearman');
        cc = cc + 1;
    end
end
cic_t = (diag(corr(QI',ileum_test','type','spearman')));
ccc_t = (diag(corr(QC',cecum_test','type','spearman')));

cis_t = (diag(corr(QI,ileum_test,'type','spearman')));
ccs_t = (diag(corr(QC,cecum_test,'type','spearman')));

 [mean(cic_t) mean(cic_0)]
 [mean(ccc_t) mean(ccc_0)]


%
load ../OTUTables

[nSamp nOs] = size(swabOTU);
[nSamp nOi] = size(ileumOTU);
[nSamp nOc] = size(cecumOTU);
% %

nC = 20;
filen = strcat('../Model Fits/combined_',num2str(nC),'.mat')
load(filen)
QSt    = exp(-[Zcom]*thetS);QSt = normalize(QSt,2,'norm',1);
QIt    = exp(-[Zcom]*thetI);QIt = normalize(QIt,2,'norm',1);
QCt    = exp(-[Zcom]*thetC);QCt = normalize(QCt,2,'norm',1);

%
css = diag(corr(QSt,swabOTU,'type','spearman'));
cis = diag(corr(QIt,ileumOTU,'type','spearman'));
ccs = diag(corr(QCt,cecumOTU,'type','spearman'));
csc = diag(corr(QSt',swabOTU','type','spearman'));
cic = diag(corr(QIt',ileumOTU','type','spearman'));
ccc = diag(corr(QCt',cecumOTU','type','spearman'));

OTUs = unique([goodOTU_swab goodOTU_ileum goodOTU_cecum]);

b = 0.05;
subplot(3,1,1)
hold on
[f,xi] = ksdensity(cic,'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f,'b','linewidth',2)
[f,xi] = ksdensity(cic_t,'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f,'r','linewidth',2)
[f,xi] = ksdensity(cic_0,'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f,'k','linewidth',2)

subplot(3,1,2)
hold on
[f,xi] = ksdensity(ccc,'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f,'b','linewidth',2)
[f,xi] = ksdensity(ccc_t,'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f,'r','linewidth',2)
[f,xi] = ksdensity(ccc_0,'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f,'k','linewidth',2)

subplot(3,1,3)
hold on
plot(cis,cis_t,'o')
plot(ccs,ccs_t,'o')
