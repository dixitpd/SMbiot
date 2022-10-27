clc
clear
%

ctf_otu  = 0.001;% based on Brian's paper
alph     = 0.1;
nTot     = 20;

%% How good is the model fit
filen = strcat('species1_',num2str(nTot),'_',num2str(alph),'_',num2str(ctf_otu),'.mat');
load(filen)
QA = exp(-Z*thetA);QA = normalize(QA,2,'norm',1);
QB = exp(-Z*thetB);QB = normalize(QB,2,'norm',1);
QF = exp(-Z*thetF);QF = normalize(QF,2,'norm',1);
M  = Z*C;
%

k = 1;
% Species correlations
cas0 = diag(corr(QA,xs_a,'type','spearman'));
cbs0 = diag(corr(QB,xs_b,'type','spearman'));
cfs0 = diag(corr(QF,xs_f,'type','spearman'));
cms0 = diag(corr(metz,M,'type','spearman'));
cac0 = diag(corr(QA',xs_a','type','spearman'));
cbc0 = diag(corr(QB',xs_b','type','spearman'));
cfc0 = diag(corr(QF',xs_f','type','spearman'));
cmc0 = diag(corr(metz',M','type','spearman'));
%



load  ../three_kingdoms_cleaned_up

load ../TrainTest/archea_pred_20
load ../TrainTest/bacteria_pred_20
load ../TrainTest/fungi_pred_20
load ../TrainTest/metadata_pred_20

% %% Predicting community composition with respect to null model
% figure

subplot(3,1,1)
hold on
b = 0.05;
cacrand = corr(xs_a','type','spearman');cacrand(cacrand==1) = [];
[f,xi] = ksdensity(cac0(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'k','linewidth',2)
[f,xi] = ksdensity(cac(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'k--','linewidth',2)
[f,xi] = ksdensity(cacrand(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'k.','linewidth',2)
% 
subplot(3,1,2)
hold on
cbcrand = corr(xs_b','type','spearman');cbcrand(cbcrand==1) = [];
[f,xi] = ksdensity(cbc0(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'r','linewidth',2)
[f,xi] = ksdensity(cbc(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'r--','linewidth',2)
[f,xi] = ksdensity(cbcrand(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'r.','linewidth',2)

subplot(3,1,3)
hold on
cfcrand = corr(xs_f','type','spearman');cfcrand(cfcrand==1) = [];
[f,xi] = ksdensity(cfc0(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'b','linewidth',2)
[f,xi] = ksdensity(cfc(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'b--','linewidth',2)
[f,xi] = ksdensity(cfcrand(:,k),'Support',[-1 1],'BoundaryCorrection','reflection','Bandwidth',b);
plot(xi,f+4*(k-1),'b.','linewidth',2)
