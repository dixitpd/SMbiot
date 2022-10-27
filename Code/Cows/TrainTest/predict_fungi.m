clc
clear
%

load ../../Cleanup/three_kingdoms_cleaned_up


%% Only abundant OTUs in species 1 of cows
sps1 = find(cow_species==1);

ctf_otu      = 0.001;% based on Brian's paper
bacteria_xs  = bacteria_xs(sps1,:);
archea_xs    = archea_xs(sps1,:);
fungi_xs     = fungi_xs(sps1,:);
metadata_cow = metadata_cow(sps1,:);
mn_bact      = mean(bacteria_xs);
goodb        = find(mn_bact > ctf_otu);
mn_arch      = mean(archea_xs);
gooda        = find(mn_arch > ctf_otu);
mn_fung      = mean(fungi_xs);
goodf        = find(mn_fung > ctf_otu);
%
restb = 1 - sum(bacteria_xs(:,goodb)')';
xs_b  = [bacteria_xs(:,goodb) restb];
resta = 1 - sum(archea_xs(:,gooda)')';
xs_a  = [archea_xs(:,gooda) resta];
restf = 1 - sum(fungi_xs(:,goodf)')';
xs_f  = [fungi_xs(:,goodf) restf];
%
mu_met = mean(metadata_cow);st_met = std(metadata_cow);
metz   = (metadata_cow-mu_met)./st_met;

%% Filter samples based on large deviations
badsamps         = find(sum(abs(metz > 10)'));
mcow             = metadata_cow;
mcow(badsamps,:) = [];
xs_a(badsamps,:) = [];
xs_f(badsamps,:) = [];
xs_b(badsamps,:) = [];
%% Testing data
load testing_data
nTest     = length(testing_data);
xs_a_test = xs_a(testing_data,:);
xs_b_test = xs_b(testing_data,:);
xs_f_test = xs_f(testing_data,:);
mcow_test = mcow(testing_data,:);

xs_a(testing_data,:) = [];
xs_b(testing_data,:) = [];
xs_f(testing_data,:) = [];
mcow(testing_data,:) = [];
mu_met           = mean(mcow);
st_met           = std(mcow);
metz             = (mcow-mu_met)./st_met;
metz_test        = (mcow_test-mu_met)./st_met;

%% Train models

nTot     = 20;
ctf_grad = 0.01;
alph     = 0.1;

filen = strcat('species1_',num2str(nTot),'_',num2str(alph),'_',num2str(ctf_otu),'.mat')
load(filen)

Ztst = randn(nTest,nTot);

etaZ = 0.001;etaT = 0.001;etaC = 0.001;
QB = exp(-Ztst*thetB);QB = normalize(QB,2,'norm',1);
QA = exp(-Ztst*thetA);QA = normalize(QA,2,'norm',1);


grdnorm = 1;
iter = 1;
while grdnorm > ctf_grad
    QB    = exp(-Ztst*thetB);QB = normalize(QB,2,'norm',1);
    QA    = exp(-Ztst*thetA);QA = normalize(QA,2,'norm',1);
    deltB = xs_b_test-QB;
    deltA = xs_a_test-QA;
    
    % gradients
    
    grzQB   = (1-alph)*deltB*thetB';
    grzQA   = (1-alph)*deltA*thetA';
    grzM    = -2*alph*(metz_test-Ztst*C)*C';
    grz     = grzQB + grzQA + grzM;
        
    % update the variables
    Ztst   = Ztst - etaZ*grz;    
    % errors    
    grdnorm = norm(grz)/norm(Z);
    
    % Output
    if mod(iter,500) == 0
        grdnorm
    end
    iter = iter + 1;
    
end

QF = exp(-Ztst*thetF);QF = normalize(QF,2,'norm',1);

cfc = (diag(corr(QF',xs_f_test','type','spearman')));
cfs = (diag(corr(QF,xs_f_test,'type','spearman')));
[mean(cfc) std(cfc)]
[mean(cfs) std(cfs)]
%save fungi_pred_20_0.05.mat cfc cfs