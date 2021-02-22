% generation simulation data with different trainingIndex of metabolites
% and different powers

clear all
close all
clc
addpath(genpath(pwd));

%% define trainingIndex range

% relaxation time for water (coretx)
R1w = 1/1.28;
R2w = 1/0.013;
cestFreq = linspace(-6, 6, 81);
tp = 4; %duration time
B1 = 2.0;

numbers = 1000;

%% define trainingIndex
trainingIndex = zeros(1, numbers)+2;

%% generate simulation data n files
simData = zeros(length(B1), length(cestFreq), length(trainingIndex));
concentration = zeros(1,6);
for mm = 1:length(trainingIndex)
     concentration(1) = 150 + (160-150).*rand([1 1]);
%      normrnd(155.4,5); %Cr
     concentration(2) = 7600 + (11600-7600).*rand([1 1]);
%      normrnd(11100,500); %MT
     concentration(3) = 195 + (205-195).*rand([1 1]);
%      normrnd(200,5); %amide
     concentration(4) = 85 + (95-85).*rand([1 1]);
%      normrnd(89,5);%PCr
     concentration(5) = 95 + (195-95).*rand([1 1]);
%      normrnd(144,50);  %NOE
     concentration(6) = 5400 + (6400-5400).*rand([1 1]);
%      normrnd(5883,500);  %mt
     simData(:,:,mm) = Simulation_muscle(R1w, R2w, concentration/111000, B1, tp);

end
% sd_noise = [0.001, 0.003, 0.005, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03];
sd_noise = [0.005, 0.008, 0.01, 0.012, 0.015, 0.017, 0.02, 0.025, 0.03];
[x_data,y_data] = size(squeeze(simData));
data = simData;

sim_data_noise = zeros(length(sd_noise),81,1000);
for ij = 1:length(sd_noise)
    for ii = 1:y_data
        noise = rician_noise(length(cestFreq), sd_noise(ij));
        sim_data_noise(ij, :, ii) = squeeze(data(:,:,ii)) + noise;
    end    
end
fprintf('-------finish!-------\n')

%% save data
savedir = '/Volumes/CX/JHU/Machine_Learning_BreastTumor/results/20210104/simulation_matlab';
saveSimData = [savedir, '/training_simDataMuscle_noise_2.0uT_20210201.mat'];
save(saveSimData,'sim_data_noise');
saveIndex = [savedir, '/training_IndexMuscle_noise_2.0uT_20210201.mat'];
save(saveIndex,'trainingIndex');

%% compare with in vivo data

figure(1)
plot(cestFreq, sim_data_noise(6,:,7),'k-'); hold on
set(gca,'XDir','reverse')

path = "/Volumes/CX/JHU/Machine_Learning_BreastTumor/results/20210104/old_mice/M3-2/M3_2_Ztab_MCF_7.mat";
expdata = load(path);
expdata = (expdata.Ztab_MCF_7)';
expfreq = linspace(-6, 6, 81);
for i = 3: 3
    figure(1),plot(expfreq, expdata(:,i),'ro')
end
set(gca,'XDir','reverse')
xlim([-6, 6])

ylim([0,1])
title('MDA-MB-231')
xlabel('Saturation Frequency (ppm)') 
ylabel('S/S_0 (%)')


path = "/Volumes/CX/JHU/Machine_Learning_BreastTumor/results/20210104/old_mice/M3-2/M3_2_Ztab_Muscle_back.mat";
expdata = load(path);
expdata = (expdata.Ztab_Muscle)';
expfreq = linspace(-6, 6, 81);
for i = 3:3
    figure(1),plot(expfreq, expdata(:,i),'ko')
end
set(gca,'XDir','reverse')
xlim([-6, 6])





