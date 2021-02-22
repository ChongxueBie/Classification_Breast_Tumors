% set initial values of CEST and MT pools (required for numeric solution)

% do the following to initialize, or call 
close all
clear all
addpath(genpath(pwd));
Sim=init_Sim(struct());


%%

Sim.n_cest_pool=4;

Sim.dwA=0.0;
Sim.R1A=1/1.45;
Sim.R2A=1/0.033;
  
% Cr pool B (Cr)
Sim.n_cest_pool=1;
Sim.fB=0.0017;  % rel. conc 100mM/111M
Sim.kBA=150;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwB=2.;     % ppm  relative to dwA
Sim.R1B=1;      % R1B relaxation rate [Hz]
Sim.R2B=1/1e-3;     % R2B relaxation rate [Hz]

%MT
Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.MT_lineshape      = 'Lorentzian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.R1C=1;
Sim.fC=0.043;
Sim.dwC=0;
Sim.R2C=1/5e-5;  % 1/9.1µs
Sim.kCA=15; 
Sim.kAC=Sim.kCA*Sim.fC;

% amide pool D
Sim.n_cest_pool=2;
Sim.fD=0.0021;  % rel. conc 190mM/111M
Sim.kDA=50;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwD=3.5;     % ppm  relative to dwA
Sim.R1D=1;      % R1B relaxation rate [Hz]
Sim.R2D=1/1e-3;     % R2B relaxation rate [Hz]

% NOE pool E
Sim.n_cest_pool=4;
Sim.fE=0.0052;  % rel. conc 35mM/111M
Sim.kEA=45;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwE=-3.5;     % ppm  relative to dwA
Sim.R1E=1;      % R1B relaxation rate [Hz]
Sim.R2E=1/5.6e-4;     % R2B relaxation rate [Hz]

% 5 CEST pool F (-2.3 MT)
Sim.n_cest_pool=5;
Sim.fF=0.01;  % rel. conc 35mM/111M
Sim.kFA=15;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwF=-2.3;     % ppm  relative to dwA
Sim.R1F=1;      % R1B relaxation rate [Hz]Sim.R2G=1/7e-5;     % R2B relaxation rate [Hz]
Sim.R2F=1/7e-5;     % R2B relaxation rate [Hz]

Sim.FREQ=500;
Sim.Zi=1.;
Sim.B1=0.27;		 % the saturation B1 in µT
Sim.tp=4;		 % [s]
Sim.n=1;	
Sim.DC=1;
Sim.shape='block';		 
Sim.pulsed=0;		

Sim.xZspec = linspace(-6, 6, 81);%-5:0.1:5;

Sim.n_cest_pool=4;
Sim.analytic          = 1;                % calculate analtical solution? 1=yes, 0=no
Sim.numeric           = 0;                % calculate numerical solution? 1=yes, 0=no
Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.Rex_sol           = 'Hyper';          % solution for Rex - cases: 'Hyper', 
Sim.MT_lineshape      = 'Lorentzian';       % MT lineshape -SuperLorentzian, Gaussian, Lorentzian
Sim.MT_sol_type       = 'Rex_MT';         % Rex_MT solution type - cases: 'Rex_MT'
Sim.B1cwpe_quad   = -1;                     %XX
Sim.c     = 1;                              %XX
Sim.dummies=1; 
Sim.flipangle=5; 
% Sim.readout='FID'; 
Sim.spoilf=0; 
 %% plot multi-power
B1 = [0.5, 1.0, 2.0];
% B1 = [2.0];
simData = zeros(3, length(Sim.xZspec));
for i = 1: length(B1)
    Sim.B1 = B1(i);
    ana = ANALYTIC_SIM(Sim);
    simData(i,:) = ana.zspec;
    figure(1), plot(ana.x,ana.zspec,'k.-'); hold on
end
% ana = ANALYTIC_SIM(Sim);
% figure(1), plot(ana.x,ana.zspec,'k.-'); hold on

%%
path = "/Users/cbie1/Documents/CEST_Machine_learning_BreastTumor/Results/double_tumor_mouse/20200824/M3_2/M3_2_Ztab_MDA_MB_231.mat";
expdata = load(path);
expdata = (expdata.Ztab_MDA_MB_231)';
expfreq = linspace(-6, 6, 81);
for i = 1: 3
    figure(1),plot(expfreq, expdata(:,i),'ko')
end
set(gca,'XDir','reverse')
xlim([-6, 6])

ylim([0,1])
title('MDA-MB-231')
xlabel('Saturation Frequency (ppm)') 
ylabel('S/S_0 (%)')
% legend('0.5\muT')

% path = "/Users/cbie1/Documents/CEST_Machine_learning_BreastTumor/Results/double_tumor_mouse/20200824/M2/M2_Ztab_MCF_7.mat";
% expdata = load(path);
% expdata = (expdata.Ztab_MCF_7)';
% expfreq = linspace(-6, 6, 81);
% for i = 1: 3
%     figure(1),plot(expfreq, expdata(:,i),'ko')
% end
% set(gca,'XDir','reverse')
% xlim([-6, 6])

ylim([0,1])
title('MDA-MB-231')
xlabel('Saturation Frequency (ppm)') 
ylabel('S/S_0 (%)')
legend('0.5\muT', '1.0\muT', '2.0\muT')


%%
% savedir = '/Users/cbie1/Documents/CEST_Machine_learning_BreastTumor/Results/double_tumor_mouse/simulation_matlab';
% saveSimData = [savedir, '/sim_MDA_MB_231.mat'];
% save(saveSimData,'simData');