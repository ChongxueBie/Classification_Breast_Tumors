% set initial values of CEST and MT pools (required for numeric solution)

% do the following to initialize, or call 
close all
clear all
Sim=init_Sim(struct());


%%

Sim.n_cest_pool=5;

Sim.dwA=0;
Sim.R1A=1/0.8;
Sim.R2A=1/0.03;

Sim.fB=0.006;
Sim.kBA=2000;
Sim.dwB=0.9;
Sim.R1B=1; 
Sim.R2B=1/0.55;

% CEST pool 'D'
Sim.fD=0.0005;     
Sim.kDA=150;            
Sim.dwD=2;          
Sim.R2D=1/9.7e-3;           
Sim.R1D=1;            

% CEST pool 'E'
Sim.fE=0.004;    
Sim.kEA=5500;           
Sim.dwE=2.8;          
Sim.R2E=1/4e-3;           
Sim.R1E=1;            

% CEST pool 'F'
Sim.fF=0.002;     
Sim.kFA=30;            
Sim.dwF=3.7;          
Sim.R2F=1/2e-3;          
Sim.R1F=1;            
    
% CEST pool 'G'
Sim.fG=0.02;     
Sim.kGA=16;            
Sim.dwG=-3.5;          
Sim.R2G=1/4.4e-4;          
Sim.R1G=1;            

    
%MT
Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.MT_lineshape      = 'Lorentzian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.R1C=1;
Sim.fC=0.12;
Sim.dwC=-0.6;
Sim.R2C=1/7.2e-05;  % 1/9.1µs
Sim.kCA=20; 
Sim.kAC=Sim.kCA*Sim.fC;

Sim.FREQ=500;
Sim.Zi=1.;
Sim.B1=1;		 % the saturation B1 in µT
Sim.tp=2;		 % [s]
Sim.n=1;	
Sim.DC=1;
Sim.shape='block';		 
Sim.pulsed=0;		

Sim.xZspec = -5:0.1:5;


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
 %%
ana = ANALYTIC_SIM(Sim);
figure(1), plot(ana.x,ana.zspec,'k.-'); hold on

%%
path = "/Users/cbie1/Documents/CEST_Machine_Learning_Brain/results/brain/brain_Ztab_cortex_2pts.mat";
expdata = load(path);
expdata = (expdata.Ztab_cortex_2pts)';
figure(1),plot(ana.x, expdata(:,3),'ko')
% xlim([5, -5]);