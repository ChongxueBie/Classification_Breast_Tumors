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
Sim.dwB=0.9; % ppm  relative to dwA
Sim.R1B=1; 
Sim.R2B=55e-3;

% CEST pool 'D'
Sim.fD=0.00146;     % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
Sim.kDA=96;            % exchange rate [s^-1]     % corresponds to creatine at ~ 22°C, pH=6.4 in PBS (Goerke et al.)
Sim.dwD=2;          % chemical shift of the CEST pool in [ppm] 
Sim.R2D=1/2.7e-3;           % transversal relaxation rate 1/T2 of pool b  [s^-1] 
Sim.R1D=1;            % longitudinal relaxation rate 1/T1 of pool b  [s^-1] 

% CEST pool 'E'
Sim.fE=0.0125;     % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
Sim.kEA=5500;            % exchange rate [s^-1]     % corresponds to creatine at ~ 22°C, pH=6.4 in PBS (Goerke et al.)
Sim.dwE=2.8;          % chemical shift of the CEST pool in [ppm] 
Sim.R2E=1/1.65e-4;           % transversal relaxation rate 1/T2 of pool b  [s^-1] 
Sim.R1E=1;            % longitudinal relaxation rate 1/T1 of pool b  [s^-1] 

% CEST pool 'F'
Sim.fF=0.0195;     % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
Sim.kFA=30;            % exchange rate [s^-1]     % corresponds to creatine at ~ 22°C, pH=6.4 in PBS (Goerke et al.)
Sim.dwF=3.4;          % chemical shift of the CEST pool in [ppm] 
Sim.R2F=1/2.2e-4;           % transversal relaxation rate 1/T2 of pool b  [s^-1] 
Sim.R1F=1;            % longitudinal relaxation rate 1/T1 of pool b  [s^-1] 
    
% CEST pool 'G'
Sim.fG=0.026;     % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
Sim.kGA=16;            % exchange rate [s^-1]     % corresponds to creatine at ~ 22°C, pH=6.4 in PBS (Goerke et al.)
Sim.dwG=-3.2;          % chemical shift of the CEST pool in [ppm] 
Sim.R2G=1/2e-4;           % transversal relaxation rate 1/T2 of pool b  [s^-1] 
Sim.R1G=1;            % longitudinal relaxation rate 1/T1 of pool b  [s^-1] 

    
%MT
Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.MT_lineshape      = 'Lorentzian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.R1C=1;
Sim.fC=0.016;
Sim.R2C=7.2e-05;  % 1/9.1µs
Sim.kCA=40; 
Sim.kAC=Sim.kCA*Sim.fC;

Sim.FREQ=500;
Sim.Zi=0.;
Sim.B1=1.0;		 % the saturation B1 in µT
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
Sim.dummies=0; 
Sim.flipangle=5; 
% Sim.readout='FID'; 
Sim.spoilf=0; 



num = NUMERIC_SIM(Sim);
figure(1), plot(num.x,num.zspec,'.'); hold on;
ana = ANALYTIC_SIM(Sim)
figure(1), plot(ana.x,ana.zspec); 


% plot(ana.x,ana.zspec,ana.x,ana.zspec(end:-1:1)-ana.zspec ); hold on;


% expname='PARACEST'
% 
% if ~exist('expn')
%     expn=1;
% end;

% T = table(num.x,num.zspec,num.zspec(end:-1:1)-num.zspec);
% T.Properties.VariableNames = {'offsets' 'Z' 'Asym'}
% filename = 'simulation_data.xlsx';
% writetable(T,filename,'Sheet',sprintf('exp%d_%s',expn,expname));
% expn=expn+1;
