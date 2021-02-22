%% this file allows to create fit Z-spectrum data with BM simulation
%%it also allows to fit single offset QUESP and QUEST formulaes
clear all
close all
addpath(genpath(pwd));

%% 1.1A   load XLS_qCEST file 
% this creates a stack of Z-spectra Z_x and teh parameters struct P
[w_x, Z_x, P, rowname, Ztab]=LOAD_xls_qCEST_2_Ztab({'MDA-MB-231'});
expname=Ztab{rowname,'exp'}{1};
P

%% 1.2: create parameter struct "Sim" for simulation/fit model

% simulation parameters
Sim=init_Sim(struct()); % initializes all pools with zeroes and sets some standard values
Sim.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
Sim.MT            = 1;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
Sim.MT_lineshape  = 'Lorentzian';    % ssMT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.n_cest_pool   = 6;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)

% MR and sequence parameters
Sim.FREQ          = 11.7*gamma_;           % frequency (=B0[T] * gamma)
Sim.B1            = 0.3;                   % B1 value in µT
Sim.Trec          = 10;                    % recover time in s

Sim.Zi            = 0;        % initial magnetisation (should be between -1 and +1)
Sim.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation

if Sim.pulsed
    Sim.shape = 'seq_gauss';              % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
    Sim.n     = 20;                      % number of saturation pulses
    Sim.tp    = 0.1;                      % saturation time per pulse in s
    Sim.DC    = 0.5;                   % duty cycle
else
    Sim.tp    = 4;                       % saturation time in s
    Sim.n     = 1;                        % choose n=1 for cw saturation
    Sim.shape = 'block';                  % choose 'block' for cw saturation
    Sim.DC    = 1.0;                      % choose DC=1 for cw saturation
end;

% Pool system parameters  
% water pool A
Sim.dwA=0;
Sim.R1A=1/2.;
Sim.R2A=1/0.017;


% 1 CEST pool B (Cr)
Sim.n_cest_pool=1;
Sim.fB=0.0005;  % rel. conc 190mM/111M
Sim.kBA=150;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwB=2;     % ppm  relative to dwA
Sim.R1B=1;      % R1B relaxation rate [Hz]
Sim.R2B=1/1e-3;     % R2B relaxation rate [Hz]

% 2 CEST pool D (amide)
Sim.n_cest_pool=2;
Sim.fD=0.001;  % rel. conc 190mM/111M
Sim.kDA=50;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwD=3.5;     % ppm  relative to dwA
Sim.R1D=1;      % R1B relaxation rate [Hz]
Sim.R2D=1/1e-3;     % R2B relaxation rate [Hz]

% 3 CEST pool E (NOE, fat)
Sim.n_cest_pool=3;
Sim.fE=0.0003;  % rel. conc 35mM/111M
Sim.kEA=45;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwE=-3.6;     % ppm  relative to dwA
Sim.R1E=1;      % R1B relaxation rate [Hz]
Sim.R2E=1/3e-4;     % R2B relaxation rate [Hz]

% 4 CEST pool F (-2.3 MT)
Sim.n_cest_pool=4;
Sim.fF=0.01;  % rel. conc 35mM/111M
Sim.kFA=15;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwF=-2.3;     % ppm  relative to dwA
Sim.R1F=1;      % R1B relaxation rate [Hz]
Sim.R2F=1/7e-5;     % R2B relaxation rate [Hz]

% MT pool C
Sim.MT = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.MT_lineshape = 'Lorentzian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.R1C=1;
Sim.R2C=1/5e-05;  % 1/9.1µs
Sim.fC=0.01;
Sim.dwC=0;     % ppm  relative to dwA
Sim.kCA=15; 
Sim.kAC=Sim.kCA*Sim.fC;


%% 1.3 set and optimize start values and boundaries
Sim.n_cest_pool=4;
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL A XXXXXXXXXX
T.varyA       = [   1           1              1            ];
T.dep_varsA   = {   'dwA',      'R1A',          'R2A'          };                  
T.startA      = [   Sim.dwA       Sim.R1A       Sim.R2A        ];                  
T.lowerA      = [   Sim.dwA-0.05     0.4   Sim.R2A/2    ];
T.upperA      = [   Sim.dwA+0.05     1     Sim.R2A*2  ];

[T.dep_varsA, T.startA, T.lowerA, T.upperA] = selectVars( T.varyA, T.dep_varsA, T.startA, T.lowerA, T.upperA );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL B XXXXXXXXXX
T.varyB       = [   1              1           0             0      ];
T.dep_varsB   = {   'dwB',        'fB',       'kBA',        'R2B'    };     
T.startB      = [   Sim.dwB       Sim.fB      Sim.kBA       Sim.R2B  ];
T.lowerB      = [   Sim.dwB-0.05    Sim.fB*0.01  Sim.kBA/5     Sim.R2B/10      ];
T.upperB      = [   Sim.dwB+0.05    Sim.fB*100    Sim.kBA*5    Sim.R2B*10     ];

[T.dep_varsB, T.startB, T.lowerB, T.upperB] = selectVars( T.varyB, T.dep_varsB, T.startB, T.lowerB, T.upperB );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL C XXXXXXXXXX  always reserved for MT in numeric BM, thus commented here
T.varyC       = [ 0         1           0           0           0           ];
T.dep_varsC   = {'dwC',     'fC',       'kCA',      'R1C',      'R2C'       };     
T.startC      = [Sim.dwC    Sim.fC      Sim.kCA    Sim.R1C     Sim.R2C       ];
T.lowerC      = [-3.5          0.0001      Sim.kCA/5          Sim.R1C/2        1/10*Sim.R2C        ];
T.upperC      = [0        0.1         Sim.kCA*2       Sim.R1C*2       100*Sim.R2C       ];
[T.dep_varsC, T.startC, T.lowerC, T.upperC] = selectVars( T.varyC, T.dep_varsC, T.startC, T.lowerC, T.upperC );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL D XXXXXXXXXX
T.varyD       = [ 1            1           0                0            ];
T.dep_varsD   = {'dwD',       'fD',       'kDA',          'R2D',         };     
T.startD      = [Sim.dwD      Sim.fD        Sim.kDA       Sim.R2D        ];
T.lowerD      = [Sim.dwD-0.05   Sim.fD*0.01    Sim.kDA/2     Sim.R2D/10            ];
T.upperD      = [Sim.dwD+0.05   Sim.fD*100     Sim.kDA*2    Sim.R2D*100            ];
 
[T.dep_varsD, T.startD, T.lowerD, T.upperD] = selectVars( T.varyD, T.dep_varsD, T.startD, T.lowerD, T.upperD );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL E XXXXXXXXXX
T.varyE       = [ 1            1           0                1            ];
T.dep_varsE   = {'dwE',       'fE',       'kEA',          'R2E',         };     
T.startE      = [Sim.dwE      Sim.fE        Sim.kEA       Sim.R2E        ];
T.lowerE      = [Sim.dwE-0.05   Sim.fE*0.01    Sim.kEA/3     Sim.R2E/10            ];
T.upperE      = [Sim.dwE+0.05   Sim.fE*100     Sim.kEA*3    Sim.R2E*100            ];

[T.dep_varsE, T.startE, T.lowerE, T.upperE] = selectVars( T.varyE, T.dep_varsE, T.startE, T.lowerE, T.upperE );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL F XXXXXXXXXX
T.varyF       = [ 0            1           0                0            ];
T.dep_varsF   = {'dwF',       'fF',       'kFA',          'R2F',         };     
T.startF      = [Sim.dwF      Sim.fF        Sim.kFA       Sim.R2F        ];
T.lowerF      = [Sim.dwF-0.1   Sim.fF*0.01    Sim.kFA/2     Sim.R2F/10            ];
T.upperF      = [Sim.dwF+0.1   Sim.fF*100     Sim.kFA*2    Sim.R2F*100            ];
 
[T.dep_varsF, T.startF, T.lowerF, T.upperF] = selectVars( T.varyF, T.dep_varsF, T.startF, T.lowerF, T.upperF );

FIT.T=T;   clear T;  % store start values and boundaries in FIT struct
FIT.Sim=Sim; % store Sim values in FIT struct

figure(2002), multiZplot(P,Sim,FIT.T,w_x,Z_x); % plot guess  with data

% Zknown=0.95; % this can be added to guess the initial magnetization Zi
% P.Zi= 1 - (1-Zknown)*exp(-P.R1A*P.Trec);

%% 1.4  multi-Z-BMfitting of the Z_x data (wherever it comes from)
% RUN full BM OPTIMIZATION 
% you need a startvalue, run 1.3 first!
Sim.analytic=1;  % set1 this to 1 if analytic fit should be used, numeric =0 can take forever
Sim.n_cest_pool=4;

% fit-options
[FIT] =multiZfit(P,Sim,FIT.T,w_x,Z_x);

figure,
FIT.Z_fit=multiZplot(P,Sim,FIT.T,w_x,Z_x,FIT.popt,FIT.pci);
title([expname ' : ' rowname ', R2=' num2str(FIT.R2,4)]);
% savefig(['BM_FIT_3p' rowname '.fig']);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

FIT.P=P;
FIT.w_x=w_x;
FIT.Z_x=Z_x;

if exist('Ztab') % save fitresult in Ztable
Ztab(rowname,'FIT3p')={{FIT}}; % name the fit for saving in Ztab
end;



