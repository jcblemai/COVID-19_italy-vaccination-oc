clear all
close all

addpath('./input_data')
addpath('./functions')

rng(0)

%% READ POPULATION DATA

load Population_data.mat

% Content

% N    : Population of each Italian provice
% n    : number of Italian provinces
% n_reg: number of Italian regions
% p    : mobile population fraction for each Italian province
% q    : n*n mobility matrix. Element qij represents the fraction of the
%        mobile population living in province i that moves to province j
% prov_name  : name of the provinces
% reg_name   : name of the regions
% prov_IDprov: Census ID of the province
% prov_IDreg : Census ID of the region each province belongs to
% prov2reg   : n_reg*n matrix. Element ij is equal to 1 if province j 
%              belongs to region i

% Transfer data into V structure
V=Population_data;
clear Population_data

%% READ EPIDEMIOLOGICAL DATA
load('Epi_data.mat')

% Content

% Date     : date of epidemiologica data. Matlab datenum format
% prov_Hnew: n*length(Date) Matrix, reconstructed daily new hospitalized
%            cases for each province
% prov_Hnew: n_reg*length(Date) Matrix, reconstructed daily new 
%            hospitalized cases for each region

% Transfer data into V structure
V=append_struct(V,Epi_data);
clear Epi_data

%% READ MOBILITY REDUCTION DATA

load MobilityReduction_data.mat

%Content

% Information are extracted from preprint
% https://www.medrxiv.org/content/10.1101/2020.03.22.20039933v2 and
% represents the reduction in extraprovice mobility during the lockdown
% in Italy.

% tmob: Date (Matlab datenum format) for which values of mobility are 
%       provided
% mob : n*length(tmob) matrix. Element ij represents the ratio of
%       extra-province mobility of provice i at time tmob(j) with respect 
%       to the pre-COVID mobility
%
% Timeserire of mobility is then obtained through linear interpolation of
% the values provided in tmob and mob.

% Transfer data into V structure
V=append_struct(V,MobilityReduction_data);
clear MobilityReduction_data

%% READ MODEL PARAMETERS

load Parameters.mat

% Content

% PAR_names        : Parameter names (see preprint text)
% Posteriror_sample: Posterior samples of PAR_names parameters (note that 
%                    some parameters have a fixed values, see main text)
% nPAR_model       : number of structural model parameters (first 
%                    nPAR_model parameters of PAR_names)
% x0               : initial conditions (1 exposed individual in the 
%                    province of Lodi)
% seeding          : indexes of province for which initial conditions are 
%                    also estimated

% Additional fixed parameters
% gammaA_over_gammaQ: ratio \gamma_A/\gamma_Q
% gammaQ_over_gammaH: ratio \gamma_Q/\gamma_H
% zeta: fraction of cases that are quarantined (see preprint text)

% Transfer data into V structure
V=append_struct(V,Parameters);
clear Parameters


%% SIMULATION SET UP

NSample=1; %number of sample from  posterior distribution of parameters
NNoise4sample=1; %number of sample from the negative binomial distribution for each posterior sample
NReal=NSample*NNoise4sample; %total number of realizations
quant=[0.025 0.25 0.5 0.75 0.975]; %quantiles of interests

V.time_model_final=datenum('6/30/2020')+1; %simulation horizon
betaInc=1;  %increase in trasmission after May 3, 2020 (Figure 2 of the preprint is produced using value 1, 1.2 and 1.4)

% Apply increase in trasmission in all provinces
V.ProvBetaInc=ones(V.n,1);
V.ProvBetaInc=betaInc;

% Define schedule of variation for beta transmission parameters
V.tbeta=zeros(1,9);
V.tbeta(1)=datenum('1/01/2020');
V.tbeta(2)=datenum('2/24/2020');
V.tbeta(3)=datenum('2/26/2020');
V.tbeta(4)=datenum('3/08/2020');
V.tbeta(5)=datenum('3/11/2020');
V.tbeta(6)=datenum('3/22/2020');
V.tbeta(7)=datenum('5/04/2020');
V.tbeta(8)=datenum('5/07/2020');
V.tbeta(9)=datenum('12/31/2020');

% String command to be used to define beta variations
V.beta_string='[on, on, betaP1P0*on, betaP1P0*on, betaP1P0*betaP2P1*on, betaP1P0*betaP2P1*betaP3P2.*on, betaP1P0*betaP2P1*betaP3P2.*on, betaP1P0*betaP2P1*betaP3P2.*on.*V.ProvBetaInc, betaP1P0*betaP2P1*betaP3P2.*on.*V.ProvBetaInc ]';

% Time_resample: time at which output are provided
time_resample=V.Date(1):V.time_model_final-1;

% Define colorscheme
V.col=get(groot,'defaultAxesColorOrder');

%% MAIN

% Preallocation of variables
reg_Hnew_real=zeros(V.n_reg,length(time_resample)-1,NReal); %simulated new hospitalized cases for each region (1st dimension), day (2nd dimension) and realization (third dimension)

cont_real=0; %inizialization of realization counter

cont_sample = 1;

display([' Posterior sample ',num2str(cont_sample),' of ',num2str(NSample)]);
% Run model for this posterior sample
    PAR=V.Posterior_sample(ceil(rand*length(V.Posterior_sample)),:); %select one random posterior sample
    [x,V]=SEPIA(PAR,V); %run model
    
    % from SEPIA.m
n=V.n;  %number of nodes
deltaE  = 1/PAR(2);
deltaP  = 1/PAR(3);
sigma   = PAR(4);
eta     = 1/PAR(5);
gammaI  = 1/PAR(6);
alphaI  = 1/PAR(7);
alphaH  = 1/PAR(7);
gammaH  = 1/PAR(6);
epsilonA= PAR(8);  %for brevity the ratio \beta_A/\beta_P is termed epsilonA
r       = PAR(9);  %the model assumes that clesses S E P A R have the same r values
Deltat0 = PAR(10);
epsilonI= PAR(11)*epsilonA; %for brevity the ratio \beta_I/\beta_P is termed epsilonI
betaP1P0= PAR(13);
betaP2P1= PAR(14);
gammaQ  = V.gammaQ_over_gammaH*gammaH;
gammaA  = V.gammaA_over_gammaQ*gammaQ;

% Parameter betaP0 is expressed as a function of the local reproductive
% number R_0 and the other relevnt parameters
betaP0 = PAR(1)/(1/deltaP + epsilonI*sigma/(gammaI + alphaI + eta) + epsilonA*(1-sigma)/gammaA);

% Add initial conditions in the exposed class for the seeding nodes
V.x0(n+V.seeding) = 10.^(PAR(V.nPAR_model+V.n_reg+1:V.nPAR_model+V.n_reg+length(V.seeding)));

% Define simulation time (time_model)
V.time_model = (V.Date(1)-(Deltat0)):V.time_model_final;

% Calculate mobility ratio for each node (1st dimension) and day of
% V.time_model (second dimension) through linear interpolation of V.mob
mob_ratio = interp1(V.tmob, V.mob',V.time_model)';

% Calculate  trasmission ratio (beta_ratio) for each node (1st dimension) and day of
% V.time_model (second dimension) through linear interpolation of V.mob
on = ones(V.n,1);
betaP3P2_reg = PAR(V.nPAR_model+1:V.nPAR_model+V.n_reg)'; %regional value
betaP3P2 = betaP3P2_reg(V.prov_IDreg); %province value
beta_p = eval(V.beta_string);
beta_ratio = interp1(V.tbeta,beta_p',V.time_model)';

%[t,x] = ode45(@eqs,V.time_model,V.x0);

% END SEPIA.m
    
    prov_cumH=interp1(V.time_model,x(:,10*V.n+1:11*V.n),time_resample)'; %cumulative hospitalized cases for each province (1st dimension) and day (2nd dimension) of time_resample
    prov_Hnew=diff(prov_cumH,1,2); %new hospitalized cases for each province (1st dimension) and day (2nd dimension) of time_resample
    
