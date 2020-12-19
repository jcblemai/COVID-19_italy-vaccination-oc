
load('results_true/state_par20210101.mat')
% Show dates
datestr(time_b)
% choose realization:
i = 30;


%%% things to check:
PAR = PAR_real(i,:);
Vt = V;
V.gammaQ_over_gammaH = Vt.gammaQgammaH;
V.gammaA_over_gammaQ = Vt.gammaAgammaQ;
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


%% Things to do:
ndays = 164;
% ...
x = zeros(ndays, 1177); % x should contains, for each days: S1 S2 S3 ... S107 E1 ... E107 ... of this realization i
V.zeta = 0; % This should contains zeta;
V.x0 = x0_real_out(i,:);  % contains x0 of this i, in the same format S1 .. S107 E1 ...E107 ...
% Beta: IN SEPIA.m there was:
% ---> Calculate  trasmission ratio (beta_ratio) for each node (1st dimension) and day of
% ---> V.time_model (second dimension) through linear interpolation of V.mob
% ---> on = ones(V.n,1);
% ---> betaP3P2_reg = PAR(V.nPAR_model+1:V.nPAR_model+V.n_reg)'; %regional value
% ---> betaP3P2 = betaP3P2_reg(V.prov_IDreg); %province value
% ---> beta_p = eval(V.beta_string);
% ---> beta_ratio = interp1(V.tbeta,beta_p',V.time_model)';
% i'd need this betaratio here, but how to get it.

%% Questions:
% is V.p a vector I shoud reapeat to the end of time ? 


