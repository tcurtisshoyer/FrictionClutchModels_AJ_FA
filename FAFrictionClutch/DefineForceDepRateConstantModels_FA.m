% This script defines ubinding rate constant models for bonds in the
% integrin-based linkages used in the FA clutch model versions.

%% Global Parameters
kT = 4.114; % [pN*nm] at T = 298 K

%% Two Bound State Model for Vcl_Factin Pulling Factin in (-) Direction
% Bond model and parameters from Huang et al., Science, 2017, T12 Vinculin
Vcl_Factin_NegPull_TwoBoundStateParams.Type = 'Two Bound State Catch';
Vcl_Factin_NegPull_TwoBoundStateParams.k0_10 = 5.3; % k0_ij for 1 --> 0
Vcl_Factin_NegPull_TwoBoundStateParams.k0_20 = 5.5E-3; % k0_ij for 2 --> 0
Vcl_Factin_NegPull_TwoBoundStateParams.k0_12 = 6.1; % k0_ij for 1 --> 2
Vcl_Factin_NegPull_TwoBoundStateParams.k0_21 = 43; % k0_ij for 2 --> 1
Vcl_Factin_NegPull_TwoBoundStateParams.x10 = 0; % xij for 1 --> 0
Vcl_Factin_NegPull_TwoBoundStateParams.x20 = 1.2; % xij for 2 --> 0
Vcl_Factin_NegPull_TwoBoundStateParams.x12 = 0.4; % xij for 1 --> 2
Vcl_Factin_NegPull_TwoBoundStateParams.x21 = -3.4; % xij for 2 --> 1   

%% Bell Model Slip Bond w/ Upper Lifetime Limit for Vcl_Talin
% Based fit of a 1 Pathway Bell Model to the SMFS Data from Le et al.,
% Science Advances, 2019, adding a max lifetime cutoff. Max lifeteime
% cutoff informed by the observation that Bell Model over-predicts lifetime
% at low force.
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.k01 = 5.5773E-07;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.x1 = 2.3164;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.k02 = 0;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.x2 = 0;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.maxLifetime = 5311.0881;

%% Bell Model (Single Exponential) Catch Bond w/ Upper Lifetime Limit for TalinABS3_Factin Pulling Factin in (-) Direction
% Bond model and parameters from Owens et al., PNAS, 2022. Max lifetime
% cutoff informed by the observation that the bond lifetime plateaus after
% 10 pN.
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.k01 = 21.0084;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.x1 = -2.7412;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.k02 = 0;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.x2 = 0;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.maxLifetime = 89.3450;

%% 2-Pathway Catch-Slip Bond Fit for Intg_alphaVbeta3_FN_Mn
% Fit of 2-Pathway Catch-Slip Bond Model to data from Elosegui-Artola et al., Nat. Cell Biol., 2016.
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.k01 = 0.347120807;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.x1 = -0.037418289;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.k02 = 2.86866E-08;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.x2 = 1.566854962;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.maxLifetime = Inf;
