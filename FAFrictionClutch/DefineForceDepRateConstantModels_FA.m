% This script defines ubinding rate constant models for bonds in the
% integrin-based linkages used in the FA clutch model versions. Control
% versions of the bond models (forced to be permanent or Bell Model Slip)
% are also defined for use in Validation Simulations.

%% Global Parameters
kT = 4.114; % [pN*nm] at T = 298 K

%% Two Bound State Model for Vcl_Factin Pulling Factin in (-) Direction
% Huang et al. 2017, T12 Vinculin
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
% Based on my own fit of a 1 Pathway Bell Model to the SMFS Data from Le et
% al. 2019. Max lifeteime cutoff informed by the observation that Bell
% Model over-predicts lifetime at low force. Max lifetime chosen as the
% experimental lifetime at the lowest force (~1 pN), which is conservative.
% The detials for this bond a low force do not matter very much because
% this bond is much longer lived than all other bonds in question.
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.k01 = 5.5773E-07;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.x1 = 2.3164;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.k02 = 0;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.x2 = 0;
Vcl_Talin_BellModelWithMaxLifetimeLimitParams.maxLifetime = 5311.0881;

%% Bell Model (Single Exponential) Catch Bond w/ Upper Lifetime Limit for TalinABS3_Factin Pulling Factin in (-) Direction
% Based on single exponential fit to bond lifetimes from 0 to 10 pN given
% in the caption to Fig. S9 from Owens et al. PNAS 2022. Max lifetime
% cutoff informed by the observation that the bond lifetime plateaus after
% 10 pN. Max lifetime chosen as the lifetime at the max force probed by
% Owens et al., i.e. {16 pN, 89.3450 sec}.
% a = 0.0476 [s] --> k01 = 1/a = 21.0084
% b = 0.6663 [1/pN] --> x1 = -kBT/b = -4.114 pN*nm * 0.6663/pN = -2.7412 [nm]
% Note: The detials for this bond at high force do not matter very much because
% this bond is much longer lived than all other bonds in question.
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.k01 = 21.0084;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.x1 = -2.7412;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.k02 = 0;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.x2 = 0;
TalinABS3_Factin_NegPull_BellModelWithMaxLifetimeLimitParams.maxLifetime = 89.3450;

%% Intg_alphaVbeta3_FN_Mn
% My 2-Pathway Catch-Slip Bond Fit to data from Elosegui-Artola et al. NCB 2016, Fig 3A
% Performed in: MAIN_Fit2PathwayBellModel.m
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.k01 = 0.347120807;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.x1 = -0.037418289;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.k02 = 2.86866E-08;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.x2 = 1.566854962;
Intg_alphaVbeta3_FN_BellModelWithMaxLifetimeLimitParams.maxLifetime = Inf;

%% XXX Below Control Versions Of The Different Bond Model Types Are Defined XXX
% Control versions of the bond models (forced to be permanent or Bell Model
% Slip) are defined for use in Validation Simulations.

%% Two Bound State Model: Permanent Bond
% This is used in validation simulations. It should behave as a permanent
% bond that cannot unbind.
PermanentBond_TwoBoundStateParams.Type = 'Two Bound State Catch';
PermanentBond_TwoBoundStateParams.k0_10 = 0; % k0_ij for 1 --> 0
PermanentBond_TwoBoundStateParams.k0_20 = 0; % k0_ij for 2 --> 0
PermanentBond_TwoBoundStateParams.k0_12 = 0; % k0_ij for 1 --> 2
PermanentBond_TwoBoundStateParams.k0_21 = 0; % k0_ij for 2 --> 1
PermanentBond_TwoBoundStateParams.x10 = 0; % xij for 1 --> 0
PermanentBond_TwoBoundStateParams.x20 = 0; % xij for 2 --> 0
PermanentBond_TwoBoundStateParams.x12 = 0; % xij for 1 --> 2
PermanentBond_TwoBoundStateParams.x21 = 0; % xij for 2 --> 1

%% Two Bound State Model: Bell Model Slip Bond Control
% This is used in validation simulations. It should reproduce a Bell Model
% slip bond with parameters k0 = 0.1 and Fb = 2 pN.
ControlSlipBond_TwoBoundStateParams.Type = 'Two Bound State Catch';
ControlSlipBond_TwoBoundStateParams.k0_10 = 0.1; % k0_ij for 1 --> 0
ControlSlipBond_TwoBoundStateParams.k0_20 = 0; % k0_ij for 2 --> 0
ControlSlipBond_TwoBoundStateParams.k0_12 = 0; % k0_ij for 1 --> 2
ControlSlipBond_TwoBoundStateParams.k0_21 = 0; % k0_ij for 2 --> 1
ControlSlipBond_TwoBoundStateParams.x10 = 1/2*kT; % xij for 1 --> 0
ControlSlipBond_TwoBoundStateParams.x20 = 0; % xij for 2 --> 0
ControlSlipBond_TwoBoundStateParams.x12 = 0; % xij for 1 --> 2
ControlSlipBond_TwoBoundStateParams.x21 = 0; % xij for 2 --> 1

%% Dudko-Hummer-Szabo Model: Permanent Bond
% This is used in validation simulations. It should behave as a permanent
% bond that cannot unbind.
PermanentBond_DudkoHummerSzaboParams.Type = 'Dudko-Hummer-Szabo Slip Bond';
PermanentBond_DudkoHummerSzaboParams.nu = 1/2; % scaling factor that specifies the free-energy profile. nu=0.5 for harmonic energy well.
PermanentBond_DudkoHummerSzaboParams.x = 0.46; % distance between the bound state and the transition state
PermanentBond_DudkoHummerSzaboParams.DeltaG_ActivationNoForce = 20.6; % free energy of activation in absence of external force
PermanentBond_DudkoHummerSzaboParams.k0 = 0; % intrinsic off rate

%% Dudko-Hummer-Szabo Model: Bell Model Slip Bond Control
% This is used in validation simulations. It should reproduce a Bell Model
% slip bond with parameters k0 = 0.1 and Fb = 2 pN.
% Dudko-Hummer-Szabo Model --> Bell Model as v-->1 and DeltaG-->Infinity
ControlSlipBond_DudkoHummerSzaboParams.Type = 'Dudko-Hummer-Szabo Slip Bond';
ControlSlipBond_DudkoHummerSzaboParams.nu = .99999; % scaling factor that specifies the free-energy profile. nu=0.5 for harmonic energy well.
ControlSlipBond_DudkoHummerSzaboParams.x = 1/2*kT; % distance between the bound state and the transition state
ControlSlipBond_DudkoHummerSzaboParams.DeltaG_ActivationNoForce = 1E8; % free energy of activation in absence of external force
ControlSlipBond_DudkoHummerSzaboParams.k0 = 0.1; % intrinsic off rate

%% Bell Model Slip Bond w/ Upper Lifetime Limit: Permanent Bond
% This is used in validation simulations. It should behave as a permanent
% bond that cannot unbind.
PermanentBond_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
PermanentBond_BellModelWithMaxLifetimeLimitParams.k01 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.x1 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.k02 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.x2 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.maxLifetime = Inf;