% This script defines ubinding rate constant models for bonds in the
% cadherin-based linkages used in the AJ clutch model versions. Control
% versions of the bond models (forced to be permanent or Bell Model Slip)
% are also defined for use in Validation Simulations.

%% Global Parameters
kT = 4.114; % [pN*nm] at T = 298 K

%% Two Bound State Model for CadCatCom_Factin
% Bond model and parameters from Buckley et al., Science, 2014
aCatFactin_TwoBoundStateParams.Type = 'Two Bound State Catch';
aCatFactin_TwoBoundStateParams.k0_10 = 11; % k0_ij for 1 --> 0
aCatFactin_TwoBoundStateParams.k0_20 = .14; % k0_ij for 2 --> 0
aCatFactin_TwoBoundStateParams.k0_12 = 3; % k0_ij for 1 --> 2
aCatFactin_TwoBoundStateParams.k0_21 = 20; % k0_ij for 2 --> 1
aCatFactin_TwoBoundStateParams.x10 = 0; % xij for 1 --> 0
aCatFactin_TwoBoundStateParams.x20 = 0.4; % xij for 2 --> 0
aCatFactin_TwoBoundStateParams.x12 = 0.2; % xij for 1 --> 2
aCatFactin_TwoBoundStateParams.x21 = -4; % xij for 2 --> 1

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

%% Dudko-Hummer-Szabo Model for Ecad_Ecad_WT at 3.0 sec contact time and
% Bond model and parameters from Rakshit et al., PNAS, 2012
Ecad_Ecad_WT_DudkoHummerSzaboParams.Type = 'Dudko-Hummer-Szabo Slip Bond';
Ecad_Ecad_WT_DudkoHummerSzaboParams.nu = 1/2; % scaling factor that specifies the free-energy profile. nu=0.5 for harmonic energy well.
Ecad_Ecad_WT_DudkoHummerSzaboParams.x = 0.46; % distance between the bound state and the transition state
Ecad_Ecad_WT_DudkoHummerSzaboParams.DeltaG_ActivationNoForce = 20.6; % free energy of activation in absence of external force
Ecad_Ecad_WT_DudkoHummerSzaboParams.k0 = 1/0.63; % intrinsic off rate

%% Bell Model Slip Bond w/ Upper Lifetime Limit for Vcl_aCat
% Fit of 1 Pathway Bell Model to the SMFS Data from Le et al.,
% Science Advances, 2019, adding a max lifetime cutoff. Max lifeteime
% cutoff informed by the observation that Bell Model over-predicts lifetime
% at low force.
Vcl_aCat_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
Vcl_aCat_BellModelWithMaxLifetimeLimitParams.k01 = 5.1360e-08;
Vcl_aCat_BellModelWithMaxLifetimeLimitParams.x1 = 3.6771;
Vcl_aCat_BellModelWithMaxLifetimeLimitParams.k02 = 0;
Vcl_aCat_BellModelWithMaxLifetimeLimitParams.x2 = 0;
Vcl_aCat_BellModelWithMaxLifetimeLimitParams.maxLifetime = 1.7691E4;

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

%% Dudko-Hummer-Szabo Model: Permanent Bond
% This is used in validation simulations. It should behave as a permanent
% bond that cannot unbind.
PermanentBond_DudkoHummerSzaboParams.Type = 'Dudko-Hummer-Szabo Slip Bond';
PermanentBond_DudkoHummerSzaboParams.nu = 1/2; % scaling factor that specifies the free-energy profile. nu=0.5 for harmonic energy well.
PermanentBond_DudkoHummerSzaboParams.x = 0.46; % distance between the bound state and the transition state
PermanentBond_DudkoHummerSzaboParams.DeltaG_ActivationNoForce = 20.6; % free energy of activation in absence of external force
PermanentBond_DudkoHummerSzaboParams.k0 = 0; % intrinsic off rate

%% Bell Model Slip Bond w/ Upper Lifetime Limit: Permanent Bond
% This is used in validation simulations. It should behave as a permanent
% bond that cannot unbind.
PermanentBond_BellModelWithMaxLifetimeLimitParams.Type = 'Bell Model With Max Lifetime Limit';
PermanentBond_BellModelWithMaxLifetimeLimitParams.k01 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.x1 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.k02 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.x2 = 0;
PermanentBond_BellModelWithMaxLifetimeLimitParams.maxLifetime = Inf;

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
