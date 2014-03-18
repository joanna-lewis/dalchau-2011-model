%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix random numbers for generating data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randn('state', 1);
rand('twister', 1);


%%%%%%%%%%%%%%%
% Set options %
%%%%%%%%%%%%%%%

Options.Burnin                = 1000;
Options.NumOfPosteriorSamples = 1;


% Set name for saving results
Options.EquationName            = 'dalchau';

Options.NumOfParameters         = 19;
Options.ParameterOrderOfMagnitude = transpose([-4, -6, -7, ...
                                    4, 4, 4, ...
                                    -6, 4, -9, ...
                                    -8, -1, 2, ...
                                    -1, 2, -9, ...
                                    -5, -3, -5, ...
                                    3, ...
                                    0]); % scale factor and sd (error) for epitope data
Options.ObservedViralSpecies    = [1];
Options.UnobservedViralSpecies  = [2];
Options.SelfSpecies             = [3];
Options.OtherSpecies            = [4:16];

Options.SDNoiseAdded            = 10;

Options.SaveMetricTensors       = true;

D             = 25;
StartTime     = 0;
EndTime       = 12.5;
TimePoints    = (StartTime:EndTime/(D-1):EndTime)';

%%%%%%%%%%%%%%%%%
% Generate Data %
%%%%%%%%%%%%%%%%%

n_o   = length(Options.ObservedViralSpecies);
n_u   = length(Options.UnobservedViralSpecies);
n_s   = length(Options.SelfSpecies);
n_tot = n_o + n_u + n_s; % total number of peptide species

Parameters    = transpose([8.764e-4, 5.658e-6, 4.177e-7,... %(u)
                    2.093e4, 1.759e4, 1.064e4,...    %(g)
                    1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
                    8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
                    0.1142, 150.5, 1.663e-9,...            %(e, g_M, b_T)
                    7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
                    1505]);                          % g_T

Options.InitialValues = dalchau_model_findss(...
                    Parameters, ...
                    n_o, ...
                    n_u, ...
                    n_s);
                
%%%%%%%%%%%%%%%%%%%%%%%%%        
% Read in data          %
%%%%%%%%%%%%%%%%%%%%%%%%%


NoisyData = csvread('simulate_data_111213_2.csv'); % same as simulate_data_111213_2.csv, but tryptic peptide levels set to zero

Options.gdata = [TimePoints'; NoisyData(1:3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%        
% Call sampling routine %
%%%%%%%%%%%%%%%%%%%%%%%%%

Options.MathParToInfer = [1, 2, 3]; % indices of the parameters in the mathematical model you're interested in 
Options.StatParToInfer = [1]; % indices of the parameters in the statistical model you're interested in 



ODE_mMALASimp_scaled(NoisyData(4,:), TimePoints, Options);


