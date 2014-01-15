function [] = ODE_mMALASimp( Y, TimePoints, Options )

% Start the timer
tic

rand('twister', sum(100*clock))
randn('state', sum(100*clock))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise user options if not already specified                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get N (number of chemical species) and D (number of time points)
[N, D] = size(Y);
NumOfSpecies = N;
NumOfObs     = D;

% Get specified options

MathParToInfer        = Options.MathParToInfer;
StatParToInfer        = Options.StatParToInfer;
ParametersToInfer     = [MathParToInfer, numel(MathParToInfer) + StatParToInfer];

NumOfParameters       = numel(ParametersToInfer);

StepSize = 10^(-2);

SpeciesObserved       = Options.ObservedViralSpecies;

SDNoiseAdded          = Options.SDNoiseAdded;


Burnin                = Options.Burnin;
NumOfPosteriorSamples = Options.NumOfPosteriorSamples;

EquationName          = Options.EquationName;

ODEoptions         = odeset('RelTol',1e-6,'AbsTol',1e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise non changeable stuff                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup initial parameters - *including* proportionality constant for
% epitope data
Parameters = transpose([8.764e-4, 5.658e-6, 4.177e-7,... %(u)
                    2.093e4, 1.759e4, 1.064e4,...    %(g)
                    1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
                    8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
                    0.1142, 150.5, 1.663e-9,...            %(e, g_M, b_T)
                    7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
                    1505, ...
                    1]);               % proportionality for epitope data

% Set up proposal counters
AcceptedMutation  = zeros(size(Parameters));
AttemptedMutation = zeros(size(Parameters));

% Set up parameter history variable
ParaHistory         = zeros(NumOfPosteriorSamples, NumOfParameters);
LLHistory           = zeros(NumOfPosteriorSamples, NumOfSpecies);


% Set up initial noise for likelihood
% Fix noise - CurrentNoise is the variance
CurrentNoise = (1.9768)^2;

% Set monitor rate for adapting step sizes
MonitorRate = 100;

% Set up converged flag
ContinueIterations = true;
Converged          = false;

% Precalculate the log likelihood
    
    n_o   = length(Options.ObservedViralSpecies);
    n_u   = length(Options.UnobservedViralSpecies);
    n_s   = length(Options.SelfSpecies);
    n_tot = n_o + n_u + n_s; % total number of peptide species
    
disp('Initialisation Completed..');

% Initialise iteration number
IterationNum = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate gradient, metric, etc. at initial point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = ParametersToInfer;
        % calculate stuff with old parameters - could speed this up by
        % remembering from previous iteration
        [OldXEstimates, OldXSens] = ode_model_sol(Parameters(1:end-1), ...
                        n_tot, ...
                        TimePoints, ...
                        Options.InitialValues, ...
                        a);
        
        OldLL       = LogNormPDF(Y(2:end), Parameters(end)*OldXEstimates(3*n_tot+1,:), CurrentNoise);
        OldLogPrior = log(ModelParameterPrior(a, Parameters(a)));
        
        OldGradL = sum( ...
                     (Y(2:end) - OldXEstimates(3*n_tot+1,:)) ...
                     .* OldXSens(3*n_tot+1,:) ...
                     /CurrentNoise ...
                     );
        OldGradL = OldGradL + ModelParameterLogPriorDerivative(a, Parameters(a));
        
        OldG = zeros(NumOfParameters);
        for i = 1:NumOfParameters
            for j = i:NumOfParameters
            
                OldG(i,j) = (1/CurrentNoise)*...
                    (OldXSens(3*n_tot+1,:,i)*OldXSens(3*n_tot+1,:,j)');
                
                OldG(j,i) = OldG(i,j);
                
            end
        end
        % add prior to Fisher Information
        
        OldG    = OldG - diag(ModelParameterLogPriorDerivative2(ParametersToInfer, Parameters(ParametersToInfer))); 
        
        OldGInv = inv(OldG + eye(NumOfParameters)*1e-6);
        
        OldMean = Parameters(a) + OldGInv*OldGradL*StepSize/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
while ContinueIterations
    
    % Increment iteration number
    IterationNum = IterationNum + 1;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate parameter values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each parameter, sample new value and accept with metropolis ratio
    for a = [ParametersToInfer]
        
        % make a proposal
        AttemptedMutation(a)  = AttemptedMutation(a) + 1;
        
        NewParas    = Parameters;
        NewParas(a) = OldMean + randn(1)*chol(StepSize*OldGInv);
        
       % NewParas(a) - Parameters(a)
        
        % calculate stuff with proposed parameters
        try
            [NewXEstimates, NewXSens] = ode_model_sol(NewParas(1:end-1), ...
                        n_tot, ...
                        TimePoints, ...
                        Options.InitialValues, ...
                        a);
        
            NewLL = LogNormPDF(Y(2:end), NewParas(end)*NewXEstimates(3*n_tot+1,:), CurrentNoise);
            
            NewGradL = sum( ...
                         (Y(2:end) - NewXEstimates(3*n_tot+1,:)) ...
                         .* NewXSens(3*n_tot+1,:) ...
                         /CurrentNoise ...
                         );
            NewGradL = NewGradL + ModelParameterLogPriorDerivative(a, NewParas(a));
            
            NewG = zeros(NumOfParameters);
            for i = 1:NumOfParameters
                for j = i:NumOfParameters

                    NewG(i,j) = (1/CurrentNoise)*...
                        (NewXSens(3*n_tot+1,:,i)*NewXSens(3*n_tot+1,:,j)');

                    NewG(j,i) = NewG(i,j);

                end
            end
        % add prior to Fisher Information
        NewG = NewG - diag(ModelParameterLogPriorDerivative2(ParametersToInfer, NewParas(ParametersToInfer))); 
        
        NewGInv = inv(NewG + eye(NumOfParameters)*1e-6);
            
        NewMean = NewParas(a) + NewGInv*NewGradL*StepSize/2;
        
        catch
            NewLL = -1e300;
            NewMean = 0;
            disp('bad proposal')
        end
        
        NewLogPrior = log(ModelParameterPrior(a, NewParas(a)));
        
        % and p(old|new) and p(new|old)
      
        LogProbNewGivenOld = -sum(log(diag(chol(StepSize*OldGInv)))) - 0.5*(OldMean-NewParas(a))'*(OldG/StepSize)*(OldMean-NewParas(a));
        
        LogProbOldGivenNew = -sum(log(diag(chol(StepSize*NewGInv)))) - 0.5*(NewMean-Parameters(a))'*(NewG/StepSize)*(NewMean-Parameters(a));
        
        
        % accept or reject
        
%         NewLL  - OldLL
%         NewLogPrior - OldLogPrior
%         LogProbOldGivenNew - LogProbNewGivenOld
        
         Ratio = NewLL + NewLogPrior + LogProbOldGivenNew ...
                    - OldLL - OldLogPrior - LogProbNewGivenOld;
            
            if Ratio > 0 || (log(rand) < Ratio)
                % Accept proposal
                % Update variables
                Parameters                 = NewParas;
                AcceptedMutation(a)        = AcceptedMutation(a) + 1;
                OldLL                      = NewLL;
                OldLogPrior                = NewLogPrior;
                OldG                       = NewG;
                OldGInv                    = NewGInv;
                OldMean                    = NewMean;
            end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%
    % Save parameters %
    %%%%%%%%%%%%%%%%%%%
%     if Converged
%         ParaHistory(IterationNum-ConvergenceIterationNum, :) = Parameters(Options.ParametersToInfer);
%         LLHistory(IterationNum-ConvergenceIterationNum, :)   = CurrentLL;
%     end
%    if Converged
        ParaHistory(IterationNum, :) = Parameters(ParametersToInfer);
        LLHistory(IterationNum, :)   = OldLL;
%    end
    
    
    % If not yet converged...
    if Converged == false
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust proposal widths %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Adjust parameter proposal widths
        if mod(IterationNum, MonitorRate) == 0
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['Iteration ' num2str(IterationNum)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(' ')
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust proposal width for parameter value inference %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     if AcceptedMutation/AttemptedMutation < 0.4
%                         StepSize = StepSize * 0.9;
%                     elseif AcceptedMutation(a)/AttemptedMutation(a) > 0.7
%                         StepSize = StepSize * 1.1;
%                     end
                    
                    disp([num2str(100*AcceptedMutation(a)/AttemptedMutation(a)) '% mutation acceptance for parameter ' num2str(a)]);
                            
            disp(' ')
            disp('Parameter step size:')
            disp(StepSize)
            
            % Reset counters
            AttemptedMutation = zeros(size(AttemptedMutation));
            AcceptedMutation  = zeros(size(AcceptedMutation));
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate convergence every 1000 steps %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Change converged tab if converged
        if IterationNum >= Burnin
            Converged               = true;
            ConvergenceIterationNum = IterationNum;
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(['Converged at iteration number ' num2str(ConvergenceIterationNum)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            
            BurnInTime = toc;
            tic;
            
        end
        
        
        
    else % Converged so decide how long to sample from posteriors
        
        if IterationNum == ConvergenceIterationNum + NumOfPosteriorSamples
            % 5000 posterior samples have been collected so stop
            ContinueIterations = false;
        end
        
    end
end   


PosteriorTime = toc;

CurTime = fix(clock);
RandTime = ceil(rand*10000);

% Save posterior
FileName = ['ODE_mMALASimp_' EquationName '_' num2str(D) 'DPS_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '_' num2str(RandTime)];
save(['./Results/' FileName], 'ParaHistory', 'LLHistory', 'StepSize', 'BurnInTime', 'PosteriorTime', 'Y', 'TimePoints');



end
