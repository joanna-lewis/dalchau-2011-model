function [] = ODE_mMALASimp_onepeptide_ss( Y, TimePoints, Options )

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

% Get specified options
n_o   = length(Options.ObservedViralSpecies);
n_u   = length(Options.UnobservedViralSpecies);
n_s   = length(Options.SelfSpecies);
n_tot = n_o + n_u + n_s; % total number of peptide species


MathParToInfer        = Options.MathParToInfer;
StatParToInfer        = Options.StatParToInfer;
ParametersToInfer     = [MathParToInfer, 2*n_tot+13 + StatParToInfer];

NumOfParameters       = numel(ParametersToInfer);

StepSize = 10^(-2); %(epsilon^2)

Burnin                = Options.Burnin;
NumOfPosteriorSamples = Options.NumOfPosteriorSamples;

EquationName          = Options.EquationName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise non changeable stuff                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup initial parameters (scaled so all of order 1) - *including* proportionality constant for
% epitope data
Parameters = transpose([1.3582, 1.3582, 1.3582,... %(u)
                    1.759, 1.759, 1.759,...    %(g)
                    0, 0, 1.663,...   %(u_T, q, b_B2705)
                    0, 1.3, 0,...              %(c, d_p, v)
                    1, 1.505, 0,...            %(e, g_M, b_T)
                    1, 1, 1, ...      % (d_M, d_T, d_Me)
                    0, ...
                    1]);         % sd(error) for epitope data
                
% Set up proposal counters
AcceptedMutation  = 0;
AttemptedMutation = 0;

% Set up parameter history variable
ParaHistory         = zeros(NumOfPosteriorSamples, NumOfParameters);
true_para_history   = zeros(NumOfPosteriorSamples, NumOfParameters);
LLHistory           = zeros(NumOfPosteriorSamples, NumOfSpecies);

% Set monitor rate for adapting step sizes
MonitorRate = 10;

% Set up converged flag
ContinueIterations = true;
Converged          = false;
    
disp('Initialisation Completed..');

% Initialise iteration number
IterationNum = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate gradient, metric, etc. at initial point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% actual parameters (with correct order of magnitude)
true_parameters = Parameters.*10.^Options.ParameterOrderOfMagnitude;

% calculate gradient and metric on true scale...
[OldLL, OldGradL, OldG] = gradient_metric(n_tot, TimePoints, Y, true_parameters, ParametersToInfer, Options, {'log-likelihood', 'gradient', 'metric'});
% then scale, to do sampling on rescaled parameters
OldGradL = OldGradL.*(10.^(Options.ParameterOrderOfMagnitude(ParametersToInfer)'));
OldG = OldG.*(10.^Options.ParameterOrderOfMagnitude(ParametersToInfer) * 10.^Options.ParameterOrderOfMagnitude(ParametersToInfer)');

OldGInv = inv(OldG + eye(NumOfParameters)*1e-6);
OldMean = Parameters(ParametersToInfer) + OldGInv*OldGradL'*StepSize/2;

OldLogPrior = 0;
if Parameters(end) >0 % Prior zero if sd(error)<=0
    for i = MathParToInfer
        OldLogPrior = OldLogPrior + log(ModelParameterPrior(i, true_parameters(i)));
    end
else
    OldLogPrior = -Inf;
end
%OldLogPrior

ParaHistory(1,:)         = Parameters(ParametersToInfer);
true_para_history(1,:)   = true_parameters(ParametersToInfer);
LLHistory(1,:)           = OldLL;

disp('Initialisation Completed..');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
while ContinueIterations
    
    % Increment iteration number
    IterationNum = IterationNum + 1
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate parameter values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sample new point in parameter space, and accept with metropolis ratio
        
    % make a proposal
    AttemptedMutation  = AttemptedMutation + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parameters of mathematical model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    NewParas    = Parameters;
    NewParas(ParametersToInfer) = OldMean + chol(StepSize*OldGInv)*randn(NumOfParameters,1);

    new_true_parameters = NewParas.*(10.^Options.ParameterOrderOfMagnitude)

   
    % calculate stuff with proposed parameters

%     Parameters(1:end-2)
%     NewParas(1:end-2)
%     Parameters(1:end-2)./NewParas(1:end-2)

try
    
    [NewLL, NewGradL, NewG] = gradient_metric(...
                    n_tot, ...
                    TimePoints, ...
                    Y, ...
                    new_true_parameters, ...
                    ParametersToInfer, ...
                    Options, ...
                    {'log-likelihood', 'gradient', 'metric'});
                
    % scale
    NewGradL = NewGradL.*10.^(Options.ParameterOrderOfMagnitude(ParametersToInfer)');
    NewG     = NewG.*(10.^Options.ParameterOrderOfMagnitude(ParametersToInfer) * 10.^Options.ParameterOrderOfMagnitude(ParametersToInfer)');

    NewGInv = inv(NewG + eye(NumOfParameters)*1e-6);
    NewMean = NewParas(ParametersToInfer) + NewGInv*NewGradL'*StepSize/2;

    NewLogPrior = 0;

    if Parameters(end) > 0 % Prior zero if sd(error)<=0 
        for i = MathParToInfer % priors for noise will cancel
            NewLogPrior = NewLogPrior + log(ModelParameterPrior(i, new_true_parameters(i)));
        end
    else
        NewLogPrior = -Inf;
    end
%    NewLogPrior

    % and p(old|new) and p(new|old)

    LogProbOldGivenNew = -sum(log(diag(chol(StepSize*NewGInv)))) - 0.5*(NewMean-Parameters(ParametersToInfer))'*(NewG/StepSize)*(NewMean-Parameters(ParametersToInfer));

    LogProbNewGivenOld = -sum(log(diag(chol(StepSize*OldGInv)))) - 0.5*(OldMean-NewParas(ParametersToInfer))'*(OldG/StepSize)*(OldMean-NewParas(ParametersToInfer));


    catch
        NewLL = -1e300;
        NewLogPrior = -1e300;
        LogProbOldGivenNew = 0;
        LogProbNewGivenOld = 0;
        NewMean = 0;
        NewG = OldG;
        NewGInv = OldGInv;
        disp('bad proposal')
    end
    
    % accept or reject

%         NewLL  - OldLL
% %         NewLogPrior
% %         OldLogPrior
%         NewLogPrior - OldLogPrior
%         LogProbOldGivenNew - LogProbNewGivenOld

% priors are both multiplied by scalings -> cancel out
% NewLL - OldLL
% NewLogPrior - OldLogPrior
% LogProbOldGivenNew - LogProbNewGivenOld


     Ratio = NewLL + NewLogPrior + LogProbOldGivenNew ...
                - OldLL - OldLogPrior - LogProbNewGivenOld;
%     exp(Ratio)
            
        if Ratio > 0 || (log(rand) < Ratio)
            % Accept proposal
            % Update variables
            %'accept'
            Parameters                 = NewParas;
            true_parameters            = new_true_parameters;
            AcceptedMutation           = AcceptedMutation + 1;
            OldLL                      = NewLL;
            OldLogPrior                = NewLogPrior;
            OldG                       = NewG;
            OldGInv                    = NewGInv;
            OldMean                    = NewMean;
        end

% %%%%%%%%%%%%%%%%%%%%%%%
% % error model
% %%%%%%%%%%%%%%%%%%%%%%%
% 
%     if ~isempty(StatParToInfer)
%     
%         
%         
%     end
   
    
    %%%%%%%%%%%%%%%%%%%
    % Save parameters %
    %%%%%%%%%%%%%%%%%%%
%     if Converged
%         ParaHistory(IterationNum-ConvergenceIterationNum, :) = Parameters(Options.ParametersToInfer);
%         LLHistory(IterationNum-ConvergenceIterationNum, :)   = CurrentLL;
%     end
%    if Converged
        ParaHistory(IterationNum, :) = Parameters(ParametersToInfer);
        true_para_history(IterationNum, :) = true_parameters(ParametersToInfer);
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
                    
                    disp([num2str(100*AcceptedMutation/AttemptedMutation) '% mutation acceptance']);
                            
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
        
        % Adjust parameter proposal widths
        if mod(IterationNum, 1000) == 0
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['Iteration ' num2str(IterationNum)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(' ')
        end
        
    end
end   


PosteriorTime = toc;

CurTime = fix(clock);
RandTime = ceil(rand*10000);

% Save posterior
FileName = ['ODE_mMALASimp_' EquationName '_' num2str(D) 'DPS_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '_' num2str(RandTime)];
save(['./Results/' FileName], 'ParaHistory', 'true_para_history', 'LLHistory', 'StepSize', 'BurnInTime', 'PosteriorTime', 'Y', 'TimePoints');



end
