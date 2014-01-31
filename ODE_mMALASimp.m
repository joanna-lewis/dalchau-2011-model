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
n_o   = length(Options.ObservedViralSpecies);
n_u   = length(Options.UnobservedViralSpecies);
n_s   = length(Options.SelfSpecies);
n_tot = n_o + n_u + n_s; % total number of peptide species


MathParToInfer        = Options.MathParToInfer;
StatParToInfer        = Options.StatParToInfer;
ParametersToInfer     = [MathParToInfer, 2*n_tot+13 + StatParToInfer];

NumOfParameters       = numel(ParametersToInfer);

StepSize = 10^(0);

SpeciesObserved       = Options.ObservedViralSpecies;

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
                    1, ...      % scale factor
                    1.9768]);         % sd(error) for epitope data
                
% Set up proposal counters
AcceptedMutation  = 0;
AttemptedMutation = 0;

% Set up parameter history variable
ParaHistory         = zeros(NumOfPosteriorSamples, NumOfParameters);
LLHistory           = zeros(NumOfPosteriorSamples, NumOfSpecies);

% Set monitor rate for adapting step sizes
MonitorRate = 100;

% Set up converged flag
ContinueIterations = true;
Converged          = false;
    
disp('Initialisation Completed..');

% Initialise iteration number
IterationNum = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate gradient, metric, etc. at initial point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate stuff with old parameters - could speed this up by
        % remembering from previous iteration
        [OldXEstimates, OldXSens] = ode_model_sol_vargen(Parameters(1:end-2), ...
                        n_tot, ...
                        1, ...
                        TimePoints, ...
                        Options.InitialValues, ...
                        Options.gdata, ...
                        ParametersToInfer(ParametersToInfer<=2*n_tot+13));
        
        OldLL       = LogNormPDF(Y(2:end), OldXEstimates(3*n_tot+1,:), Parameters(end)^2);
        
        OldLogPrior = 0;
        for i = ParametersToInfer
            OldLogPrior = OldLogPrior + log(ModelParameterPrior(i, Parameters(i)));
        end
        
        % wrt model parameters
        OldGradL = ((Y(2:end) - OldXEstimates(3*n_tot+1,:)) ...
                     *permute(OldXSens(3*n_tot+1,:,:),[2,3,1])) ...
                     /Parameters(end)^2;
        
        if length(StatParToInfer)>0
        % wrt noise
        OldGradL = cat(2, OldGradL, - ...
            (length(Y)-1)/Parameters(end) + sum((Y(2:end) - OldXEstimates(3*n_tot+1,:)).^2))/Parameters(end)^3;
        end
        
        % add in prior
        for i = 1:NumOfParameters
            OldGradL(i) = OldGradL(i) + ModelParameterLogPriorDerivative(i, Parameters(i));
        end
        
        OldG = zeros(NumOfParameters);
        % mathematical parameters
        for i = 1:length(MathParToInfer)
            for j = i:length(MathParToInfer)
            
                OldG(i,j) = (1/Parameters(end)^2)*...
                    (OldXSens(3*n_tot+1,:,i)*OldXSens(3*n_tot+1,:,j)');
                
                OldG(j,i) = OldG(i,j);
                
            end
        end
        if length(StatParToInfer)>0
        % error in epitope expression
        for i = length(MathParToInfer)+1:NumOfParameters
            for j = i:NumOfParameters
            
                OldG(i,j) = 2*length(Y(2:end))/Parameters(end)^2   ;             
                OldG(j,i) = OldG(i,j);
                
            end
        end
        end

% prior=zeros(NumOfParameters);
%         % add prior to Fisher Information
%         for i = 1:NumOfParameters
%             prior(i,i) = - ModelParameterLogPriorDerivative2(ParametersToInfer(i), Parameters(ParametersToInfer(i))); 
%         end
        % add prior to Fisher Information
        for i = 1:NumOfParameters
            ModelParameterLogPriorDerivative2(ParametersToInfer(i), Parameters(ParametersToInfer(i)))
            OldG(i,i) = OldG(i,i) - ModelParameterLogPriorDerivative2(ParametersToInfer(i), Parameters(ParametersToInfer(i))); 
        end
        OldG
        OldGInv = inv(OldG + eye(NumOfParameters)*1e-6);

        OldMean = Parameters(ParametersToInfer) + OldGInv*OldGradL'*StepSize/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
while ContinueIterations
    
    % Increment iteration number
    IterationNum = IterationNum + 1;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate parameter values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sample new point in parameter space, and accept with metropolis ratio
        
    % make a proposal
    AttemptedMutation  = AttemptedMutation + 1;

    NewParas    = Parameters;
    NewParas(ParametersToInfer) = OldMean + chol(StepSize*OldGInv)*randn(NumOfParameters,1);
    NewParas(ParametersToInfer) = OldMean + chol(StepSize*OldGInv)*randn(NumOfParameters,1);
    
   % NewParas(a) - Parameters(a)

    % calculate stuff with proposed parameters
try
%     Parameters(1:end-2)
%     NewParas(1:end-2)
%     Parameters(1:end-2)./NewParas(1:end-2)
    
        [NewXEstimates, NewXSens] = ode_model_sol_vargen(NewParas(1:end-2), ...
                    n_tot, ...
                    1, ...
                    TimePoints, ...
                    Options.InitialValues, ...
                    Options.gdata, ...
                    ParametersToInfer(ParametersToInfer<=2*n_tot+13));
%                size(NewXSens)

        NewLL = LogNormPDF(Y(2:end), NewXEstimates(3*n_tot+1,:), NewParas(end)^2);

                
        % wrt parameters
        NewGradL = ((Y(2:end) - NewXEstimates(3*n_tot+1,:)) ...
                     * permute(NewXSens(3*n_tot+1,:,:),[2,3,1])) ...
                     /NewParas(end)^2;
        
        if length(StatParToInfer)>0
        % wrt noise
        NewGradL = cat(2, NewGradL, - ...
            ((length(Y)-1)/NewParas(end) + sum((Y(2:end) - NewXEstimates(3*n_tot+1,:)).^2))/NewParas(end)^3);
        end
        
         % add in prior
        for i = 1:NumOfParameters
            NewGradL(i) = NewGradL(i) + ModelParameterLogPriorDerivative(i, NewParas(i));
        end
        
       
        NewG = zeros(NumOfParameters);
        
        % mathematical parameters
        for i = 1:length(MathParToInfer)
            for j = i:length(MathParToInfer)
  
                NewG(i,j) = (1/Parameters(end)^2)*...
                    (NewXSens(3*n_tot+1,:,i)*NewXSens(3*n_tot+1,:,j)');

                NewG(j,i) = NewG(i,j);

            end
        end
        
        if length(StatParToInfer)>0
        % error in epitope expression
        for i = length(MathParToInfer)+1:NumOfParameters
            for j = i:NumOfParameters
            
                NewG(i,j) = 2*length(Y(2:end))/Parameters(end)^2  ;              
                NewG(j,i) = NewG(i,j);
                
            end
        end
        end
            
        % add prior to Fisher Information
        for i = 1:NumOfParameters
            NewG(i,i) = NewG(i,i) - ModelParameterLogPriorDerivative2(ParametersToInfer(i), NewParas(ParametersToInfer(i)));
        end

    NewGInv = inv(NewG + eye(NumOfParameters)*1e-6);
    
    NewMean = NewParas(ParametersToInfer) + NewGInv*NewGradL'*StepSize/2;

    catch
        NewLL = -1e300;
        NewMean = 0;
        NewG = OldG;
        NewGInv = OldGInv;
        disp('bad proposal')
    end

    NewLogPrior = 0;
    for i = ParametersToInfer
        NewLogPrior = NewLogPrior + log(ModelParameterPrior(i, NewParas(i)));
    end

    % and p(old|new) and p(new|old)

    LogProbNewGivenOld = -sum(log(diag(chol(StepSize*OldGInv)))) - 0.5*(OldMean-NewParas(ParametersToInfer))'*(OldG/StepSize)*(OldMean-NewParas(ParametersToInfer));

    LogProbOldGivenNew = -sum(log(diag(chol(StepSize*NewGInv)))) - 0.5*(NewMean-Parameters(ParametersToInfer))'*(NewG/StepSize)*(NewMean-Parameters(ParametersToInfer));


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
            AcceptedMutation           = AcceptedMutation + 1;
            OldLL                      = NewLL;
            OldLogPrior                = NewLogPrior;
            OldG                       = NewG;
            OldGInv                    = NewGInv;
            OldMean                    = NewMean;
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
save(['./Results/' FileName], 'ParaHistory', 'LLHistory', 'StepSize', 'BurnInTime', 'PosteriorTime', 'Y', 'TimePoints');



end
