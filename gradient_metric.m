% function to calculate gradient and/or metric, from true parameters

function[LL, gradient, metric] = gradient_metric(n_tot, ...
                                             TimePoints, ...
                                             Y, ...
                                             true_parameters, ...
                                             ParametersToInfer, ...
                                             Options, ...
                                             req)


% req indicates which of log-likelihood ('log-likelihood'), gradient ('gradient'), metric ('metric')

LL = [];
gradient = [];
metric = [];

NumOfParameters = length(ParametersToInfer);
MathParToInfer = Options.MathParToInfer;

[InitialValues, InitialJacobian] = dalchau_model_findss(...
                    true_parameters(1:(2*n_tot+13)), ...
                    length(Options.ObservedViralSpecies), ...
                    length(Options.UnobservedViralSpecies), ...
                    length(Options.SelfSpecies), ...
                    ParametersToInfer(ParametersToInfer<=2*n_tot+13));

[XEstimates, XSens] = ode_model_sol_vargen(true_parameters(1:2*n_tot+13), ...
                    n_tot, ...
                    1, ...
                    TimePoints, ...
                    InitialValues, ...
                    Options.gdata, ...
                    ParametersToInfer(ParametersToInfer<=2*n_tot+13));


                
% [XEstimates, XSens] = ode_model_sol_vargen(true_parameters(1:2*n_tot+13), ...
%                     n_tot, ...
%                     1, ...
%                     TimePoints, ...
%                     InitialValues, ...
%                     Options.gdata, ...
%                     ParametersToInfer(ParametersToInfer<=2*n_tot+13))
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember('log-likelihood', req) || ismember('gradient', req)
    
    LL = LogNormPDF(Y(2:end), XEstimates(2*n_tot+1,:), true_parameters(end)^2);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember('gradient', req)
           
        
%     LogPrior = 0;
%     for i = ParametersToInfer
%         LogPrior = LogPrior + log(ModelParameterPrior(i, true_parameters(i)));
%     end
    
    % wrt model parameters
    
    gradient = ((Y(2:end) - XEstimates(2*n_tot+1,:)) ...
                 *permute(XSens(2*n_tot+1,:,:),[2,3,1])) ...
                 /true_parameters(end)^2;

    if ismember(length(true_parameters), ParametersToInfer)
    % wrt noise
    gradient = cat(2, gradient, - ...
        (length(Y)-1)/true_parameters(end) + sum((Y(2:end) - XEstimates(2*n_tot+1,:)).^2)/true_parameters(end)^3);
    end
    % will have to put wrt scaling here, in due course

    % add in prior
    for i = 1:NumOfParameters
        gradient(i) = gradient(i) + ...
            ModelParameterLogPriorDerivative(ParametersToInfer(i), true_parameters(ParametersToInfer(i)));
    end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if ismember('metric', req)
       
        metric = zeros(NumOfParameters);
        
        % mathematical parameters
        for i = 1:length(MathParToInfer)
            for j = i:length(MathParToInfer)
                
                    metric(i,j) = (1/true_parameters(end)^2)*...
                        (XSens(2*n_tot+1,:,i)*XSens(2*n_tot+1,:,j)');
                
                metric(j,i) = metric(i,j);
                
            end
        end
        if ismember(length(true_parameters), ParametersToInfer)
        % error in epitope expression
            dim = size(metric,1);
            metric(dim, dim) = 2*length(Y(2:end))/true_parameters(end)^2;                
        end
        
%
% prior=zeros(NumOfParameters);
%         % add prior to Fisher Information
%         for i = 1:NumOfParameters
%             prior(i,i) = - ModelParameterLogPriorDerivative2(ParametersToInfer(i), Parameters(ParametersToInfer(i))); 
%         end
        %add prior to Fisher Information
                
        for i = 1:NumOfParameters
            metric(i,i) = metric(i,i) - ModelParameterLogPriorDerivative2(ParametersToInfer(i), true_parameters(ParametersToInfer(i))); 
        end
                        
    end
    
end