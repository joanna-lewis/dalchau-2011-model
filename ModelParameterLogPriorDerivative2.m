function PP = ModelParameterLogPriorDerivative2(ParaNum, Value)

% Return the partial derivative of the log(prior) w.r.t. the parameter

% % Lognormal Priors
% 
% % define means and standard deviations of prior distributions
% modes = [1.3534e-05, 1.3582e-04, 1.3582e-04,... %(u)
%                     0.1359, 1.359e9, 1.359e9,...    %(g)
%                     1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
%                     8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
%                     0.1142, ...
%                     1.359e9, ...
%                     1.663e-9,...            %(e, g_M, b_T)
%                     7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
%                     1505, ...
%                     1.9768];         % sd(error) for epitope data
% 
% % sds = 0.1*abs(modes);
% sds = [1.4142, 0.0739, 0.0739, ...
%     1.4126, 1.414, 1.414, ...
%     sqrt(log(10)/2)*ones(1,7),...
%     10, ...
%     sqrt(log(10)/2)*ones(1,numel(modes)-14)];
% %sds = abs(modes).*(log10(modes)+9);
% %sds = 0.56770116784365*ones(1,numel(modes));
% %sds = 0.4755198067258725*ones(1,numel(modes));
% 
% means = log(modes) + sds.^2;
% 
% 
% if (Value <= 0)
%     PP = 0;
% else
%     if ParaNum < numel(modes)
%     %PP = 0;
%         PP = (1 - (1 - log(Value) + means(ParaNum))/sds(ParaNum).^2)./Value.^2;
%     else
%         PP = 0;  % gradient = 0 for sd(error) (uniform distribution)
%     end
% end


% Gamma priors

% define means and standard deviations of prior distributions
A = 2*ones(1, 20);         
A(1) = 1.2462;
A(4) = 1.2462;
A(end) = 5.3547;

B = 2*ones(1, 20);         
B(1) = 0.0002;
B(4) = 17.0412;
B(end) = 5.4816;


if (Value <= 0)
    PP = 0;
else
    %if ParaNum < numel(A)
    %PP = 0;
        PP = (1 - A(ParaNum))/Value^2 ;
    %else
    %    PP = 0; % gradient = 0 for sd(error) (uniform distribution)
    %end
end

end
