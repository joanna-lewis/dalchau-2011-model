function PP = ModelParameterPrior(ParaNum, Value)

%%%%%%%%%
% Lognormal priors %
%%%%%%%%%

% define means and standard deviations of prior distributions
means = log([8.764e-4, 5.658e-6, 4.177e-7,... %(u)
                    2.093e4, 1.759e4, 1.064e4,...    %(g)
                    1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
                    8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
                    0.1142, 150.5, 1.663e-9,...            %(e, g_M, b_T)
                    7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
                    1505, ...
                    1.4920]);
sds = 0.1*abs(means);

if Value == 'random'
    % for the moment, initialise chains at true values
    PP = exp(means(ParaNum));
    % Produce random value from the prior
    %PP = random('LogNormal',means(ParaNum), sds(ParaNum));
else
    if (Value < 0)
        PP = 0;
    else
        % Calculate probability of value from the prior
        PP = lognpdf(Value, means(ParaNum), sds(ParaNum));
    end
end

        
       
end