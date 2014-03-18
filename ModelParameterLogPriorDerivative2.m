function PP = ModelParameterLogPriorDerivative2(ParaNum, Value)

% Return the partial derivative of the log(prior) w.r.t. the parameter

% Lognormal Priors

% define means and standard deviations of prior distributions
modes = [8.764e-4, 5.658e-6, 4.177e-7,... %(u)
                    2.093e4, 1.759e4, 1.064e4,...    %(g)
                    1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
                    8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
                    0.1142, 150.5, 1.663e-9,...            %(e, g_M, b_T)
                    7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
                    1505, ...
                    1, ...      % proportionality for tryptic data
                    1.9768];         % sd(error) for epitope data

%sds = 0.1*abs(modes);
sds = abs(modes).*(log10(modes)+9);
means = log(modes) + sds.^2;


if (Value <= 0)
    PP = 0;
else
    %PP = 0;
    PP = (1 - (1 - log(Value) + means(ParaNum))/sds(ParaNum).^2)./Value.^2;
end


end
