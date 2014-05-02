% plot prior

fitparms = [];

for parnum = 1:size(true_para_history,2)
    
    figure()
    plot(true_para_history(:,parnum))

    xs = floor(log10(min(true_para_history(:,parnum)))):0.01:ceil(log10(max(true_para_history(:,parnum))));
    
    mpp = [];
    for i = xs
        mpp =  cat(1,mpp,ModelParameterPrior(Options.ParametersToInfer(parnum), 10^i));
    end

    figure
    plot(10.^(xs)', mpp)
    [n, xout] = hist(true_para_history(:,parnum),50);
    binwidth = xout(2) - xout(1);

    hold on
    bar(xout, n/(binwidth*size(true_para_history,1)))

    fitparms = cat(1, fitparms, lognfit(true_para_history(:,parnum)));
    post = [];
    for i = xs
        post =  cat(1,post,lognpdf(10^i, fitparms(parnum,1), fitparms(parnum,2)));
    end
    plot(10.^xs', post, 'r')
    
end

% % plot simulation with initial parameters
% Parameters = transpose([1.3582, 1.3582, 1.3582,... %(u)
%                     2*1.359, 2*1.359, 2*1.359,...    %(g)
%                     0, 0, 2.3256/2,...   %(u_T, q, b_B2705)
%                     0, 1.3, 0,...              %(c, d_p, v)
%                     1.142, 1.505, 0,...            %(e, g_M, b_T)
%                     1, 1, 9.329, ...      % (d_M, d_T, d_Me)
%                     0]).*10.^...
%              transpose([-4, -4, -4, ...
%                                     1, 1, 1, ...
%                                     0, 0, -5, ...
%                                     0, -1, 0, ...
%                                     -1, -1, 0, ...       %(e, g_M, b_T)
%                                     -4, 0, -5, ...       % (d_M, d_T, d_Me)
%                                     0]);
                                
% plot simulation with initial parameters
Parameters = transpose([1.3582, 1.3582, 1.3582,... %(u)
                    1.759, 1.759, 1.759,...    %(g)
                    0, 0, 1.663,...   %(u_T, q, b)
                    0, 1.3, 0,...              %(c, d_p, v)
                    1, 1.505, 0,...            %(e, g_M, b_T)
                    1, 1, 1, ...      % (d_M, d_T, d_Me)
                    0]).*10.^...
             transpose([-4, -4, -4, ...
                                    1, 1, 1, ...
                                    0, 0, -6, ...   %(u_T, q, b)
                                    0, -1, 0, ...          %(c, d_p, v)
                                    -1, -1, 0, ...       %(e, g_M, b_T)
                                    -4, 0, -4, ...       % (d_M, d_T, d_Me)
                                    0]);
                                
                                
                
TimePoints    = 3600*([0, 0.5, 1:9])';
NoisyData = csvread('A6_data_310314.csv'); % same as simulate_data_111213_2.csv, but tryptic peptide levels set to zero
gdata = [TimePoints'; NoisyData(1:3,:)];


InitialValues = dalchau_model_findss(...
    Parameters, 1, 1, 1,1);
                
[solution, solS] = ode_model_sol_vargen(Parameters, ...
    n_tot, ...
    1, ...
    TimePoints, ...
    InitialValues, ...
    gdata, ...
    [1]);

simfig=figure();
hold on
plot(TimePoints'/3600, NoisyData(4,:), 'o k')
xlabel('Time since infection / hours', 'FontSize', 20)
ylabel('Epitope copies per cell', 'FontSize', 20)
plot(TimePoints'/3600, [InitialValues(2*n_tot+1), solution(2*n_tot+1,:)])

set(gca, 'fontsize',20)

% now add simulation with posterior parameters
for parnum = 1:size(Options.MathParToInfer,2)
    map = exp(fitparms(parnum,1) - fitparms(parnum,2)^2)
    Parameters(Options.MathParToInfer(parnum)) = map;
end

InitialValues = dalchau_model_findss(...
    Parameters, 1, 1, 1, 1);
                
solution = ode_model_sol_vargen(Parameters, ...
    n_tot, ...
    1, ...
    TimePoints, ...
    InitialValues, ...
    gdata, ...
    [1]);
plot(TimePoints'/3600, [InitialValues(2*n_tot+1), solution(2*n_tot+1,:)],'r')

save2pdf('simplot.pdf')

figure()
plot(true_para_history(:,1), true_para_history(:,2),'.')