% plot prior

mpp = [];
for i=-8:0.1:-2
    mpp =  cat(1,mpp,ModelParameterPrior(1, 10^i));
end

figure
semilogx(10.^(-8:0.1:-2)', mpp)
[n, xout] = hist(ParaHistory*10^-7,50);
binwidth = xout(2) - xout(1);

hold on
bar(xout, n/(binwidth*length(ParaHistory)))

fitparms = lognfit(ParaHistory*10^-7);
post = [];
for i=-8:0.1:-2
    post =  cat(1,post,lognpdf(10^i, fitparms(1), fitparms(2)));
end
plot(10.^(-8:0.1:-2)', post, 'r')

% plot simulation with initial parameters
Parameters = transpose([8.764, 5.658, 4.177,... %(u)
                    2.093, 1.759, 1.064,...    %(g)
                    1.185, 2.104, 1.945,...   %(u_T, q, b_B2705)
                    8.303, 1.3, 9.363,...              %(c, d_p, v)
                    1.142, 1.505, 1.663,...            %(e, g_M, b_T)
                    7.989, 1.726, 9.329, ...      % (d_M, d_T, d_Me)
                    1.505]).*10.^...
             transpose([-7, -9, -10, ...
                    2, 2, 2, ...
                    -9, 4, -13, ...
                    -12, -4, 2, ...
                    -4, 0, -13, ...
                    -8, -6, -8, ...
                    1]);
                
TimePoints    = 3600*([0, 0.5, 1:9])';
NoisyData = csvread('B8_data_260314.csv'); % same as simulate_data_111213_2.csv, but tryptic peptide levels set to zero
gdata = [TimePoints'; NoisyData(1:3,:)];


InitialValues = dalchau_model_findss(...
    Parameters, 1, 1, 1);
                
solution = ode_model_sol_vargen(Parameters, ...
    n_tot, ...
    1, ...
    TimePoints, ...
    InitialValues, ...
    gdata, ...
    [1]);

figure
hold on
plot(TimePoints', NoisyData(4,:), 'o k')
plot(TimePoints(2:end)', solution(10,:))

% now add simulation with posterior parameters
map = exp(fitparms(1) - fitparms(2)^2);
Parameters(1) = map

InitialValues = dalchau_model_findss(...
    Parameters, 1, 1, 1);
                
solution = ode_model_sol_vargen(Parameters, ...
    n_tot, ...
    1, ...
    TimePoints, ...
    InitialValues, ...
    gdata, ...
    [1]);
plot(TimePoints(2:end)', solution(10,:),'r')