function[sol, solS]=ode_model_sol(parms, n, solveat, initial, sens_required)

    % state vector is: (P_is, MP_is, TMP_is, MeP_is, M, T, TM, Me)
    % parms is vector of parameters: (u, g, u_T, q, b, c, d_p, v, e, g_M, b_T, d_M, d_T, d_Me)
    % solveat is the timepoints in the experiment (including the initial tp)   
    % initial is the state of the system at t=0 
    % sens_required is a vector giving the indices of the parameters for
    % which sensitvity analysis is required
    
    % values from the paper (for B2705) were 
    % parms=transpose([8.764e-4, 5.658e-6, 4.177e-7,... %(u)
    %        2.093e4, 1.759e4, 1.064e4,...    %(g)
    %        1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
    %        8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
    %        0.1142, 150.5, 1.663e-9,...            %(e, g_M, b_T)
    %        7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
    %        1505]) % (g_T)
   
    % put parameters into data structure for cvodes
    data.p = parms;
    data.n = n;
    options = CVodeSetOptions(  'UserData', data,...
                                'RelTol',1.e-8,...
                                'AbsTol',1.e-6,...
                                'LinearSolver','Dense');
    %%%%%%%%%%%%%%%%%%%%%%%
    % Initialise cvodes
    %%%%%%%%%%%%%%%%%%%%%%%
    
    t0 = solveat(1);   % first timepoint
    tps = solveat(2:end);    % timepoints required
    y0 = initial;     % initial value

    CVodeInit(@model, 'BDF', 'Newton', t0, y0, options);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Initialise cvodes
    %%%%%%%%%%%%%%%%%%%%%%%
    
    Ns = numel(sens_required); % just look at one parameter for the moment
    yS0 = zeros(4*n+4,Ns);  %Initial conditions for sensitivity variables. YS0 must 
                        % be a matrix with N rows and Ns columns, where N is the 
                        % problem dimension and Ns the number of sensitivity systems.
    FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                 'ErrControl', true,...
                                 'ParamField', 'p',...
                                 'ParamScales', data.p(sens_required));
                             
    CVodeSensInit(Ns, [], yS0, FSAoptions);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % problem solution
    %%%%%%%%%%%%%%%%%%%%%%%
    
    [status, t, sol, solS] = CVode(solveat(2:end),'Normal');
    solS                   = permute(solS, [1 3 2]);

    CVodeFree;

return

function[y_dot, flag, new_data] = model(t, y, data)

    parms = data.p;
    n = data.n;

    u=parms(1:n);
    g=parms(n+1:2*n);
    u_T=parms(2*n+1);
    q=parms(2*n+2);
    b=parms(2*n+3);
    c=parms(2*n+4);
    d_p=parms(2*n+5);
    v=parms(2*n+6);
    e=parms(2*n+7);
    g_M=parms(2*n+8);
    b_T=parms(2*n+9);
    d_M=parms(2*n+10);
    d_T=parms(2*n+11);
    d_Me=parms(2*n+12);
    g_T=parms(2*n+13);
    
    
            % ode model from paper
            % y is current state vector:
            % U for unlabelled, L for labelled
            % n is the number of peptides being considered
            % y_dot is dy/dt at this y
            P_is=y(1:n);
            MP_is_U=y(n+1:2*n) ;
            TMP_is_U=y(2*n+1:3*n);
            MeP_is_U=y(3*n+1:4*n); 
            M_U=y(4*n+1); 
            T=y(4*n+2); 
            TM_U=y(4*n+3);
            Me_U=y(4*n+4);
            
            
 %           y_dot = zeros(size(y));
 
            y_dot(1:n) = (u.*MP_is_U + q*(u.*TMP_is_U) - (b*M_U+c*TM_U+d_p)*P_is)+ g; % Eqn 8: d(P_is)/dt
            
            y_dot(n+1:2*n) = b*M_U*P_is + u_T*v*TMP_is_U - (u+e).*MP_is_U; % Eqn 5 d(MP_is)/dt (unlabelled)
            
            y_dot(2*n+1:3*n) = c*TM_U*P_is - (u*q+u_T*v).*TMP_is_U; % Eqn 7 d(TMP_is)/dt (unlabelled)
            
            y_dot(3*n+1:4*n) = e*MP_is_U - u.*MeP_is_U; % Eqn 9 d(MeP_is)/dt (unlabelled)
                        
            y_dot(4*n+1) = dot(u,MP_is_U) + u_T*TM_U + g_M - (b*sum(P_is)+b_T*T +d_M)*M_U; % Eqn 3 dM/dt (unlabelled)
            
            y_dot(4*n+2) = u_T*TM_U + g_T +u_T*v*sum(TMP_is_U) - (b_T*M_U + d_T)*T; % Eqn 4 dT/dt

            y_dot(4*n+3) = b_T*M_U*T + q*dot(u,TMP_is_U) - (u_T+c*sum(P_is))*TM_U; % Eqn 6 d(TM)/dt (unlabelled)
            
            y_dot(4*n+4) = dot(u,MeP_is_U) - d_Me*Me_U; % Eqn 10 d(Me)/dt (unlabelled)
                        
            flag=0;
            new_data = [];

return