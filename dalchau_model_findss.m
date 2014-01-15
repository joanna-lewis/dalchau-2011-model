function[init]=dalchau_model_findss(parms, n_o, n_u, n_s) %u, u_T, q, g, b, c, d_p, v, e, g_M, b_T, d_M, d_T, d_Me)
    
    % n_o is no. observed viral peptides
    % n_u is no. unobserved viral peptides
    % n_s is no. unobserved self peptides
    
    % state vector is: (P_is, MP_is, TMP_is, MeP_is, M, T, TM, Me)
    % parms is vector of parameters: (u, g, u_T, q, b, c, d_p, v, e, g_M, b_T, d_M, d_T, d_Me)
    % t_label is a two-element vector: [time labelling began, time
    % labelling ended]
    % values from the paper (for B2705) were 
    % parms=transpose([8.764e-4, 5.658e-6, 4.177e-7,... %(u)
    %        2.093e4, 1.759e4, 1.064e4,...    %(g)
    %        1.185e-6, 2.104e4, 1.945e-9,...   %(u_T, q, b_B2705)
    %        8.303e-8, 0.13, 936.3,...              %(c, d_p, v)
    %        0.1142, 150.5, 1.663e-9,...            %(e, g_M, b_T)
    %        7.989e-5, 1.726e-3, 9.329e-5, ...      % (d_M, d_T, d_Me)
    %        1505]) % (g_T)
   
    % only self peptide species present initially
    n = n_o+n_u+n_s;

    u    = parms(n_o+n_u+1:n);
    g    = parms(n+n_o+n_u+1:2*n);
    u_T  = parms(2*n+1);
    q    = parms(2*n+2);
    b    = parms(2*n+3);
    c    = parms(2*n+4);
    d_p  = parms(2*n+5);
    v    = parms(2*n+6);
    e    = parms(2*n+7);
    g_M  = parms(2*n+8);
    b_T  = parms(2*n+9);
    d_M  = parms(2*n+10);
    d_T  = parms(2*n+11);
    d_Me = parms(2*n+12);
    g_T  = parms(2*n+13);
    
    
        function[y_dot] = model(y)
            % ode model from paper
            % y is current state vector:
            % n is the number of peptides being considered
            % y_dot is dy/dt at this y
       
            P_is   = y(1:n_s);
            MP_is  = y(n_s+1:2*n_s) ;
            TMP_is = y(2*n_s+1:3*n_s) ;
            MeP_is = y(3*n_s+1:4*n_s); 
            M      = y(4*n_s+1) ;
            T      = y(4*n_s+2); 
            TM     = y(4*n_s+3);
            Me     = y(4*n_s+4);
            
            
            y_dot = zeros(size(y));
            
            y_dot(1:n_s)       = (u.*MP_is + q*(u.*TMP_is) - (b*M+c*TM+d_p).*P_is)+ g; % Eqn 8: d(P_is)/dt
            
            y_dot(n_s+1:2*n_s)   = b*M*P_is + u_T*v*TMP_is - (u+e).*MP_is; % Eqn 5 d(MP_is)/dt (unlabelled)
            
            y_dot(2*n_s+1:3*n_s) = c*TM*P_is - (u*q+u_T*v).*TMP_is; % Eqn 7 d(TMP_is)/dt (unlabelled)
            
            y_dot(3*n_s+1:4*n_s) = e*MP_is - u.*MeP_is; % Eqn 9 d(MeP_is)/dt (unlabelled)
                        
            y_dot(4*n_s+1)     = dot(u,MP_is) + u_T*TM + g_M - (b*sum(P_is)+b_T*T +d_M)*M; % Eqn 3 dM/dt (unlabelled)
            
            y_dot(4*n_s+2)     = u_T*TM + g_T +u_T*v*sum(TMP_is) - (b_T*M + d_T)*T; % Eqn 4 dT/dt

            y_dot(4*n_s+3)     = b_T*M*T + q*dot(u,TMP_is) - (u_T+c*sum(P_is))*TM; % Eqn 6 d(TM)/dt (unlabelled)
            
            y_dot(4*n_s+4)     = dot(u,MeP_is) - d_Me*Me; % Eqn 10 d(Me)/dt (unlabelled)
          
        end
        
 ss = fsolve(@model,cat(1,...
        6*10^2*ones(n_s,1), ...     % self P_is
        2*10^2*ones(n_s,1), ...     % self MP_is
        2*10^4*ones(n_s,1), ...     % self TMP_is
        2*10^6*ones(n_s,1), ...     % self MeP_is ...
        [10^3; 10^4; 10^8; 10^1]),...
        optimset('MaxFunEvals', 10000000,'MaxIter',10000000, 'Display', 'off')...
    );
    
init = [zeros(n_o+n_u, 1); ss(1:n_s); ...       % self P_is
        zeros(n_o+n_u, 1); ss(n_s+1:2*n_s);...  % self MP_is
        zeros(n_o+n_u, 1); ss(2*n_s+1:3*n_s);...% self TMP_is
        zeros(n_o+n_u, 1); ss(3*n_s+1:4*n_s);...% self MeP_is
        ss(4*n_s+1:end)];
end
 
           