function [EE, F, psi] = EE_U_BF_hyb(B, Gamma, Gamma_s, Sigma, L_set, P_n, eta, P_fix, P_t, R_min, B_W, psi)   

    [N_B, ~, L, K] = size(B);  

    F = zeros(N_B, K);
    
    Iter_alt = 100;
    iter_alt = 1;
    eta_alt = 5 * 10^(-2);
    EE = 2 * 10^(-10);
    EE_last = 10^(-10);

    while ( EE / EE_last - 1 > eta_alt ) && ( iter_alt <= Iter_alt )

        EE_last = EE;

        [F] = EE_U_BF_act(B, Gamma_s, F, psi, Sigma, L_set, P_n, eta, P_fix, P_t, R_min, B_W);
        [EE, psi] = EE_U_BF_pass(B, Gamma_s, F, psi, Sigma, L_set, P_n, eta, P_fix, R_min, B_W);

        iter_alt = iter_alt + 1;

    end

    G = zeros(L, K);
    for k = 1:K
        G(:, k) = diag( (Sigma(:, k)).^(0.5) ) * Gamma(:, k);
    end  

    A = zeros(N_B, L, K);
    for k = 1:K
        for l = 1:L
            A(:, l, k) = B(:, :, l, k) * conj( psi );
        end
    end     

    H = zeros(N_B, K);
    for k = 1:K
        H(:, k) = A(:, :, k) * conj( G(:, k) );
    end

    R = zeros(K,1);           
    for k = 1:K
        IN = P_n - abs( H(:, k)' * F(:,k) ) ^ 2 ;
        for m = 1:K
            IN = IN + abs( H(:, k)' * F(:,m) ) ^ 2;
        end
        R(k) = B_W * log2( 1 + abs( H(:, k)' * F(:,k) ) ^ 2 / IN );
    end
    R_sum = sum(R);
    EE = R_sum / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );   


end
