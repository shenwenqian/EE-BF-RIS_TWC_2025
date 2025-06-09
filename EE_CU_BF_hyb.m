function [EE, R_sum] = EE_CU_BF_hyb(B, Gamma, Gamma_s, Sigma, L_set, L_RB, L_RU, P_n, eta, P_fix, P_t, R_min, B_W, psi)   

    [N_B, ~, L, K] = size(B);

    F = zeros(N_B, K);

    alpha = (randn(L_RB,1) + 1i * randn(L_RB,1)) / sqrt(2);
    Beta = zeros(L_RU,K);
    Gamma_rand = zeros(L, K);
    for k = 1:K
        Beta(:,k) = (randn(L_RU,1) + 1i * randn(L_RU,1)) / sqrt(2);
        Gamma_rand(:,k) = kron( Beta(:,k), alpha );
        if isempty(L_set) ~= 1
            Gamma_rand(L_set(:,k), k) = Gamma_s(:, k);
        end
    end    

    Iter_alt = 100;
    iter_alt = 1;
    eta_alt = 5 * 10^(-2);
    EE = 2 * 10^(-10);
    EE_last = 10^(-10);

    while ( EE / EE_last - 1 > eta_alt ) && ( iter_alt <= Iter_alt )

        EE_last = EE;

        [F] = EE_CU_BF_act(B, Gamma_rand, F, psi, Sigma, P_n, eta, P_fix, P_t, R_min, B_W);
        [EE, psi] = EE_CU_BF_pass(B, Gamma_rand, F, psi, Sigma, P_n, eta, P_fix, R_min, B_W);

        iter_alt = iter_alt + 1;

    end

    if ( norm(F, 'fro') == 0 ) && ( norm(psi, 'fro') == 0 )

        EE = 0;
        R_sum = 0;

    else

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

end
