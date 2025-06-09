function [EE, psi] = EE_CU_BF_pass(B, Gamma, F, psi, Sigma, P_n, eta, P_fix, R_min, B_W)

    [~, N_R, L, K] = size(B);

    G = zeros(L, K);
    for k = 1:K      
        G(:, k) = diag( (Sigma(:, k)).^(0.5) ) * Gamma(:, k);
    end

    D = zeros(N_R, K, K);
    for k = 1:K
        for m = 1:K
            for l = 1:L
                D(:, m, k) = D(:, m, k) + conj( G(l, k) ) * B(:, :, l, k).' * conj( F(:, m) );
            end
        end
    end

    R = zeros(K,1);
    mu = zeros(K,1);
    xi = P_n * ones(K,1);
    for k = 1:K
        Index_except_k = 1:K;
        Index_except_k(Index_except_k == k) = [];
        for m = Index_except_k        
            xi(k) = xi(k) + abs( D(:,m,k)' * psi ) ^ 2;
        end
        mu(k) = abs( D(:,k,k)' * psi ) ^ 2 / xi(k);
        R(k) = B_W * log2( 1 + mu(k) );
    end
    EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

    if all( R >= R_min )
        EE_save = EE;
        psi_save = psi;
        Sign_cons = 1;
    else                   
        EE_save = 0;
        psi_save = zeros(N_R, 1);
        Sign_cons = 0;
    end
   
    Iter_SCA = 10;  
    iter_SCA = 1;
    eta_SCA = 10 ^ (-2);
    EE = EE_save;
    EE_last = 10^(-10);
    C = 10^(0);    

    while (EE / EE_last - 1 >= eta_SCA) && (iter_SCA <= Iter_SCA) && (Sign_cons == 1)

        EE_last = EE;
        psi_last = psi;
        mu_last = mu;
        xi_last = xi;

        cvx_clear
        cvx_begin quiet
            variable psi(N_R) complex;
            variable mu(K);
            variable xi(K);
            expressions Goal_max Cons_1_l(K) Cons_1_r(K) Cons_2_l(K) Cons_4_l(K) complex;
    
            for k = 1:K
                Goal_max = Goal_max + 1 / (2 * log(2)) * mu(k) / (1 + mu_last(k));
                Cons_1_l(k) = sqrt( mu_last(k) * xi_last(k) ) ...
                              + sqrt( mu_last(k) / xi_last(k) ) * ( xi(k) - xi_last(k) ) / 2 ...
                              + sqrt( xi_last(k) / mu_last(k) ) * ( mu(k) - mu_last(k) ) / 2;
                Cons_1_r(k) = D(:,k,k)' * psi;
                Index_except_k = 1:K;
                Index_except_k(Index_except_k == k) = [];
                for m = Index_except_k        
                    Cons_2_l(k) = Cons_2_l(k) + square_pos( norm( D(:,m,k)' * psi ) );
                end
                Cons_2_l(k) = Cons_2_l(k) + P_n;
            end    

            for n_R = 1:N_R
                Goal_max = Goal_max + 2 * C * conj( psi_last(n_R) ) * psi(n_R);
                Cons_4_l(n_R) = conj(psi(n_R)) * psi(n_R);
            end

            maximize ( real( Goal_max ) )
            subject to
                real( Cons_1_l ) <= real( Cons_1_r );
                real( Cons_2_l ) <= xi;           
                mu >= 2.^(R_min / B_W) - ones(K,1);
                real( Cons_4_l ) <= ones(N_R,1);
                mu >= 0;
                xi >= 0;
    
        cvx_end  

        psi_test = psi ./ abs(psi);
        R = zeros(K,1);
        IN = P_n * ones(K,1);
        for k = 1:K
            Index_except_k = 1:K;
            Index_except_k(Index_except_k == k) = [];
            for m = Index_except_k        
                IN(k) = IN(k) + abs( D(:,m,k)' * psi_test ) ^ 2;
            end
            R(k) = B_W * log2( 1 + abs( D(:,k,k)' * psi_test ) ^ 2 / IN(k) );
        end
        EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

        if ( all( R >= R_min ) ) && ( EE > EE_save )
            EE_save = EE;
            psi_save = psi_test;
            Sign_cons = 1;
        else
            Sign_cons = 0;
        end

        iter_SCA = iter_SCA + 1;

    end

    EE = EE_save;
    psi = psi_save;

end
