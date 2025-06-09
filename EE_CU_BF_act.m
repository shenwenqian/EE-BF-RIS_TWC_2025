function [F] = EE_CU_BF_act(B, Gamma, F, psi, Sigma, P_n, eta, P_fix, P_t, R_min, B_W)   
    
    [N_B, ~, L, K] = size(B);

    A = zeros(N_B, L, K);
    G = zeros(L, K);
    H = zeros(N_B, K);
    for k = 1:K
        for l = 1:L
            A(:, l, k) = B(:, :, l, k) * conj( psi );
        end        
        G(:, k) = diag( (Sigma(:, k)).^(0.5) ) * Gamma(:, k);
        H(:, k) = A(:, :, k) * conj( G(:, k) );
    end

    mu = zeros(K,1);
    xi = P_n * ones(K,1);
    R = zeros(K,1);           
    for k = 1:K
        Index_except_k = 1:K;
        Index_except_k(Index_except_k == k) = [];
        for m = Index_except_k 
            xi(k) = xi(k) + abs( H(:,k)' * F(:,m) ) ^ 2;
        end 
        mu(k) = abs( H(:,k)' * F(:,k) ) ^ 2 / xi(k);
        R(k) = B_W * log2( 1 + mu(k) );
    end 
    EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

    if ( all( R >= R_min ) ) && ( norm( F, 'fro' ) ^ 2 <= P_t * ( 1 + 10^(-3) ) )
        Sign_cons = 1;
        EE_save = EE;
        F_save = F;            
    end

    Iter_SCA = 100;        
    iter_SCA = 1;
    eta_SCA = 10^(-2);    
    EE = EE_save;
    EE_last = 10^(-10);    

    while (EE / EE_last - 1 >= eta_SCA) && (Sign_cons == 1) && (iter_SCA <= Iter_SCA)

        EE_last = EE;
        mu_last = mu;
        xi_last = xi;

        eta_Dink = 10 ^ (-2);
        Iter_Dink = 100;         
        iter_Dink = 1;
        q = sum( log2( ones(K,1) + mu_last ) ) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );
        q_last = -Inf;        

        while (q - q_last >= eta_Dink) && (iter_Dink <= Iter_Dink)

            q_last = q;

            cvx_clear
            cvx_begin quiet
            cvx_solver SeDuMi
                variable F(N_B, K) complex;
                variable mu(K);
                variable xi(K);
                expressions Obj Cons_1_l(K) Cons_1_r(K) Cons_2_l(K) complex;
    
                for k = 1:K
                    Cons_1_l(k) = sqrt( mu_last(k) * xi_last(k) ) ...
                                  + sqrt( mu_last(k) / xi_last(k) ) * ( xi(k) - xi_last(k) ) / 2 ...
                                  + sqrt( xi_last(k) / mu_last(k) ) * ( mu(k) - mu_last(k) ) / 2;
                    Cons_1_r(k) = H(:,k)' * F(:,k);
                    Index_except_k = 1:K;
                    Index_except_k(Index_except_k == k) = [];
                    for m = Index_except_k        
                        Cons_2_l(k) = Cons_2_l(k) + F(:,m)' * H(:,k) * H(:,k)' * F(:,m);
                    end
                    Cons_2_l(k) = Cons_2_l(k) + P_n;
                end

                Obj = sum( log2( ones(K,1) + mu_last ) + ( mu - mu_last ) ./ ( log(2) * ( ones(K,1) + mu_last ) ) ) ...
                      - q * ( square_pos( norm( F, 'fro') ) / eta + P_fix );

                maximize ( real( Obj ) )
    
                subject to
    
                    real( Cons_1_l ) <= real( Cons_1_r );
                    real( Cons_2_l ) <= xi;
                    mu >= 2 .^ (R_min / B_W) - ones(K,1);
                    square_pos( norm( F, 'fro') ) <= P_t;
                    mu >= 0;
                    xi >= 0;
    
           cvx_end

           q = sum( log2( ones(K,1) + mu_last ) + ( mu - mu_last ) ./ ( log(2) * ( ones(K,1) + mu_last ) ) ) ...
               / ( trace( F' * F ) / eta + P_fix );

           iter_Dink = iter_Dink + 1;

        end

        R = zeros(K,1);           
        for k = 1:K
            IN = P_n - abs( H(:,k)' * F(:,k) ) ^ 2;
            for m = 1:K
                IN = IN + abs( H(:,k)' * F(:,m) ) ^ 2;
            end
            R(k) = B_W * log2( 1 + abs( H(:,k)' * F(:,k) ) ^ 2 / IN );
        end 
        EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

        if ( all( R >= R_min ) ) && ( norm( F, 'fro' ) ^ 2 <= P_t * ( 1 + 10^(-3) ) )
            Sign_cons = 1;
        else
            Sign_cons = 0;
        end

        if (Sign_cons == 1) && (EE > EE_save)
            EE_save = EE;
            F_save = F;
        end

        iter_SCA = iter_SCA + 1;

    end

    F = F_save;

end
