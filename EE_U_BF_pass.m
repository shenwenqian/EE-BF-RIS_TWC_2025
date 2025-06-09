function [EE, psi] = EE_U_BF_pass(B, Gamma_s, F, psi, Sigma, L_set, P_n, eta, P_fix, R_min, B_W)

    [N_B, N_R, L, K] = size(B);
    L_s = size(L_set, 1);
    L_r = L - L_s;

    if L_r < L

        L_rest = zeros(L_r, K);
        for k = 1:K
            L_rest(:,k) = setdiff( (1:1:L).', L_set(:,k) );
        end
    
        Sigma_s = zeros(L_s,K);
        Sigma_r = zeros(L_r,K);
        for k = 1:K
            Sigma_s(:,k) = Sigma(L_set(:,k), k);
            Sigma_r(:,k) = Sigma(L_rest(:,k), k);
        end   
    
        B_s = zeros(N_B,N_R,L_s,K);
        B_r = zeros(N_B,N_R,L_r,K);
        for k = 1:K
            B_s(:,:,:,k) = B(:, :, L_set(:,k), k);
            B_r(:,:,:,k) = B(:, :, L_rest(:,k), k);
        end
        
        J_I = zeros(L_r + 1, N_R, K, K);
        for k = 1:K   
            for m = 1:K
                for l_s = 1:L_s
                    J_I(1, 1:N_R, m, k) = J_I(1, 1:N_R, m, k) + Gamma_s(l_s, k) * Sigma_s(l_s, k)^(0.5) * F(:,m).' * conj( B_s(:, :, l_s, k) );
                end
                for l_r = 1:L_r
                    J_I(l_r + 1, 1:N_R, m, k) = Sigma_r(l_r, k)^(0.5) * F(:,m).' * conj( B_r(:, :, l_r, k) );
                end           
            end
        end
    
        J_S = zeros(L_r + 1, N_R, K);
        for k = 1:K
            J_S(:, :, k) = J_I( :, :, k, k );
        end

    else

        Sigma_r = Sigma;   
        B_r = B;

        J_I = zeros(L_r, N_R, K, K);
        for k = 1:K   
            for m = 1:K
                for l_r = 1:L_r
                    J_I(l_r, 1:N_R, m, k) = Sigma_r(l_r, k)^(0.5) * F(:,m).' * conj( B_r(:, :, l_r, k) );
                end           
            end
        end
    
        J_S = zeros(L_r, N_R, K);
        for k = 1:K
            J_S(:, :, k) = J_I( :, :, k, k );
        end

    end

    R = zeros(K,1);
    IN = P_n * ones(K,1);
    for k = 1:K
        Index_except_k = 1:K;
        Index_except_k(Index_except_k == k) = [];
        for m = Index_except_k        
            IN(k) = IN(k) + norm( J_I(:,:,m,k) * psi, 'fro' ) ^ 2;
        end
        R(k) = B_W * log2( 1 + norm( J_S(:,:,k) * psi, 'fro' ) ^ 2 / IN(k) );
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

    Iter_FP = 5000;
    iter_FP = 1;
    eta_FP = 10^(-5);
    EE = EE_save;
    EE_last = 10^(-10);

    Z = kron(ones(1, K + 1), psi);
    Lambda = zeros(N_R, K + 1);
    C = 10^(2);      
  
    while (EE / EE_last - 1 >= eta_FP) && (iter_FP <= Iter_FP) && (Sign_cons == 1)

        EE_last = EE;

        x = zeros(K,1);
        Y = zeros(size(J_S,1), K);
        for k = 1:K
            Inter_0 = P_n;
            Index_except_k = 1:K;
            Index_except_k(Index_except_k == k) = [];
            for m = Index_except_k        
                Inter_0 = Inter_0 + norm( J_I(:,:,m,k) * psi, 'fro' ) ^ 2;
            end
            x(k) = norm( J_S(:,:,k) * psi, 'fro' ) ^ 2 / Inter_0;
            Y(:,k) = sqrt( 1 + x(k) ) / ( norm( J_S(:,:,k) * psi, 'fro' ) ^ 2 + Inter_0 ) * J_S(:,:,k) * psi;
        end
            
        Inter_1 = C * (K + 1) * eye(N_R);
        Inter_3 = C * ( Z(:, K + 1) + Lambda(:, K + 1) );
        for k = 1:K
            Inter_2 = J_S(:,:,k)' * J_S(:,:,k);
            Index_except_k = 1:K;
            Index_except_k(Index_except_k == k) = [];
            for m = Index_except_k        
                Inter_2 = Inter_2 + J_I(:,:,m,k)' * J_I(:,:,m,k);
            end 
            Inter_1 = Inter_1 + Y(:,k)' * Y(:,k) * Inter_2;
            Inter_3 = Inter_3 + sqrt( 1 + x(k) ) * J_S(:,:,k)' * Y(:,k) + C * ( Z(:, k) + Lambda(:, k) );
        end
        psi = Inter_1 \ Inter_3;   

        for k = 1:K
            J_S_e = J_S(:,:,k)' * J_S(:,:,k);
            Inter_4 = zeros(N_R, N_R);
            Index_except_k = 1:K;
            Index_except_k(Index_except_k == k) = [];
            for m = Index_except_k        
                Inter_4 = Inter_4 + J_I(:,:,m,k)' * J_I(:,:,m,k);
            end                 
            J_S_e = J_S_e - ( 2 ^ ( R_min(k) / B_W ) - 1 ) * Inter_4;
            Sigma_e = ( 2 ^ ( R_min(k) / B_W ) - 1 ) * P_n;

            if ( psi - Lambda(:, k) )' * J_S_e * ( psi - Lambda(:, k) ) >= Sigma_e
                Z(:, k) = psi - Lambda(:, k);
            else
                [Q_k, Xi] = eig(J_S_e);
                Lambda_e = Q_k' * ( psi - Lambda(:, k) );

                if Xi(1, 1) < 0    
                    if Xi(N_R, N_R) > 0    
                        mu_max = - 1 / Xi(1, 1);
                        mu_min = - 1 / Xi(N_R, N_R);
                    else
                        mu_max = - 1 / Xi(N_R, N_R);
                        mu_min = 1 / Xi(N_R, N_R);
                    end
                else
                    mu_max = 1 / Xi(1, 1);
                    mu_min = - 1 / Xi(1, 1);
                end

                Iter_BIS = 100;
                iter_BIS = 1;                    
                eta_BIS = 10 ^ (-2);
                while (mu_max - mu_min >= eta_BIS) && (iter_BIS <= Iter_BIS)  
                    mu = (mu_max + mu_min) / 2;
                    f_mu = sum( diag(Xi) .* (abs( Lambda_e )) .^ 2 ...
                        ./ ( ones(N_R, 1) + mu * diag(Xi) ) .^ 2 ) - Sigma_e;
                    if f_mu > 0
                        mu_min = mu;
                    else
                        mu_max = mu;
                    end
                    iter_BIS = iter_BIS + 1;
                end
                mu = (mu_max + mu_min) / 2;
                Z_e = ( eye(N_R) + mu * Xi ) \ Lambda_e;
                Z(:, k) = Q_k * Z_e;
            end
        end

        Z(:, K + 1) = ( psi - Lambda(:, K + 1) ) ./ abs( psi - Lambda(:, K + 1) );

        Lambda = Lambda + Z - kron(ones(1, K + 1), psi);

        psi_test = psi ./ abs(psi);
        R = zeros(K,1);
        IN = P_n * ones(K,1);
        for k = 1:K
            Index_except_k = 1:K;
            Index_except_k(Index_except_k == k) = [];
            for m = Index_except_k        
                IN(k) = IN(k) + norm( J_I(:,:,m,k) * psi_test, 'fro' ) ^ 2;
            end
            R(k) = B_W * log2( 1 + norm( J_S(:,:,k) * psi_test, 'fro' ) ^ 2 / IN(k) );
        end
        EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

        if ( all( R >= R_min ) ) && ( EE > EE_save )
            EE_save = EE;
            psi_save = psi_test;
            Sign_cons = 1;
        else
            Sign_cons = 0;
        end
        
        iter_FP = iter_FP + 1;

    end

    EE = EE_save;
    psi = psi_save;

end
