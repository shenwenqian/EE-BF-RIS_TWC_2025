function [F] = EE_U_BF_act(B, Gamma_s, F, psi, Sigma, L_set, P_n, eta, P_fix, P_t, R_min, B_W)   
    
    [N_B, ~, L, K] = size(B);
    L_s = size(L_set, 1);    
    L_r = L - L_s;

    A = zeros(N_B,L,K);
    for k = 1:K
        for l = 1:L
            A(:,l,k) = B(:,:,l,k) * conj(psi);
        end
    end

    if L_r < L

        L_rest = zeros(L_r,K);
        for k = 1:K
            L_rest(:,k) = setdiff( (1:1:L).', L_set(:,k) );
        end
    
        Sigma_s = zeros(L_s,K);
        Sigma_r = zeros(L_r,K);
        for k = 1:K
            Sigma_s(:,k) = Sigma(L_set(:,k), k);
            Sigma_r(:,k) = Sigma(L_rest(:,k), k);
        end
        
        A_s = zeros(N_B,L_s,K);
        A_r = zeros(N_B,L_r,K);
        for k = 1:K
            A_s(:,:,k) = A(:, L_set(:,k), k);
            A_r(:,:,k) = A(:, L_rest(:,k), k);
        end
        
        O = zeros(L_r + 1, N_B, K);
        for k = 1:K
            O(1, 1:N_B, k) = Gamma_s(:,k).' * diag( (Sigma_s(:,k)).^(0.5) ) * A_s(:,:,k)';
            O(2:L_r + 1, 1:N_B, k) = diag( (Sigma_r(:,k)).^(0.5) ) * A_r(:,:,k)';
        end

    else

        Sigma_r = Sigma;
        A_r = A;
        
        O = zeros(L_r, N_B, K);
        for k = 1:K
            O(:, 1:N_B, k) = diag( (Sigma_r(:,k)).^(0.5) ) * A_r(:,:,k)';
        end

    end

    u = zeros(K,1);
    R = zeros(K,1);           
    IN = P_n * ones(K,1);
    for k = 1:K
        Index_except_k = 1:K;
        Index_except_k(Index_except_k == k) = [];
        for m = Index_except_k 
            IN(k) = IN(k) + norm( O(:,:,k) * F(:,m) ) ^ 2;
        end                    
        u(k) = 1 + norm( O(:,:,k) * F(:,k) ) ^ 2 / IN(k);
        R(k) = B_W * log2( u(k) );
    end  
    EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

    if ( all( R >= R_min ) ) && ( norm( F, 'fro' ) ^ 2 < P_t * ( 1 + 10^(-3) ) )
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
        F_last = F;
        u_last = u;

        if L_r < L

            cvx_clear
            cvx_begin quiet
                variable F_e(N_B, K) complex;
                variable u_e(K);
                variable t_e(K);
                variable theta;
                expressions e_1((L_r + 1) * (K - 1) + 2, K) e_2((L_r + 1) * (K - 1) + 3, K) e_3(2, K) e_4(N_B * K + 2) e_5(N_B * K)...
                            Cons_1_l(K) Cons_1_r(K) Cons_2_l(K) Cons_2_r(K) Cons_3_l(K) Cons_3_r(K)...
                            func_f(K) Val_fOfe(K) Arr_Ofe((L_r + 1) * (K - 1), K) F_e_vec(N_B * K) complex;
                for k = 1:K
    
                    func_f(k) = P_n / u_last(k);
                    for m = 1:K
                        func_f(k) = func_f(k) + norm( O(:,:,k) * F_last(:,m) ) ^ 2 / u_last(k);
                        Val_fOfe(k) = Val_fOfe(k) + F_last(:,m)' * O(:,:,k)' * O(:,:,k) * F_e(:,m);
                        if m < k
                            Arr_Ofe((L_r + 1) * (m - 1) + 1:(L_r + 1) * m, k) = O(:,:,k) * F_e(:,m);
                        end
                        if m > k
                            Arr_Ofe((L_r + 1) * (m - 2) + 1:(L_r + 1) * (m - 1), k) = O(:,:,k) * F_e(:,m);
                        end
                    end
    
                    e_1(1:(L_r + 1) * (K - 1), k) = Arr_Ofe(:,k);
                    e_1((L_r + 1) * (K - 1) + 1, k) = ( - func_f(k) * u_e(k) + 2 * Val_fOfe(k) ) / u_last(k) - theta / 4;
                    e_1((L_r + 1) * (K - 1) + 2, k) = theta * sqrt( P_n * ( u_last(k) - 2 ) / u_last(k) );
                    Cons_1_l(k) = norm( e_1(:,k) );
                    Cons_1_r(k) = ( - func_f(k) * u_e(k) + 2 * Val_fOfe(k) ) / u_last(k) + theta / 4;
    
                    e_2(1:(L_r + 1) * (K - 1), k) = 2 ^ ( R_min(k) / 2 / B_W ) * Arr_Ofe(:,k);
                    e_2((L_r + 1) * (K - 1) + 1, k) = ( theta * P_n + Val_fOfe(k) ) / u_last(k) - u_e(k) / 2;
                    e_2((L_r + 1) * (K - 1) + 2, k) = sqrt( func_f(k) / u_last(k) ) * u_e(k);
                    e_2((L_r + 1) * (K - 1) + 3, k) = 2 ^ ( R_min(k) / 2 / B_W ) * theta * sqrt(P_n);
                    Cons_2_l(k) = norm( e_2(:,k) );
                    Cons_2_r(k) = ( theta * P_n + Val_fOfe(k) ) / u_last(k) + u_e(k) / 2;
    
                    e_3(1,k) = theta * log2( exp(1) * u_last(k) ) - t_e(k) - u_e(k) / 4;
                    e_3(2,k) = theta * sqrt( log2( exp(1) ) * u_last(k) );
                    Cons_3_l(k) = norm( e_3(:,k) );
                    Cons_3_r(k) = theta * log2( exp(1) * u_last(k) ) - t_e(k) + u_e(k) / 4;
    
                    F_e_vec(N_B * (k - 1) + 1: N_B * k) = F_e(:,k);
    
                end
    
                e_4(1:N_B * K) = F_e_vec / sqrt(eta);
                e_4(N_B * K + 1) = theta - 1 / 4;
                e_4(N_B * K + 2) = theta * sqrt(P_fix);
    
                e_5 = F_e_vec;
    
                maximize ( sum(t_e) )
    
                subject to
    
                    real( Cons_1_l ) <= real( Cons_1_r );
                    real( Cons_2_l ) <= real( Cons_2_r );
                    real( Cons_3_l ) <= real( Cons_3_r );
                    norm( e_4 ) <= theta + 1 / 4;
                    norm( e_5 ) <= theta * sqrt(P_t);
                    u_e >= 0;
                    t_e >= 0;
                    theta >= 0; 
    
           cvx_end

        else

            cvx_clear
            cvx_begin quiet
                variable F_e(N_B, K) complex;
                variable u_e(K);
                variable t_e(K);
                variable theta;
                expressions e_1((L_r) * (K - 1) + 2, K) e_2((L_r) * (K - 1) + 3, K) e_3(2, K) e_4(N_B * K + 2) e_5(N_B * K)...
                            Cons_1_l(K) Cons_1_r(K) Cons_2_l(K) Cons_2_r(K) Cons_3_l(K) Cons_3_r(K)...
                            func_f(K) Val_fOfe(K) Arr_Ofe((L_r) * (K - 1), K) F_e_vec(N_B * K) complex;
                for k = 1:K
    
                    func_f(k) = P_n / u_last(k);
                    for m = 1:K
                        func_f(k) = func_f(k) + norm( O(:,:,k) * F_last(:,m) ) ^ 2 / u_last(k);
                        Val_fOfe(k) = Val_fOfe(k) + F_last(:,m)' * O(:,:,k)' * O(:,:,k) * F_e(:,m);
                        if m < k
                            Arr_Ofe(L_r * (m - 1) + 1:L_r * m, k) = O(:,:,k) * F_e(:,m);
                        end
                        if m > k
                            Arr_Ofe(L_r * (m - 2) + 1:L_r * (m - 1), k) = O(:,:,k) * F_e(:,m);
                        end
                    end
    
                    e_1(1:L_r * (K - 1), k) = Arr_Ofe(:,k);
                    e_1(L_r * (K - 1) + 1, k) = ( - func_f(k) * u_e(k) + 2 * Val_fOfe(k) ) / u_last(k) - theta / 4;
                    e_1(L_r * (K - 1) + 2, k) = theta * sqrt( P_n * ( u_last(k) - 2 ) / u_last(k) );
                    Cons_1_l(k) = norm( e_1(:,k) );
                    Cons_1_r(k) = ( - func_f(k) * u_e(k) + 2 * Val_fOfe(k) ) / u_last(k) + theta / 4;
    
                    e_2(1:L_r * (K - 1), k) = 2 ^ ( R_min(k) / 2 / B_W ) * Arr_Ofe(:,k);
                    e_2(L_r * (K - 1) + 1, k) = ( theta * P_n + Val_fOfe(k) ) / u_last(k) - u_e(k) / 2;
                    e_2(L_r * (K - 1) + 2, k) = sqrt( func_f(k) / u_last(k) ) * u_e(k);
                    e_2(L_r * (K - 1) + 3, k) = 2 ^ ( R_min(k) / 2 / B_W ) * theta * sqrt(P_n);
                    Cons_2_l(k) = norm( e_2(:,k) );
                    Cons_2_r(k) = ( theta * P_n + Val_fOfe(k) ) / u_last(k) + u_e(k) / 2;
    
                    e_3(1,k) = theta * log2( exp(1) * u_last(k) ) - t_e(k) - u_e(k) / 4;
                    e_3(2,k) = theta * sqrt( log2( exp(1) ) * u_last(k) );
                    Cons_3_l(k) = norm( e_3(:,k) );
                    Cons_3_r(k) = theta * log2( exp(1) * u_last(k) ) - t_e(k) + u_e(k) / 4;
    
                    F_e_vec(N_B * (k - 1) + 1: N_B * k) = F_e(:,k);
    
                end
    
                e_4(1:N_B * K) = F_e_vec / sqrt(eta);
                e_4(N_B * K + 1) = theta - 1 / 4;
                e_4(N_B * K + 2) = theta * sqrt(P_fix);
    
                e_5 = F_e_vec;
    
                maximize ( sum(t_e) )
    
                subject to
    
                    real( Cons_1_l ) <= real( Cons_1_r );
                    real( Cons_2_l ) <= real( Cons_2_r );
                    real( Cons_3_l ) <= real( Cons_3_r );
                    norm( e_4 ) <= theta + 1 / 4;
                    norm( e_5 ) <= theta * sqrt(P_t);
                    u_e >= 0;
                    t_e >= 0;
                    theta >= 0; 
    
           cvx_end

        end

        F = F_e / theta;

        u = zeros(K,1);
        R = zeros(K,1);           
        IN = P_n * ones(K,1);
        for k = 1:K
            Index_except_k = 1:K;
            Index_except_k(Index_except_k == k) = [];
            for m = Index_except_k 
                IN(k) = IN(k) + norm( O(:,:,k) * F(:,m) ) ^ 2;
            end                    
            u(k) = 1 + norm( O(:,:,k) * F(:,k) ) ^ 2 / IN(k);
            R(k) = B_W * log2( u(k) );
       end  
       EE = sum(R) / ( norm( F, 'fro' ) ^ 2 / eta + P_fix );

       if ( all( R >= R_min ) ) && ( norm( F, 'fro' ) ^ 2 < P_t * ( 1 + 10^(-3) ) )
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
