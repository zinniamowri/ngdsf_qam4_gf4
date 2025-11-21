function [seqgf,failed_init,l] = decodeQamMinDis_GF4(code_seq, hard_d_cmplx, hard_d_gf4, ...
    qam4, gf4, y, h, N, M, T, w, add_mat, mul_mat, div_mat, ...
    CN_lst, nse_std, qam_binary_map, flip_num, No)

   
    d_deco_0   = hard_d_cmplx;        
    seqgf      = hard_d_gf4;
    d_dec_f    = hard_d_cmplx;        % complex reference for metric E
    l          = 0;
    d_gf4      = hard_d_gf4;

    S  = decod_prod(hard_d_gf4, h, CN_lst, mul_mat, add_mat);
    Sb = double(S == 0);
    num_satisfied = sum(Sb);
    failed_init   = M - num_satisfied;

    % Track best-so-far
    min_failed  = failed_init;
    best_d_gf4  = d_gf4;

    Hb = double(h > 0);

    % --- Iterative flipping 
    while (l < T)
        l = l + 1;

        Sb       = double(S == 0);
        WSH      = w * Sb * Hb;                 
        E        = (-abs(y - d_dec_f).^2 /No) + WSH + nse_std * randn(1, N);
        %E        = (-abs(y - d_dec_f).^2 /No) + WSH ; %without noise
        
        % pick "flip_num" most unreliable VNs
        [~, idx] = mink(E, flip_num);

        temp_gf4 = d_gf4;

        % Local random search in nearest constellation neighbors
        for j = 1:length(idx)
            % current complex point around which to pick neighbors
            current_symbol_c = d_dec_f(idx(j));

            % nearest neighbors in QPSK (there are only 4 points total)
            distance          = abs(qam4 - current_symbol_c);
            [~, sorted_idx]   = mink(distance, 4);   % all 4 are valid
            closest_symbols   = gf4(sorted_idx);

            % choose randomly among the 2 closest 
            cand_idx = randi([1, 2]);
            new_symbol = closest_symbols(cand_idx);

            if new_symbol == temp_gf4(idx(j))
                cand_idx = 3 - cand_idx; % flip 1 <-> 2
                new_symbol = closest_symbols(cand_idx);
            end

            temp_gf4(idx(j)) = new_symbol;
        end

        % Check new syndrome
        S_temp  = decod_prod(temp_gf4, h, CN_lst, mul_mat, add_mat);
        Sb_temp = double(S_temp == 0);
        num_satisfied_temp = sum(Sb_temp);
        failed_temp        = M - num_satisfied_temp;

        if failed_init == 0
            break;
        end

        if failed_temp < failed_init
            d_gf4      = temp_gf4;
            failed_init = failed_temp;

            if failed_temp < min_failed
                min_failed = failed_temp;
                best_d_gf4 = temp_gf4;
            end

            % Move forward: update current S as well
            S = S_temp;
        end
    end

    seqgf = best_d_gf4;

    % --------------------------
    % Post-processing 
    % --------------------------
    S  = decod_prod(seqgf, h, CN_lst, mul_mat, add_mat);
    Sb = double(S == 0);

    unsatisfied_indices = find(Sb == 0);

    % Collect VN candidates touching unsatisfied CNs
    VN_candidates = [];
    for ii = 1:length(unsatisfied_indices)
        idx = unsatisfied_indices(ii);
        VN_list = CN_lst{idx};
        VN_candidates = [VN_candidates, VN_list]; 
    end
    VN_candidates = unique(VN_candidates);

    best_seq    = seqgf;
    best_failed = M - sum(Sb);

    % Try all alternatives 0..3 for each VN candidate
    for ii = 1:length(VN_candidates)
        vn_idx = VN_candidates(ii);
        current_symbol = seqgf(vn_idx);

        for test_sym = 0:3
            if test_sym ~= current_symbol
                temp_seq = seqgf;
                temp_seq(vn_idx) = test_sym;

                S_check  = decod_prod(temp_seq, h, CN_lst, mul_mat, add_mat);
                Sb_check = double(S_check == 0);
                failed_now = M - sum(Sb_check);

                if failed_now < best_failed
                    best_failed = failed_now;
                    best_seq    = temp_seq;
                end
            end
        end
    end

    seqgf = best_seq;
end
