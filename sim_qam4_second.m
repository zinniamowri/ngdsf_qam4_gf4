clear
rng(0)

%% === Simulation Parameters ===
Eb_No_db = 6:1:11;
T = 60;                    % max decoder iteration
eta = .3;                 % noise perturbation factor (placeholder)
w = 2;                    % decoder weighting parameter (placeholder)
flip_num = 5;              % max flipping per iteration

p = 2;                     % bits per GF(4) symbol (QPSK)
q = 2^p;

%% === Load GF(4) Arithmetic and H ===
load('arith_4.mat');       % add_mat, mul_mat, div_mat
load('204.102.3.6.16.mat'); % LDPC H matrix from GF(16 version)

% Convert GF(16) H → GF(4) structure
h_GF16 = full(h);
h_GF4 = h_GF16;
nz = h_GF4 ~= 0;
h_GF4(nz) = randi([1 q-1], 1, nnz(nz));  % assign random GF(4) labels (1–3)
h = h_GF4;

%% === Code Parameters ===
N = size(h,2);   % codeword length
M = size(h,1);   % number of parity checks
K = N - M;       % message length
R = K / N;       % code rate

fprintf('GF(4) LDPC code: N=%d, K=%d, R=%.3f\n', N, K, R);

% Parity-check adjacency list
CN_lst = cell(M,1);
for i = 1:M
    CN_lst{i,1} = find(h(i,:));
end

% Generate G from H
[G, Hsys, col_perm] = generate_G_from_H_GFq(h_GF4, add_mat, mul_mat);

 % --- Encode Random Info ---
 info_seq = randi([0 q-1], 1, K);
 code_seq = gf_mat_mul(info_seq, G, add_mat, mul_mat); % GF(4) codeword
 %Synd = gf_mat_mul(code_seq, h', add_mat, mul_mat);
 %assert(all(Synd == 0), 'Invalid codeword detected');
 Syndromes = decod_prod(code_seq,h,CN_lst, mul_mat, add_mat);

 %initializing vector size
 y=zeros(size(code_seq));
 hard_d_cmplx=zeros(size(code_seq));
 hard_d_gf4=zeros(size(code_seq));

%% === QPSK Constellation (Gray-coded, normalized) ===
qam4 = [-1-1i; -1+1i; 1+1i; 1-1i];
%qam4 = qam4 / sqrt(mean(abs(qam4).^2));   % Es = 1
gf4 = (0:q-1);
qam_binary_map = [0 0; 0 1; 1 1; 1 0];

avg_pow = qam4'*qam4/q; %sum of squared magnitudes of all symbols(the total power)/q.
nrm_fct=sqrt(avg_pow); %normalization factor
gf4 = (0:q-1); %GF field symbols
alph_bin =  logical(fliplr(dec2bin(gf4, p) - 48)); % symbols in binary

%avg_pow = mean(abs(qam4).^2);
%nrm_fct = sqrt(avg_pow);

%% === Initialize Storage ===
FE = zeros(length(Eb_No_db),1);
genFrame = zeros(length(Eb_No_db),1);
iters_cnr = zeros(length(Eb_No_db),1);
BE_Coded = zeros(length(Eb_No_db),1);
BE_unCoded = zeros(length(Eb_No_db),1);

targetFE = 500;          % stop criterion by frame error
max_gen  = 1e6;          % safety stop

%% === Main Eb/No Loop ===
for i = 1:length(Eb_No_db)

                    % for flipping metric (unused here)

    while (FE(i) < targetFE && genFrame(i) < max_gen)
        genFrame(i) = genFrame(i) + 1;

        c(1,1:N) = qam4(code_seq'+1,1); % codeword in complex
        avg_symbol_energy = 1;


        % compute noise std for this Eb/No
       Eb_No_linear = 10.^(Eb_No_db(i)/10);
       No = 1 / (p * Eb_No_linear * R);       % noise spectral density
       sigma0 = sqrt(No/2) * nrm_fct;
       nse_std = eta * sigma0;

        n = sigma0*(randn(1,N) + 1i*randn(1,N));
        y = c + n;

        for j=1:N
            distance=abs(qam4-y(j));
            [~,min_idx]= min(distance);
            hard_d_cmplx(j) = qam4(min_idx); % hard decision in complex
            hard_d_gf4(j) = gf4(min_idx);
        end
        
        errors = hard_d_cmplx ~= c; 
        n_errors_hard =sum(errors); %total number of symbol errors after hard decision made
              
        errors_uncoded_bit = zeros(1, K);

        %bit error calculation in hard decision
        for e = 1 : K          
                if hard_d_gf4(e)~=code_seq(e) 
                    s1 = dec2bin(code_seq(e),p);
                    s2 = dec2bin(hard_d_gf4(e),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_uncoded_bit(e) = errors_uncoded_bit(e) + num_diff_bit;                
                end
        end

        un_bit_error = sum(errors_uncoded_bit); % no of bit error in each frame
        BE_unCoded(i)=BE_unCoded(i)+un_bit_error; % total no. of bit error in total frame generated in current Eb/No
        

        % --- Call GF(4) decoder ---
        [seqgf, failed, l] = decodeQamMinDis_GF4( ...
            code_seq, hard_d_cmplx, hard_d_gf4, ...
            qam4, gf4, y, h, N, M, T, w, add_mat, mul_mat, div_mat, ...
            CN_lst, nse_std, qam_binary_map, flip_num, No);
        
        iters_cnr(i) = iters_cnr(i) + l;
        
       %bit error calculation for decoded sequence 
    
        errors_coded_bit = zeros(1, K);
            for g = 1 : K            
                if seqgf(g)~=code_seq(g)
                    s1 = dec2bin(code_seq(g),p);
                    s2 = dec2bin(seqgf(g),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_coded_bit(g) = errors_coded_bit(g) + num_diff_bit;                
                end     
            end

        bit_error = sum(errors_coded_bit); % no. of bit error in each frame
        BE_Coded(i)=BE_Coded(i)+bit_error; % total no. of bit error in total frame generated in the current Eb/No

        % --- Frame error ---
        if failed > 0
            FE(i) = FE(i) + 1;
        end

   end % end while

    fprintf('Eb/No = %.1f dB: BER (uncoded) = %.6e BER(coded)=%.6e\n', ...
        Eb_No_db(i), BE_unCoded(i) / (genFrame(i)*N*p), BE_Coded(i) / (genFrame(i)*N*p));
end

%% === Compute BERs ===
BERunCoded = BE_unCoded ./ (genFrame * N * p);
BERCoded = BE_Coded ./ (genFrame * N * p);

figure;
semilogy(Eb_No_db, BERCoded, 'gx-', 'LineWidth',1.5); hold on;
semilogy(Eb_No_db, BERunCoded, 'rx-', 'LineWidth',1.5);

grid on; box on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-7 1e-1]);
xlim([min(Eb_No_db)-0.5 max(Eb_No_db)+0.5]);
legend( ...
    sprintf('Proposed GF(4) Decoder (w=%.1f, \\eta=%.2f)', w, eta), ...
    'Uncoded', ...
    'Location','southwest');
title(sprintf('Coded vs Uncoded BER (GF(4) LDPC, Gray QPSK) — R=%.3f', R));
