
clear all; close all; clc;

fprintf('==================================================\n');
fprintf('禁忌搜索改进的人工地震动时程生成程序\n');
fprintf('==================================================\n\n');

%% 第一步：参数配置（完全按照论文设置）
% ================= 基础参数 =================
params.dt = 0.01;               % 时间步长 (s)
params.N = 3000;                % 数据点数（对应30s总时长）
params.T_total = params.N * params.dt; % 总时长 (s)
params.damping = 0.05;          % 阻尼比 (5%)

% ================= 包络函数参数 =================
params.t1 = 5.0;                % 上升段时间(s)（论文值）
params.t2 = 13.0;               % 平稳段结束时间(s)（论文值）
params.c = 0.3;                 % 衰减系数（论文值）

% ================= 目标反应谱参数 =================
params.PGA_target = 1.0;        % 目标峰值加速度 (g)（论文值）
params.n_control_points = 30;   % 反应谱控制点数

% ================= 功率谱参数 =================
params.psd_envelope_threshold = 0.8;  % 80%包络标准（论文要求）
params.psd_freq_range = [0.3, 24];    % 功率谱验证频率范围 (Hz)（论文要求）

% ================= 迭代参数 =================
params.max_iterations = 30000;     % 最大迭代次数
params.max_below_points = 5;     % 低于目标谱的最大控制点数（论文要求）
params.error_tolerance = 0.10;   % 误差容限 (10%)（论文要求）

% ================= 频域参数 =================
params.Nfft = 2^nextpow2(params.N);  % FFT点数

% ================= 禁忌搜索参数 =================
params.tabu_tenure = 15;         % 禁忌长度（禁止重复调整的次数）
params.diversification_threshold = 30; % 多样化阈值（连续无改进次数）
params.intensification_threshold = 10; % 强化阈值（连续改进次数）
params.adaptive_relaxation = true;     % 自适应松弛因子

fprintf('参数配置完成（使用禁忌搜索改进）:\n');
fprintf('  总时长: %.0f s（论文值）\n', params.T_total);
fprintf('  上升段: %.1f s，平稳段: %.1f s（论文值）\n', params.t1, params.t2-params.t1);
fprintf('  阻尼比: %.2f\n', params.damping);
fprintf('  目标PGA: %.2f g\n', params.PGA_target);
fprintf('  禁忌长度: %d 次迭代\n', params.tabu_tenure);
fprintf('  多样化阈值: %d 次无改进\n', params.diversification_threshold);

%% 第二步：生成目标反应谱和功率谱（RG1.60标准谱）
fprintf('\n==================================================\n');
fprintf('生成目标反应谱和功率谱（RG1.60标准）\n');
fprintf('==================================================\n');

% 生成时间向量
time = (0:params.N-1)' * params.dt;

% 生成目标反应谱（RG1.60标准谱）
T_target = logspace(-2, 1, params.n_control_points)';
T_target = T_target(T_target >= 0.05 & T_target <= 3.0);
f_target = 1./T_target;  % 转换为频率

% RG1.60标准反应谱形状
T0 = 0.1;  % 起始周期
Tc = 0.5;  % 特征周期
Sa_target = zeros(size(T_target));

for i = 1:length(T_target)
    T = T_target(i);
    if T <= T0
        Sa_target(i) = params.PGA_target;
    elseif T <= Tc
        Sa_target(i) = params.PGA_target * (1 + 2.0 * (T - T0)/(Tc - T0));
    else
        Sa_target(i) = params.PGA_target * 3.0 * (Tc/T);
    end
end

% 确保目标谱值大于0
Sa_target = max(Sa_target, params.PGA_target * 0.001);

% 计算目标功率谱密度（RG1.60标准功率谱）
freq_std = logspace(log10(0.1), log10(30), 200)';
psd_std = calculate_RG1_60_power_spectrum_simple(freq_std);

fprintf('目标谱生成完成:\n');
fprintf('  反应谱控制点: %d个，周期范围: %.3f - %.3f s\n', length(T_target), min(T_target), max(T_target));
fprintf('  功率谱频率点: %d个，频率范围: %.1f - %.1f Hz\n', length(freq_std), min(freq_std), max(freq_std));

%% 第三步：生成初始人工时程（按照论文方法）
fprintf('\n==================================================\n');
fprintf('生成初始人工时程（按照论文方法）\n');
fprintf('==================================================\n');

fprintf('1. 计算目标功率谱对应的傅里叶幅值（论文公式(2)）...\n');

% 频率向量
n_freq = floor(params.Nfft/2) + 1;
freq_fft = (0:n_freq-1)' / (params.Nfft * params.dt);
delta_f = freq_fft(2) - freq_fft(1);

% 插值得到目标功率谱密度在等间距频率点的值
psd_interp = interp1(freq_std, psd_std, freq_fft, 'pchip', 0);

% 计算傅里叶幅值（论文公式(2)）
A_i = sqrt(4 * psd_interp * delta_f);  % 幅值
phi_i = 2 * pi * rand(n_freq, 1);      % 随机相位

fprintf('2. 使用三角级数叠加法生成平稳加速度时程（论文公式(2)）...\n');

% 对称处理频率分量
A_full = [A_i; flipud(A_i(2:end-1))];
phi_full = [phi_i; -flipud(phi_i(2:end-1))];

% 使用IFFT生成时程
complex_spectrum = A_full .* exp(1i * phi_full);
acc_stationary_full = real(ifft(complex_spectrum, params.Nfft)) * sqrt(params.Nfft);
acc_stationary = acc_stationary_full(1:params.N);

fprintf('3. 应用非平稳强度包络函数（论文公式(3)(4)）...\n');

% 应用包络函数
envelope = get_envelope_vector_simple(time, params);
acc_initial = acc_stationary .* envelope;

% 调整初始时程的PGA
current_PGA = max(abs(acc_initial));
if current_PGA > 0
    acc_initial = acc_initial * (params.PGA_target / current_PGA);
end

% 计算初始时程的反应谱
Sa_initial = calculate_response_spectrum_simple(acc_initial, params.dt, params.damping, T_target);

% 计算初始功率谱
[freq_gen, psd_gen] = calculate_power_spectrum_simple(acc_initial, params.dt);

% 计算初始误差统计
[below_initial, max_error_initial, mean_error_initial] = compute_error_stats_simple(Sa_initial, Sa_target);

% 计算初始功率谱包络度
envelope_initial = compute_psd_envelope_simple(freq_gen, psd_gen, freq_std, psd_std, params);

fprintf('初始时程生成完成:\n');
fprintf('  低于目标谱点数: %d/%d（要求≤%d）\n', below_initial, length(T_target), params.max_below_points);
fprintf('  最大误差: %.2f%%（要求≤%.0f%%）\n', max_error_initial*100, params.error_tolerance*100);
fprintf('  平均误差: %.2f%%\n', mean_error_initial*100);
fprintf('  功率谱包络度: %.1f%%（要求≥%.0f%%）\n', envelope_initial*100, params.psd_envelope_threshold*100);
fprintf('  PGA: %.4f g\n', max(abs(acc_initial)));

%% 第四步：禁忌搜索算法改进的迭代流程
fprintf('\n==================================================\n');
fprintf('禁忌搜索算法改进的迭代流程\n');
fprintf('==================================================\n');

fprintf('禁忌搜索策略:\n');
fprintf('  1. 记录最近调整的频率成分，避免重复调整（禁忌表）\n');
fprintf('  2. 连续无改进时进行多样化搜索\n');
fprintf('  3. 连续改进时进行强化搜索\n');
fprintf('  4. 自适应调整松弛因子\n\n');

fprintf('开始禁忌搜索迭代:\n');
fprintf('迭代  低于谱点数  最大误差%%  平均误差%%  包络度%%  状态\n');
fprintf('------------------------------------------------------------\n');

% 初始化
acc_current = acc_initial;
Sa_current = Sa_initial;
[below_current, max_error_current, mean_error_current] = compute_error_stats_simple(Sa_current, Sa_target);
[freq_gen_current, psd_gen_current] = calculate_power_spectrum_simple(acc_current, params.dt);
envelope_current = compute_psd_envelope_simple(freq_gen_current, psd_gen_current, freq_std, psd_std, params);

% 计算综合评分（越低越好）
current_score = compute_composite_score(below_current, max_error_current, mean_error_current, envelope_current, params);

% 初始化禁忌表
tabu_list = zeros(n_freq, 1);  % 禁忌表，记录每个频率分量的禁忌剩余次数
tabu_list_history = cell(params.max_iterations, 1); % 记录禁忌表历史

% 初始化最佳解
best_acc = acc_current;
best_score = current_score;
best_iteration = 1;
no_improve_count = 0;  % 连续无改进次数
improve_count = 0;     % 连续改进次数

% 初始化迭代历史
iteration_history = zeros(params.max_iterations, 7); % [迭代次数, 低于点数, 最大误差, 平均误差, 包络度, 综合评分, 状态]
iteration_history(1, :) = [1, below_current, max_error_current*100, mean_error_current*100, ...
                           envelope_current*100, current_score, 0]; % 状态0:初始

% 初始化自适应松弛因子
relaxation_factor = 0.3;  % 初始松弛因子
min_relaxation = 0.1;     % 最小松弛因子
max_relaxation = 0.8;     % 最大松弛因子

% 初始化多样化搜索参数
diversification_count = 0;  % 多样化搜索计数器

for iter = 2:params.max_iterations
    t_start = tic;
    
    % 更新禁忌表（所有禁忌剩余次数减1）
    tabu_list = max(tabu_list - 1, 0);
    
    % 根据当前状态选择搜索策略
    if no_improve_count >= params.diversification_threshold
        % 多样化搜索：连续无改进次数过多，进行随机扰动
        search_strategy = 'diversification';
        strategy_code = 2;
        diversification_count = diversification_count + 1;
        
        % 增加松弛因子，扩大搜索范围
        relaxation_factor = min(max_relaxation, relaxation_factor * 1.2);
        
        fprintf('  [多样化搜索%d，松弛因子: %.2f]\n', diversification_count, relaxation_factor);
    elseif improve_count >= params.intensification_threshold
        % 强化搜索：连续改进，进行精细调整
        search_strategy = 'intensification';
        strategy_code = 3;
        
        % 减小松弛因子，进行精细调整
        relaxation_factor = max(min_relaxation, relaxation_factor * 0.8);
        
        fprintf('  [强化搜索，松弛因子: %.2f]\n', relaxation_factor);
    else
        % 正常搜索
        search_strategy = 'normal';
        strategy_code = 1;
    end
    
    % 备份当前傅里叶系数
    acc_fft = fft(acc_current, params.Nfft);
    F = acc_fft(1:n_freq);
    F_backup = F;
    
    % 根据搜索策略调整傅里叶系数
    switch search_strategy
        case 'normal'
            % 正常搜索：按照论文方法进行调整
            F = normal_search_adjustment(F, acc_current, freq_fft, T_target, f_target, ...
                                        Sa_target, params, relaxation_factor, tabu_list);
            
        case 'diversification'
            % 多样化搜索：随机调整一些频率分量
            F = diversification_search_adjustment(F, freq_fft, tabu_list, diversification_count);
            
        case 'intensification'
            % 强化搜索：针对当前误差最大的频率进行精细调整
            F = intensification_search_adjustment(F, acc_current, freq_fft, T_target, f_target, ...
                                                 Sa_target, Sa_current, params, relaxation_factor);
    end
    
    % 更新禁忌表：记录本次调整的频率分量
    adjusted_freq_indices = find(abs(F - F_backup) > 1e-6);
    if ~isempty(adjusted_freq_indices)
        tabu_list(adjusted_freq_indices) = params.tabu_tenure;
    end
    
    % 保存禁忌表历史
    tabu_list_history{iter} = tabu_list;
    
    % 傅里叶逆变换得到调整后的时程
    F_full = [F; conj(flipud(F(2:end-1)))];
    acc_new = real(ifft(F_full, params.Nfft));
    acc_new = acc_new(1:params.N);
    
    % 应用非平稳强度包络
    acc_new = acc_new .* envelope;
    
    % 检查PGA限制（论文公式(9)）
    for i = 1:length(acc_new)
        if abs(acc_new(i)) > params.PGA_target * envelope(i)
            acc_new(i) = sign(acc_new(i)) * params.PGA_target * envelope(i);
        end
    end
    
    % 计算新时程的性能指标
    Sa_new = calculate_response_spectrum_simple(acc_new, params.dt, params.damping, T_target);
    [below_new, max_error_new, mean_error_new] = compute_error_stats_simple(Sa_new, Sa_target);
    [freq_gen_new, psd_gen_new] = calculate_power_spectrum_simple(acc_new, params.dt);
    envelope_new = compute_psd_envelope_simple(freq_gen_new, psd_gen_new, freq_std, psd_std, params);
    
    % 计算综合评分
    new_score = compute_composite_score(below_new, max_error_new, mean_error_new, envelope_new, params);
    
    % 判断是否接受新解（禁忌搜索准则）
    accept_new = false;
    if new_score < current_score
        % 新解优于当前解，接受
        accept_new = true;
    elseif ~any(tabu_list(adjusted_freq_indices) > 0) && rand() < 0.1
        % 非禁忌且以一定概率接受劣解（避免陷入局部最优）
        accept_new = true;
    end
    
    if accept_new
        % 接受新解
        acc_current = acc_new;
        Sa_current = Sa_new;
        below_current = below_new;
        max_error_current = max_error_new;
        mean_error_current = mean_error_new;
        envelope_current = envelope_new;
        current_score = new_score;
        
        % 更新改进计数器
        if new_score < best_score
            % 找到新的最优解
            best_acc = acc_new;
            best_score = new_score;
            best_iteration = iter;
            no_improve_count = 0;
            improve_count = improve_count + 1;
            
            % 输出找到新最优解的信息
            fprintf('  [找到新最优解! 评分: %.2f -> %.2f]\n', iteration_history(iter-1, 6), best_score);
        else
            % 改进但未达到最优
            no_improve_count = no_improve_count + 1;
            improve_count = improve_count + 1;
        end
    else
        % 拒绝新解，保持当前解
        no_improve_count = no_improve_count + 1;
        improve_count = 0;
        
        % 如果连续拒绝次数过多，减小禁忌长度
        if no_improve_count > 20
            params.tabu_tenure = max(5, params.tabu_tenure - 1);
        end
    end
    
    % 记录迭代历史
    iteration_history(iter, :) = [iter, below_current, max_error_current*100, ...
                                 mean_error_current*100, envelope_current*100, ...
                                 current_score, strategy_code];
    
    % 输出迭代信息
    if mod(iter, 20) == 0 || iter == 2 || iter == params.max_iterations
        status_str = get_status_string(strategy_code);
        fprintf('%4d      %3d       %6.2f      %6.2f     %6.2f     %s\n', ...
            iter, below_current, max_error_current*100, mean_error_current*100, ...
            envelope_current*100, status_str);
    end
    
    % 检查是否满足论文要求
    if below_current <= params.max_below_points && ...
       max_error_current <= params.error_tolerance && ...
       envelope_current >= params.psd_envelope_threshold
        fprintf('\n✓ 同时满足反应谱和功率谱要求！\n');
        fprintf('  低于目标谱点数: %d (要求≤%d)\n', below_current, params.max_below_points);
        fprintf('  最大误差: %.2f%% (要求≤%.0f%%)\n', max_error_current*100, params.error_tolerance*100);
        fprintf('  功率谱包络度: %.2f%% (要求≥%.0f%%)\n', envelope_current*100, params.psd_envelope_threshold*100);
        break;
    end
    
    % 如果连续无改进次数过多，重置禁忌表
    if no_improve_count > params.diversification_threshold * 2
        tabu_list = zeros(n_freq, 1);
        no_improve_count = 0;
        fprintf('  [重置禁忌表]\n');
    end
    
    % 自适应调整松弛因子
    if params.adaptive_relaxation
        if improve_count > 5
            % 连续改进，减小松弛因子进行精细搜索
            relaxation_factor = max(min_relaxation, relaxation_factor * 0.95);
        elseif no_improve_count > 10
            % 连续无改进，增大松弛因子扩大搜索范围
            relaxation_factor = min(max_relaxation, relaxation_factor * 1.05);
        end
    end
end

% 使用最佳解
if best_score < current_score
    acc_current = best_acc;
    Sa_current = calculate_response_spectrum_simple(best_acc, params.dt, params.damping, T_target);
    [below_current, max_error_current, mean_error_current] = compute_error_stats_simple(Sa_current, Sa_target);
    [freq_gen_current, psd_gen_current] = calculate_power_spectrum_simple(best_acc, params.dt);
    envelope_current = compute_psd_envelope_simple(freq_gen_current, psd_gen_current, freq_std, psd_std, params);
    
    fprintf('\n使用第%d次迭代的最佳解\n', best_iteration);
end

if iter == params.max_iterations
    fprintf('\n达到最大迭代次数\n');
end

fprintf('\n禁忌搜索迭代完成，总迭代次数: %d\n', min(iter, params.max_iterations));
fprintf('最佳解在第%d次迭代获得，综合评分: %.2f\n', best_iteration, best_score);

%% 第五步：最终时程生成与评估
fprintf('\n==================================================\n');
fprintf('最终时程生成与评估\n');
fprintf('==================================================\n');

% 计算速度和位移时程
vel_final = cumtrapz(acc_current) * params.dt;
disp_final = cumtrapz(vel_final) * params.dt;
disp_final = detrend(disp_final);  % 去除基线漂移

% 计算峰值
PGA = max(abs(acc_current));
PGV = max(abs(vel_final)) * 100;  % 转换为cm/s
PGD = max(abs(disp_final)) * 100; % 转换为cm

% 计算最终反应谱
Sa_final = calculate_response_spectrum_simple(acc_current, params.dt, params.damping, T_target);
[below_final, max_error_final, mean_error_final] = compute_error_stats_simple(Sa_final, Sa_target);

% 计算最终功率谱包络度
[freq_gen, psd_gen] = calculate_power_spectrum_simple(acc_current, params.dt);
envelope_final = compute_psd_envelope_simple(freq_gen, psd_gen, freq_std, psd_std, params);

fprintf('最终结果:\n');
fprintf('  反应谱平均误差: %.3f%%\n', mean_error_final*100);
fprintf('  反应谱最大误差: %.3f%%\n', max_error_final*100);
fprintf('  低于目标谱点数: %d/%d（要求≤%d）\n', below_final, length(T_target), params.max_below_points);
fprintf('  功率谱包络度: %.1f%%（要求≥%.0f%%）\n', envelope_final*100, params.psd_envelope_threshold*100);
fprintf('  峰值加速度(PGA): %.4f g\n', PGA);
fprintf('  峰值速度(PGV): %.2f cm/s\n', PGV);
fprintf('  峰值位移(PGD): %.2f cm\n', PGD);

% 判断是否满足论文要求
if below_final <= params.max_below_points
    fprintf('  ✓ 低于目标谱点数满足要求\n');
else
    fprintf('  ✗ 低于目标谱点数不满足要求\n');
end

if max_error_final <= params.error_tolerance
    fprintf('  ✓ 最大误差满足要求\n');
else
    fprintf('  ✗ 最大误差不满足要求\n');
end

if envelope_final >= params.psd_envelope_threshold
    fprintf('  ✓ 功率谱包络度满足要求\n');
else
    fprintf('  ✗ 功率谱包络度不满足要求\n');
end

if below_final <= params.max_below_points && ...
   max_error_final <= params.error_tolerance && ...
   envelope_final >= params.psd_envelope_threshold
    fprintf('\n✓ 全部满足论文要求！\n');
else
    fprintf('\n✗ 未完全满足论文要求\n');
end

%% 第六步：结果可视化
fprintf('\n==================================================\n');
fprintf('结果可视化\n');
fprintf('==================================================\n');

create_results_plots_tabu(time, acc_current, vel_final, disp_final, ...
    T_target, Sa_target, Sa_final, ...
    freq_gen, psd_gen, freq_std, psd_std, ...
    iteration_history, min(iter, params.max_iterations), ...
    params, PGA, PGV, PGD, mean_error_final, envelope_final, best_iteration);

%% 第七步：保存结果
fprintf('\n保存结果...\n');
save_results_files_tabu(time, acc_current, vel_final, disp_final, ...
    T_target, Sa_target, Sa_final, params, ...
    mean_error_final, envelope_final, PGA, PGV, PGD, ...
    iteration_history, best_iteration);

fprintf('\n程序完成！\n');
fprintf('==================================================\n');

%% 辅助函数定义
% ============================================================

function envelope = get_envelope_vector_simple(time, params)
    % 获取包络函数向量
    envelope = zeros(size(time));
    
    for i = 1:length(time)
        t = time(i);
        if t <= params.t1
            envelope(i) = (t/params.t1)^2;
        elseif t <= params.t2
            envelope(i) = 1.0;
        else
            envelope(i) = exp(-params.c * (t - params.t2));
        end
    end
end

function Sa = calculate_response_spectrum_simple(acc, dt, damping, T_points)
    % 简化的反应谱计算方法
    nT = length(T_points);
    Sa = zeros(nT, 1);
    N = length(acc);
    
    for i = 1:nT
        T = T_points(i);
        if T <= 0
            Sa(i) = 0;
            continue;
        end
        
        omega = 2 * pi / T;
        xi = damping;
        
        % 使用简化的计算方法
        nfft = 2^nextpow2(N);
        acc_fft = fft(acc, nfft);
        freq = (0:nfft-1)' / (nfft * dt);
        
        % 频率响应函数
        omega_freq = 2 * pi * freq;
        H = (omega^2 + 2i*xi*omega*omega_freq) ./ ...
            (omega^2 - omega_freq.^2 + 2i*xi*omega*omega_freq);
        H(1) = 1;  % 直流分量
        
        % 计算响应
        resp_fft = acc_fft .* H;
        resp = real(ifft(resp_fft, nfft));
        resp = resp(1:N);
        
        % 计算最大绝对值
        Sa(i) = max(abs(resp));
    end
    
    % 确保没有NaN或Inf
    Sa(isnan(Sa)) = 0;
    Sa(isinf(Sa)) = 0;
end

function [acc_response, t_max] = calculate_sdof_response_simple(acc_input, dt, omega0, xi)
    % 简化的单自由度系统响应计算
    N = length(acc_input);
    
    % 使用频域方法计算响应
    nfft = 2^nextpow2(N);
    acc_fft = fft(acc_input, nfft);
    freq = (0:nfft-1)' / (nfft * dt);
    
    % 传递函数
    omega = 2 * pi * freq;
    H = (omega0^2 + 2i*xi*omega0*omega) ./ ...
        (omega0^2 - omega.^2 + 2i*xi*omega0*omega);
    H(1) = 1;
    
    % 计算响应
    resp_fft = acc_fft .* H;
    resp = real(ifft(resp_fft, nfft));
    acc_response = resp(1:N);
    
    % 找到最大绝对值时刻
    [~, idx_max] = max(abs(acc_response));
    t_max = (idx_max-1) * dt;
end

function [below_target, max_error, mean_error] = compute_error_stats_simple(Sa_calc, Sa_target)
    % 简化的误差统计计算
    Sa_calc = Sa_calc(:);
    Sa_target = Sa_target(:);
    
    % 避免除以0
    min_target = max(Sa_target) * 0.001;
    valid_idx = (Sa_target > min_target) & ~isnan(Sa_calc) & ~isinf(Sa_calc);
    
    if sum(valid_idx) < 3
        below_target = length(Sa_target);
        max_error = 1.0;
        mean_error = 1.0;
        return;
    end
    
    % 计算低于目标谱的点数
    below_target = sum(Sa_calc(valid_idx) < Sa_target(valid_idx));
    
    % 计算相对误差
    relative_error = abs(Sa_calc(valid_idx) - Sa_target(valid_idx)) ./ Sa_target(valid_idx);
    
    % 限制最大误差
    relative_error = min(relative_error, 1.0);
    
    max_error = max(relative_error);
    mean_error = mean(relative_error);
end

function composite_score = compute_composite_score(below_points, max_error, mean_error, envelope_ratio, params)
    % 计算综合评分（越低越好）
    % 反应谱匹配惩罚
    spectrum_penalty = 0;
    if below_points > params.max_below_points
        spectrum_penalty = spectrum_penalty + (below_points - params.max_below_points) * 5;
    end
    
    if max_error > params.error_tolerance
        spectrum_penalty = spectrum_penalty + (max_error - params.error_tolerance) * 50;
    end
    
    % 功率谱包络惩罚
    psd_penalty = 0;
    if envelope_ratio < params.psd_envelope_threshold
        psd_penalty = (params.psd_envelope_threshold - envelope_ratio) * 30;
    end
    
    % 基础评分（即使满足要求，也鼓励更小的误差）
    base_score = mean_error * 100 + below_points * 2 + (1 - envelope_ratio) * 20;
    
    % 综合评分
    composite_score = base_score + spectrum_penalty + psd_penalty;
end

function [freq, psd] = calculate_power_spectrum_simple(acc, dt)
    % 简化的功率谱密度计算方法
    N = length(acc);
    
    % 使用汉宁窗
    window = hanning(N);
    acc_windowed = acc .* window;
    
    % FFT计算
    nfft = 2^nextpow2(N);
    F = fft(acc_windowed, nfft);
    n_freq = floor(nfft/2) + 1;
    freq = (0:n_freq-1)' / (nfft * dt);
    
    % 计算功率谱密度
    amplitude = abs(F(1:n_freq)) * 2 / sum(window);
    psd = (amplitude.^2) / (2 * dt) * 980^2;  % 转换为cm²/s³
    
    % 处理直流分量
    if n_freq > 1
        psd(1) = psd(2);
    end
    
    % 确保没有NaN或Inf
    psd(isnan(psd)) = 0;
    psd(isinf(psd)) = 0;
end

function psd = calculate_RG1_60_power_spectrum_simple(freq)
    % 计算RG1.60标准功率谱（论文公式(1)）
    psd = zeros(size(freq));
    
    for i = 1:length(freq)
        f = freq(i);
        if f >= 0.3 && f < 2.5
            psd(i) = 4193.54 * (f/2.5)^0.2;
        elseif f >= 2.5 && f < 9.0
            psd(i) = 4193.54 * (2.5/f)^1.8;
        elseif f >= 9.0 && f < 16.0
            psd(i) = 418.06 * (9.0/f)^3;
        elseif f >= 16.0 && f <= 24.0
            psd(i) = 74.19 * (16.0/f)^8;
        else
            psd(i) = 0;
        end
    end
    
    % 确保没有NaN或Inf
    psd(isnan(psd)) = 0;
    psd(isinf(psd)) = 0;
end

function envelope_ratio = compute_psd_envelope_simple(freq_gen, psd_gen, freq_std, psd_std, params)
    % 简化的功率谱包络度计算
    % 在0.3-24Hz范围内检查是否包络80%的目标功率谱
    
    % 找到在指定频率范围内的点
    idx_std = freq_std >= params.psd_freq_range(1) & freq_std <= params.psd_freq_range(2);
    
    if sum(idx_std) == 0
        envelope_ratio = 0;
        return;
    end
    
    % 创建共同的频率点
    n_points = min(100, sum(idx_std));
    freq_common = linspace(params.psd_freq_range(1), params.psd_freq_range(2), n_points)';
    
    % 插值
    psd_gen_interp = interp1(freq_gen, psd_gen, freq_common, 'pchip', 0);
    psd_std_interp = interp1(freq_std(idx_std), psd_std(idx_std), freq_common, 'pchip', 0);
    
    % 确保没有负值
    psd_gen_interp = max(psd_gen_interp, 0);
    psd_std_interp = max(psd_std_interp, 0);
    
    % 计算阈值线（80%目标功率谱）
    threshold_line = params.psd_envelope_threshold * psd_std_interp;
    
    % 计算包络度
    valid_points = psd_std_interp > 0;
    if sum(valid_points) == 0
        envelope_ratio = 0;
        return;
    end
    
    envelope_points = psd_gen_interp(valid_points) >= threshold_line(valid_points);
    envelope_ratio = sum(envelope_points) / sum(valid_points);
end

function F = normal_search_adjustment(F, acc_current, freq_fft, T_target, f_target, Sa_target, params, relaxation_factor, tabu_list)
    % 正常搜索调整（按照论文方法）
    
    % 计算当前反应谱
    Sa_current = calculate_response_spectrum_simple(acc_current, params.dt, params.damping, T_target);
    
    n_freq = length(F);
    
    % 对每个反应谱控制频率进行调整
    for k = 1:length(f_target)
        f_k = f_target(k);
        Sa_target_k = Sa_target(k);
        Sa_current_k = Sa_current(k);
        
        % 如果当前点已满足要求，跳过
        if abs(Sa_current_k - Sa_target_k) / Sa_target_k < 0.01
            continue;
        end
        
        % 计算单自由度系统响应并找到最大响应时刻
        [acc_response, t_max] = calculate_sdof_response_simple(acc_current, params.dt, 2*pi*f_k, params.damping);
        
        % 找到最大响应值及其符号
        [~, idx_max] = max(abs(acc_response));
        a_max_R = acc_response(idx_max);
        sign_max = sign(a_max_R);
        
        % 计算有效带宽
        if k == 1
            f_k1 = 0;
            f_k2 = (f_target(k) + f_target(k+1)) / 2;
        elseif k == length(f_target)
            f_k1 = (f_target(k-1) + f_target(k)) / 2;
            f_k2 = f_target(k) * 1.2;
        else
            f_k1 = (f_target(k-1) + f_target(k)) / 2;
            f_k2 = (f_target(k) + f_target(k+1)) / 2;
        end
        
        % 确定调整方向
        if Sa_current_k < Sa_target_k
            % 需要增大响应
            if sign_max > 0
                adjustment_factor = 1 + relaxation_factor * (Sa_target_k/Sa_current_k - 1);
            else
                adjustment_factor = 1 - relaxation_factor * (1 - Sa_current_k/Sa_target_k);
            end
        else
            % 需要减小响应
            if sign_max > 0
                adjustment_factor = 1 - relaxation_factor * (1 - Sa_target_k/Sa_current_k);
            else
                adjustment_factor = 1 + relaxation_factor * (Sa_current_k/Sa_target_k - 1);
            end
        end
        
        % 限制调整幅度
        adjustment_factor = max(0.5, min(2.0, adjustment_factor));
        
        % 在有效带宽内调整傅里叶系数（避开禁忌的频率）
        for j = 2:n_freq-1
            f_j = freq_fft(j);
            
            % 检查是否在有效带宽内且不在禁忌表中
            if f_j >= f_k1 && f_j <= f_k2 && tabu_list(j) == 0
                % 计算该频率分量在t_max时刻的贡献
                A_j = abs(F(j));
                phi_j = angle(F(j));
                component_at_tmax = A_j * cos(2*pi*f_j*t_max + phi_j);
                
                % 判断贡献方向
                if component_at_tmax * sign_max > 0
                    % 与最大响应同方向
                    if Sa_current_k < Sa_target_k
                        F(j) = F(j) * adjustment_factor;
                    else
                        F(j) = F(j) / adjustment_factor;
                    end
                else
                    % 与最大响应反方向
                    if Sa_current_k < Sa_target_k
                        F(j) = F(j) / adjustment_factor;
                    else
                        F(j) = F(j) * adjustment_factor;
                    end
                end
            end
        end
    end
end

function F = diversification_search_adjustment(F, freq_fft, tabu_list, diversification_count)
    % 多样化搜索调整：随机调整一些频率分量
    
    n_freq = length(F);
    
    % 随机选择一些频率分量进行调整
    n_adjust = min(50, round(n_freq * 0.1));  % 调整10%的频率分量
    adjust_indices = randperm(n_freq-2, n_adjust) + 1;  % 跳过直流分量
    
    for idx = 1:n_adjust
        j = adjust_indices(idx);
        
        % 尽量避免调整禁忌的频率
        if tabu_list(j) == 0 || rand() < 0.3
            % 随机调整幅度
            random_factor = 0.7 + 0.6 * rand();  % 0.7~1.3
            F(j) = F(j) * random_factor;
            
            % 多样化搜索时可以调整相邻频率
            if j > 2 && j < n_freq-1
                neighbor_factor = 0.9 + 0.2 * rand();  % 0.9~1.1
                F(j-1) = F(j-1) * neighbor_factor;
                F(j+1) = F(j+1) * neighbor_factor;
            end
        end
    end
    
    % 多样化搜索时也可以调整相位
    if diversification_count > 2
        n_phase_adjust = min(20, round(n_freq * 0.05));
        phase_indices = randperm(n_freq-2, n_phase_adjust) + 1;
        
        for idx = 1:n_phase_adjust
            j = phase_indices(idx);
            if tabu_list(j) == 0
                % 随机调整相位
                phase_shift = 2 * pi * (rand() - 0.5) * 0.1;  % ±10%相位调整
                F(j) = F(j) * exp(1i * phase_shift);
            end
        end
    end
end

function F = intensification_search_adjustment(F, acc_current, freq_fft, T_target, f_target, Sa_target, Sa_current, params, relaxation_factor)
    % 强化搜索调整：针对当前误差最大的频率进行精细调整
    
    n_freq = length(F);
    
    % 找出误差最大的频率点
    errors = abs(Sa_current - Sa_target) ./ Sa_target;
    [~, max_error_idx] = max(errors);
    
    if max_error_idx > 0 && max_error_idx <= length(f_target)
        f_k = f_target(max_error_idx);
        Sa_target_k = Sa_target(max_error_idx);
        Sa_current_k = Sa_current(max_error_idx);
        
        % 计算单自由度系统响应
        [acc_response, t_max] = calculate_sdof_response_simple(acc_current, params.dt, 2*pi*f_k, params.damping);
        
        % 找到最大响应值及其符号
        [~, idx_max] = max(abs(acc_response));
        a_max_R = acc_response(idx_max);
        sign_max = sign(a_max_R);
        
        % 计算有效带宽（比正常搜索更窄）
        if max_error_idx == 1
            f_k1 = 0;
            f_k2 = (f_target(max_error_idx) + f_target(max_error_idx+1)) / 2;
        elseif max_error_idx == length(f_target)
            f_k1 = (f_target(max_error_idx-1) + f_target(max_error_idx)) / 2;
            f_k2 = f_target(max_error_idx);
        else
            f_k1 = (f_target(max_error_idx-1) + f_target(max_error_idx)) / 2;
            f_k2 = (f_target(max_error_idx) + f_target(max_error_idx+1)) / 2;
        end
        
        % 使用较小的调整因子进行精细调整
        fine_relaxation = relaxation_factor * 0.5;
        
        % 确定调整方向
        if Sa_current_k < Sa_target_k
            % 需要增大响应
            if sign_max > 0
                adjustment_factor = 1 + fine_relaxation * (Sa_target_k/Sa_current_k - 1);
            else
                adjustment_factor = 1 - fine_relaxation * (1 - Sa_current_k/Sa_target_k);
            end
        else
            % 需要减小响应
            if sign_max > 0
                adjustment_factor = 1 - fine_relaxation * (1 - Sa_target_k/Sa_current_k);
            else
                adjustment_factor = 1 + fine_relaxation * (Sa_current_k/Sa_target_k - 1);
            end
        end
        
        % 限制调整幅度
        adjustment_factor = max(0.9, min(1.1, adjustment_factor));
        
        % 在有效带宽内精细调整傅里叶系数
        for j = 2:n_freq-1
            f_j = freq_fft(j);
            
            % 检查是否在有效带宽内
            if f_j >= f_k1 && f_j <= f_k2
                % 计算该频率分量在t_max时刻的贡献
                A_j = abs(F(j));
                phi_j = angle(F(j));
                component_at_tmax = A_j * cos(2*pi*f_j*t_max + phi_j);
                
                % 判断贡献方向
                if component_at_tmax * sign_max > 0.5
                    % 与最大响应同方向，进行精细调整
                    if Sa_current_k < Sa_target_k
                        F(j) = F(j) * adjustment_factor;
                    else
                        F(j) = F(j) / adjustment_factor;
                    end
                elseif component_at_tmax * sign_max < -0.5
                    % 与最大响应反方向，进行精细调整
                    if Sa_current_k < Sa_target_k
                        F(j) = F(j) / adjustment_factor;
                    else
                        F(j) = F(j) * adjustment_factor;
                    end
                end
            end
        end
    end
end

function status_str = get_status_string(status_code)
    % 获取状态字符串
    switch status_code
        case 0
            status_str = '初始';
        case 1
            status_str = '正常';
        case 2
            status_str = '多样化';
        case 3
            status_str = '强化';
        otherwise
            status_str = '未知';
    end
end

function create_results_plots_tabu(time, acc, vel, disp, ...
    T_target, Sa_target, Sa_final, ...
    freq_gen, psd_gen, freq_std, psd_std, ...
    iteration_history, iter, ...
    params, PGA, PGV, PGD, error_final, envelope_ratio, best_iteration)
    
    try
        fig = figure('Position', [100, 100, 1400, 900], 'Name', '禁忌搜索改进的结果', 'Color', 'w');
        
        % 1. 反应谱对比
        subplot(3, 3, [1, 2]);
        hold on; grid on; box on;
        plot(T_target, Sa_target, 'b-', 'LineWidth', 2.5, 'DisplayName', '目标反应谱');
        plot(T_target, Sa_final, 'r-', 'LineWidth', 1.8, 'DisplayName', sprintf('生成波 (平均误差: %.2f%%)', error_final*100));
        
        % 标记低于目标谱的点
        below_idx = Sa_final < Sa_target;
        if any(below_idx)
            plot(T_target(below_idx), Sa_final(below_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2, ...
                'MarkerFaceColor', 'r', 'DisplayName', sprintf('低于目标谱点 (%d个)', sum(below_idx)));
        end
        
        % 添加误差带
        error_band = Sa_target * params.error_tolerance;
        fill([T_target; flipud(T_target)], [Sa_target.*(1-params.error_tolerance); flipud(Sa_target.*(1+params.error_tolerance))], ...
            [0.8 0.9 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', sprintf('±%.0f%%误差带', params.error_tolerance*100));
        
        xlabel('周期 T (s)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('谱加速度 Sa (g)', 'FontSize', 12, 'FontWeight', 'bold');
        title('反应谱匹配结果', 'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legend('Location', 'best', 'FontSize', 10);
        xlim([min(T_target), max(T_target)]);
        
        % 2. 功率谱对比
        subplot(3, 3, 3);
        hold on; grid on; box on;
        
        % 填充包络要求区域
        idx_fill = freq_std >= params.psd_freq_range(1) & freq_std <= params.psd_freq_range(2);
        freq_fill = freq_std(idx_fill);
        psd_fill = psd_std(idx_fill) * params.psd_envelope_threshold;
        fill([freq_fill; flipud(freq_fill)], [psd_fill; zeros(size(psd_fill))], ...
            [0.8 0.95 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', '包络要求区域');
        
        plot(freq_std, psd_std, 'b-', 'LineWidth', 2, 'DisplayName', '目标功率谱');
        plot(freq_gen, psd_gen, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('生成波 (包络度: %.1f%%)', envelope_ratio*100));
        plot(freq_std, psd_std * params.psd_envelope_threshold, 'b--', 'LineWidth', 1.2, 'DisplayName', '80%包络线');
        
        xlabel('频率 (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('功率谱密度 (cm²/s³)', 'FontSize', 12, 'FontWeight', 'bold');
        title('功率谱对比与包络验证', 'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legend('Location', 'best', 'FontSize', 9);
        xlim([0.1, 30]);
        ylim([1e1, 1e5]);
        
        % 3. 迭代历史 - 综合评分
        subplot(3, 3, [4, 5]);
        valid_iter = iteration_history(1:iter, 1) > 0;
        iters = iteration_history(valid_iter, 1);
        
        plot(iters, iteration_history(valid_iter, 6), 'k-', 'LineWidth', 2, 'DisplayName', '综合评分');
        hold on; grid on; box on;
        
        % 标记最佳解位置
        plot(best_iteration, iteration_history(best_iteration, 6), 'r*', 'MarkerSize', 15, ...
            'LineWidth', 2, 'DisplayName', sprintf('最佳解 (第%d次迭代)', best_iteration));
        
        % 标记搜索策略变化
        strategy_colors = [0 0 0; 0 0 1; 1 0.5 0; 0 0.5 0]; % 黑、蓝、橙、绿
        for i = 2:length(iters)
            strategy = iteration_history(i, 7);
            if strategy > 0 && strategy <= 4
                plot([iters(i-1), iters(i)], [iteration_history(i-1, 6), iteration_history(i, 6)], ...
                    'Color', strategy_colors(strategy, :), 'LineWidth', 1.5);
            end
        end
        
        xlabel('迭代次数', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('综合评分 (越低越好)', 'FontSize', 12, 'FontWeight', 'bold');
        title('迭代历史 - 综合评分', 'FontSize', 14, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        
        % 4. 迭代历史 - 各项指标
        subplot(3, 3, [6, 7]);
        yyaxis left
        plot(iters, iteration_history(valid_iter, 2), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        hold on; grid on; box on;
        plot([1, max(iters)], [params.max_below_points, params.max_below_points], 'r--', 'LineWidth', 1.5);
        ylabel('低于目标谱点数', 'FontSize', 12, 'FontWeight', 'bold');
        
        yyaxis right
        plot(iters, iteration_history(valid_iter, 5), 'g^-', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'g');
        plot([1, max(iters)], [params.psd_envelope_threshold*100, params.psd_envelope_threshold*100], 'g--', 'LineWidth', 1.5);
        ylabel('功率谱包络度 (%)', 'FontSize', 12, 'FontWeight', 'bold');
        
        xlabel('迭代次数', 'FontSize', 12, 'FontWeight', 'bold');
        title('迭代历史 - 关键指标', 'FontSize', 14, 'FontWeight', 'bold');
        legend({'低于谱点数', '允许最大值', '包络度', '要求值'}, 'Location', 'best', 'FontSize', 9);
        
        % 5. 加速度时程
        subplot(3, 3, 8);
        plot(time, acc, 'b-', 'LineWidth', 1.2);
        hold on; grid on; box on;
        
        % 绘制包络线
        envelope = get_envelope_vector_simple(time, params);
        plot(time, params.PGA_target * envelope, 'r--', 'LineWidth', 1.5);
        plot(time, -params.PGA_target * envelope, 'r--', 'LineWidth', 1.5);
        
        xlabel('时间 (s)', 'FontSize', 12);
        ylabel('加速度 (g)', 'FontSize', 12);
        title(sprintf('加速度时程\nPGA = %.4f g', PGA), 'FontSize', 13, 'FontWeight', 'bold');
        xlim([0, min(20, time(end))]);
        legend({'加速度时程', '包络线'}, 'Location', 'best', 'FontSize', 9);
        
        % 6. 速度和位移时程
        subplot(3, 3, 9);
        yyaxis left
        plot(time, vel*100, 'r-', 'LineWidth', 1.2);
        ylabel('速度 (cm/s)', 'FontSize', 12, 'Color', 'r');
        ylim([-max(abs(vel*100))*1.2, max(abs(vel*100))*1.2]);
        
        yyaxis right
        plot(time, disp*100, 'm-', 'LineWidth', 1.2);
        ylabel('位移 (cm)', 'FontSize', 12, 'Color', 'm');
        ylim([-max(abs(disp*100))*1.2, max(abs(disp*100))*1.2]);
        
        grid on; box on;
        xlabel('时间 (s)', 'FontSize', 12);
        title(sprintf('速度和位移时程\nPGV = %.1f cm/s, PGD = %.1f cm', PGV, PGD), ...
            'FontSize', 13, 'FontWeight', 'bold');
        xlim([0, min(20, time(end))]);
        legend({'速度', '位移'}, 'Location', 'best', 'FontSize', 9);
        
        % 添加总体信息文本
        annotation('textbox', [0.02, 0.02, 0.96, 0.05], ...
            'String', sprintf('禁忌搜索改进算法 | 总迭代: %d | 最佳迭代: %d | 平均误差: %.2f%% | 包络度: %.1f%% | PGA: %.4fg', ...
            iter, best_iteration, error_final*100, envelope_ratio*100, PGA), ...
            'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'BackgroundColor', [0.9, 0.95, 1.0], ...
            'EdgeColor', 'b', 'LineWidth', 1);
        
        % 保存图形
        set(fig, 'PaperPositionMode', 'auto');
        saveas(fig, '禁忌搜索改进算法结果.png');
        fprintf('结果图已保存: 禁忌搜索改进算法结果.png\n');
        
    catch ME
        fprintf('创建图形时出错: %s\n', ME.message);
    end
end

function save_results_files_tabu(time, acc, vel, disp, ...
    T_target, Sa_target, Sa_final, params, ...
    error_final, envelope_ratio, PGA, PGV, PGD, ...
    iteration_history, best_iteration)
    
    try
        % 1. 时程数据
        time_data = [time, acc, vel*100, disp*100];
        fid = fopen('禁忌搜索法地震波时程数据.txt', 'w');
        fprintf(fid, '时间(s)\t加速度(g)\t速度(cm/s)\t位移(cm)\n');
        for i = 1:min(3000, size(time_data, 1))
            fprintf(fid, '%.6f\t%.6e\t%.6e\t%.6e\n', time_data(i,:));
        end
        fclose(fid);
        
        % 2. 反应谱数据
        fid = fopen('禁忌搜索法反应谱数据.txt', 'w');
        fprintf(fid, '周期(s)\t目标谱(g)\t生成谱(g)\t相对误差(%%)\n');
        for i = 1:length(T_target)
            if Sa_target(i) > 0
                error = abs(Sa_final(i) - Sa_target(i)) / Sa_target(i) * 100;
            else
                error = 0;
            end
            fprintf(fid, '%.6f\t%.6e\t%.6e\t%.2f\n', ...
                T_target(i), Sa_target(i), Sa_final(i), error);
        end
        fclose(fid);
        
        % 3. 迭代历史数据
        valid_iter = iteration_history(:,1) > 0;
        iter_data = iteration_history(valid_iter, :);
        fid = fopen('禁忌搜索法迭代历史.txt', 'w');
        fprintf(fid, '迭代次数\t低于谱点数\t最大误差%%\t平均误差%%\t包络度%%\t综合评分\t搜索策略\n');
        for i = 1:size(iter_data, 1)
            strategy_str = get_status_string(iter_data(i, 7));
            fprintf(fid, '%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n', ...
                iter_data(i,1), iter_data(i,2), iter_data(i,3), iter_data(i,4), ...
                iter_data(i,5), iter_data(i,6), strategy_str);
        end
        fclose(fid);
        
        % 4. 参数文件
        fid = fopen('禁忌搜索法参数配置.txt', 'w');
        fprintf(fid, '参数配置（禁忌搜索改进）\n');
        fprintf(fid, '================================\n');
        fprintf(fid, '总时长: %.0f s\n', params.T_total);
        fprintf(fid, '上升段: %.1f s\n', params.t1);
        fprintf(fid, '平稳段: %.1f s\n', params.t2-params.t1);
        fprintf(fid, '阻尼比: %.2f\n', params.damping);
        fprintf(fid, '目标PGA: %.3f g\n', params.PGA_target);
        fprintf(fid, '最大迭代次数: %d\n', params.max_iterations);
        fprintf(fid, '允许最大低于谱点数: %d\n', params.max_below_points);
        fprintf(fid, '允许最大误差: %.1f%%\n', params.error_tolerance*100);
        fprintf(fid, '功率谱包络目标: ≥%.1f%%\n', params.psd_envelope_threshold*100);
        fprintf(fid, '禁忌长度: %d\n', params.tabu_tenure);
        fprintf(fid, '多样化阈值: %d\n', params.diversification_threshold);
        fprintf(fid, '强化阈值: %d\n', params.intensification_threshold);
        fprintf(fid, '自适应松弛因子: %s\n', mat2str(params.adaptive_relaxation));
        fclose(fid);
        
        % 5. 生成报告
        fid = fopen('禁忌搜索改进算法报告.txt', 'w');
        fprintf(fid, '禁忌搜索改进的人工地震动时程生成报告\n');
        fprintf(fid, '==========================================\n\n');
        
        fprintf(fid, '参考论文:\n');
        fprintf(fid, '王海涛，何树延，张征明，于溯源\n');
        fprintf(fid, '"同时满足反应谱匹配和功率谱密度包络要求的人工时程生成算法"\n\n');
        
        fprintf(fid, '生成时间: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        
        fprintf(fid, '一、算法参数\n');
        fprintf(fid, '--------------------------------\n');
        fprintf(fid, '总时长: %.0f s\n', params.T_total);
        fprintf(fid, '上升段: %.1f s\n', params.t1);
        fprintf(fid, '平稳段: %.1f s\n', params.t2-params.t1);
        fprintf(fid, '阻尼比: %.2f\n', params.damping);
        fprintf(fid, '目标PGA: %.3f g\n', params.PGA_target);
        fprintf(fid, '最大迭代次数: %d\n', params.max_iterations);
        fprintf(fid, '允许最大低于谱点数: %d（论文要求）\n', params.max_below_points);
        fprintf(fid, '允许最大误差: %.1f%%（论文要求）\n', params.error_tolerance*100);
        fprintf(fid, '功率谱包络目标: ≥%.1f%%（论文要求）\n', params.psd_envelope_threshold*100);
        fprintf(fid, '禁忌长度: %d\n', params.tabu_tenure);
        fprintf(fid, '多样化阈值: %d\n', params.diversification_threshold);
        fprintf(fid, '强化阈值: %d\n', params.intensification_threshold);
        fprintf(fid, '自适应松弛因子: %s\n\n', mat2str(params.adaptive_relaxation));
        
        fprintf(fid, '二、最终结果\n');
        fprintf(fid, '--------------------------------\n');
        fprintf(fid, '反应谱平均误差: %.3f%%\n', error_final*100);
        
        % 计算最大误差
        valid_idx = Sa_target > 0;
        if any(valid_idx)
            max_rel_error = max(abs(Sa_final(valid_idx)-Sa_target(valid_idx))./Sa_target(valid_idx))*100;
            fprintf(fid, '反应谱最大误差: %.3f%%\n', max_rel_error);
        end
        
        fprintf(fid, '低于目标谱点数: %d/%d\n', sum(Sa_final < Sa_target), length(T_target));
        fprintf(fid, '功率谱包络度: %.1f%%\n', envelope_ratio*100);
        fprintf(fid, '峰值加速度 (PGA): %.4f g\n', PGA);
        fprintf(fid, '峰值速度 (PGV): %.2f cm/s\n', PGV);
        fprintf(fid, '峰值位移 (PGD): %.2f cm\n\n', PGD);
        
        fprintf(fid, '三、论文要求评估\n');
        fprintf(fid, '--------------------------------\n');
        if sum(Sa_final < Sa_target) <= params.max_below_points
            fprintf(fid, '低于目标谱点数: 满足论文要求 ✓ (要求≤%d)\n', params.max_below_points);
        else
            fprintf(fid, '低于目标谱点数: 未满足论文要求 ✗ (要求≤%d)\n', params.max_below_points);
        end
        
        if max_rel_error <= params.error_tolerance*100
            fprintf(fid, '最大相对误差: 满足论文要求 ✓ (要求≤%.1f%%)\n', params.error_tolerance*100);
        else
            fprintf(fid, '最大相对误差: 未满足论文要求 ✗ (要求≤%.1f%%)\n', params.error_tolerance*100);
        end
        
        if envelope_ratio >= params.psd_envelope_threshold
            fprintf(fid, '功率谱包络度: 满足论文要求 ✓ (要求≥%.1f%%)\n', params.psd_envelope_threshold*100);
        else
            fprintf(fid, '功率谱包络度: 未满足论文要求 ✗ (要求≥%.1f%%)\n', params.psd_envelope_threshold*100);
        end
        
        fprintf(fid, '\n四、搜索过程统计\n');
        fprintf(fid, '--------------------------------\n');
        fprintf(fid, '总迭代次数: %d\n', sum(valid_iter));
        fprintf(fid, '最佳解在第 %d 次迭代获得\n', best_iteration);
        
        % 统计搜索策略使用情况
        strategies = iter_data(:, 7);
        fprintf(fid, '正常搜索次数: %d\n', sum(strategies == 1));
        fprintf(fid, '多样化搜索次数: %d\n', sum(strategies == 2));
        fprintf(fid, '强化搜索次数: %d\n', sum(strategies == 3));
        
        fprintf(fid, '\n五、生成文件清单\n');
        fprintf(fid, '--------------------------------\n');
        fprintf(fid, '1. 禁忌搜索法地震波时程数据.txt\n');
        fprintf(fid, '2. 禁忌搜索法反应谱数据.txt\n');
        fprintf(fid, '3. 禁忌搜索法迭代历史.txt\n');
        fprintf(fid, '4. 禁忌搜索法参数配置.txt\n');
        fprintf(fid, '5. 禁忌搜索改进算法报告.txt\n');
        fprintf(fid, '6. 禁忌搜索改进算法结果.png\n');
        
        fclose(fid);
        
        fprintf('结果文件已保存:\n');
        fprintf('  时程数据: 禁忌搜索法地震波时程数据.txt\n');
        fprintf('  反应谱数据: 禁忌搜索法反应谱数据.txt\n');
        fprintf('  迭代历史: 禁忌搜索法迭代历史.txt\n');
        fprintf('  参数配置: 禁忌搜索法参数配置.txt\n');
        fprintf('  详细报告: 禁忌搜索改进算法报告.txt\n');
        fprintf('  图形文件: 禁忌搜索改进算法结果.png\n');
        
    catch ME
        fprintf('保存结果时出错: %s\n', ME.message);
    end
end