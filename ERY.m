% 参数设置
S0 = 1000; % 初始敏感细菌丰度
R0 = 100;  % 初始耐药细菌丰度
antibiotic_concentration = 0:20:2000; % 抗生素浓度范围
time = 0:0.01:100; % 时间范围，0到100，步长为0.01

% 非饮用水源地参数
M_non_drink = 5; % 营养浓度
mu_Rmax_non_drink = 0.11; % 耐药细菌最大生长速率
mu_Rmin_non_drink = -0.1; % 耐药细菌最小生长速率
mu_Smax_non_drink = 0.14; % 敏感细菌最大生长速率
mu_Smin_non_drink = -0.1; % 敏感细菌最小生长速率
d_non_drink = 0.6; % 死亡率
N = 100000; % 环境容纳量
k_non_drink = 150; % 希尔系数
KM_non_drink = 200; % 营养系数
MICs_non_drink = 4000; % 敏感细菌希尔系数
MICr_non_drink = 8000; % 耐药细菌希尔系数

% 饮用水源地参数
M_drink = 2; % 营养浓度
mu_Rmax_drink = 0.08;
mu_Rmin_drink = -0.1;
mu_Smax_drink = 0.16;
mu_Smin_drink = -0.1;
d_drink = 0.5;
k_drink = 10;
KM_drink = 100;
MICs_drink = 1000;  
MICr_drink = 15000;  

% 定义耐药率计算函数
calculate_resistance_rate = @(R, S, R0, S0) (R ./ (R + S)) ./ (R0 / (R0 + S0));

% 初始化曲线数据
resistance_rate_non_drink = zeros(length(antibiotic_concentration), length(time));
resistance_rate_drink = zeros(length(antibiotic_concentration), length(time));

% 模拟不同抗生素浓度下的细菌动态
for i = 1:length(antibiotic_concentration)
    a = antibiotic_concentration(i);
    
    % 初始化细菌丰度
    R_non_drink = R0;
    S_non_drink = S0;
    R_drink = R0;
    S_drink = S0;
    
    % 迭代时间步长
    for t = 1:length(time)
        % 计算非饮用水源地的耐药细菌增长速率
        growth_rate_R_non_drink = mu_Rmax_non_drink * M_non_drink / (KM_non_drink + M_non_drink) - ...
            ((mu_Rmax_non_drink * M_non_drink / (KM_non_drink + M_non_drink) - mu_Rmin_non_drink) * (a / MICr)^k_non_drink) / ...
            ((a / MICr)^k_non_drink - mu_Rmin_non_drink / mu_Rmax_non_drink);
        
        % 计算非饮用水源地的敏感细菌增长速率
        growth_rate_S_non_drink = mu_Smax_non_drink * M_non_drink / (KM_non_drink + M_non_drink) - ...
            ((mu_Smax_non_drink * M_non_drink / (KM_non_drink + M_non_drink) - mu_Smin_non_drink) * (a / MICs)^k_non_drink) / ...
            ((a / MICs)^k_non_drink - mu_Smin_non_drink / mu_Smax_non_drink);
        
        % 更新非饮用水源地的细菌丰度
        R_non_drink = R_non_drink + (growth_rate_R_non_drink * R_non_drink * (1 - (R_non_drink + S_non_drink) / N) - d_non_drink * R_non_drink) * 0.01;
        S_non_drink = S_non_drink + (growth_rate_S_non_drink * S_non_drink * (1 - (R_non_drink + S_non_drink) / N) - d_non_drink * S_non_drink) * 0.01;
        
        % 计算非饮用水源地的耐药率
        resistance_rate_non_drink(i, t) = calculate_resistance_rate(R_non_drink, S_non_drink, R0, S0);
        
        % 计算饮用水源地的耐药细菌增长速率
        growth_rate_R_drink = mu_Rmax_drink * M_drink / (KM_drink + M_drink) - ...
            ((mu_Rmax_drink * M_drink / (KM_drink + M_drink) - mu_Rmin_drink) * (a / MICr)^k_drink) / ...
            ((a / MICr)^k_drink - mu_Rmin_drink / mu_Rmax_drink);
        
        % 计算饮用水源地的敏感细菌增长速率
        growth_rate_S_drink = mu_Smax_drink * M_drink / (KM_drink + M_drink) - ...
            ((mu_Smax_drink * M_drink / (KM_drink + M_drink) - mu_Smin_drink) * (a / MICs)^k_drink) / ...
            ((a / MICs)^k_drink - mu_Smin_drink / mu_Smax_drink);
        
        % 更新饮用水源地的细菌丰度
        R_drink = R_drink + (growth_rate_R_drink * R_drink * (1 - (R_drink + S_drink) / N) - d_drink * R_drink) * 0.01;
        S_drink = S_drink + (growth_rate_S_drink * S_drink * (1 - (R_drink + S_drink) / N) - d_drink * S_drink) * 0.01;
        
        % 计算饮用水源地的耐药率
        resistance_rate_drink(i, t) = calculate_resistance_rate(R_drink, S_drink, R0, S0);
    end
end

% 绘制结果，取最终时间点的数据
figure;
plot(antibiotic_concentration, resistance_rate_non_drink(:, end), 'b', 'LineWidth', 2);
hold on;
plot(antibiotic_concentration, resistance_rate_drink(:, end), 'r', 'LineWidth', 2);
xlabel('抗生素浓度');
ylabel('耐药率');
title('ERY的耐药率与抗生素浓度图');
legend('非饮用水源地', '饮用水源地');
grid on;
hold off;
