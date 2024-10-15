#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

// Computes the Asian option price using Binomial model simulation
double asian_option_binomial(int num_paths, int num_points, double S0, double K, double T, double r, double sigma) {
    double dt = T / num_points;                                // 每步的时间
    double upFactor = exp(sigma * sqrt(dt));                   // 上升因子
    double downFactor = 1.0 / upFactor;                        // 下降因子
    double riskNeutralProb = (exp(r * dt) - downFactor) / (upFactor - downFactor);  // 风险中性概率
    double discountFactor = exp(-r * T);                       // 折现因子
    double total_payoff = 0.0;

    for (int path = 0; path < num_paths; path++) {
        double cumulative_price = S0;                          // 初始化累计价格
        double stock_price = S0;                               // 初始股票价格

        // 模拟单条路径
        for (int step = 1; step <= num_points; step++) {
            double rand_num = (double)rand() / RAND_MAX;       // 生成随机数
            if (rand_num < riskNeutralProb) {
                stock_price *= upFactor;                       // 股票价格上涨
            } else {
                stock_price *= downFactor;                     // 股票价格下跌
            }
            cumulative_price += stock_price;                   // 累加每步的股票价格
        }

        // 计算该路径的平均价格
        double average_price = cumulative_price / (num_points + 1);
        // 计算该路径的期权收益
        double payoff = fmax(average_price - K, 0.0);          // 只计算看涨期权（call）
        total_payoff += payoff;
    }

    // 返回期权价格的折现期望值
    return discountFactor * total_payoff / num_paths;
}

int main() {
    struct timeval start_time, end_time;
    clock_t start_cpu, end_cpu;

    gettimeofday(&start_time, NULL);  // Start the wall clock timing
    start_cpu = clock();              // Start the CPU timing

    // Seed the random number generator
    srand(time(NULL));

    // Option parameters
    int num_paths = 10000000;
    int num_points = 4096;
    double S0 = 100.0;    // 初始股票价格
    double K = 100.0;     // 行权价格
    double T = 1.0;       // 到期时间（年）
    double r = 0.05;      // 无风险利率
    double sigma = 0.2;   // 波动率

    // 计算亚洲期权价格
    double price = asian_option_binomial(num_paths, num_points, S0, K, T, r, sigma);

    end_cpu = clock();                // End the CPU timing
    gettimeofday(&end_time, NULL);    // End the wall clock timing

    double cpu_time_used = ((double) (end_cpu - start_cpu)) / CLOCKS_PER_SEC;
    double elapsed_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    // 输出结果
    printf("Estimated Asian option price using Binomial model: $%.2f\n", price);
    printf("CPU execution time: %f seconds\n", cpu_time_used);
    printf("Wall clock execution time: %f seconds\n", elapsed_time);

    return 0;
}