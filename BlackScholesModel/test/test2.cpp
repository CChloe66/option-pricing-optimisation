#include <iostream>
#include <vector>
#include <cmath>
#include <random>

// 生成符合正态分布的随机数
double normalRandom(double mean = 0.0, double std_dev = 1.0) {
    static std::random_device rd;
    static std::mt19937 generator(rd());
    std::normal_distribution<double> distribution(mean, std_dev);
    return distribution(generator);
}

// Andersen QE方法简化版
double andersenQE(double v0, double kappa, double theta, double sigma, double rho, double dt) {
    double psi = sigma * sigma * dt / (4.0 * kappa);
    double m = theta + (v0 - theta) * std::exp(-kappa * dt);
    double s2 = v0 * sigma * sigma * std::exp(-kappa * dt) / kappa * (1.0 - std::exp(-kappa * dt)) + theta * sigma * sigma / (2.0 * kappa) * (1.0 - std::exp(-kappa * dt)) * (1.0 - std::exp(-kappa * dt));
    double phi = s2 / (m * m);

    double u = normalRandom();
    double A = psi * (1.0 - std::exp(-kappa * dt)) / (4.0 * kappa);
    double B = 1.0 + phi;
    double z = std::log(u) / A;
    double v = m * (1.0 + z + std::sqrt(z * z + 2.0 * z * B));

    return v;
}

// 假设的支付函数，根据资产价格计算期权价值
double payoff(double assetPrice) {
    return std::max(assetPrice - 100.0, 0.0);  // 看涨期权的支付函数
}

// Longstaff-Schwartz简化版（只演示蒙特卡洛路径生成，未包含回归部分）
std::vector<double> longstaffSchwartz(int numPaths, int numTimesteps, double S0, double dt) {
    std::vector<double> optionValues(numPaths, 0.0);
    for (int i = 0; i < numPaths; ++i) {
        double assetPrice = S0;
        for (int t = 0; t < numTimesteps; ++t) {
            assetPrice += assetPrice * normalRandom() * std::sqrt(dt);  // 简化的资产价格模拟
        }
        optionValues[i] = payoff(assetPrice);  // 计算期权的终值
    }
    return optionValues;
}

