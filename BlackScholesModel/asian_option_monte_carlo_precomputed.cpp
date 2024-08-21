#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Utility function to return a standard normal random variable (Box-Muller transform)
double gaussian_rng() {
    double u = (double)rand() / RAND_MAX;
    double v = (double)rand() / RAND_MAX;
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

// Computes the Asian option price using Monte Carlo simulation
double asian_option_monte_carlo(int num_paths, double S0, double K, double T, double r, double sigma, int num_points) {
    double dt = T / num_points;
    double drift = (r - 0.5 * sigma * sigma) * dt;   // Precomputed part of the exponent for S update
    double vol = sigma * sqrt(dt);                   // Precomputed volatility term
    double total_payoff = 0.0;

    for (int i = 0; i < num_paths; i++) {
        double path_sum = 0.0;
        double S = S0;
        for (int j = 0; j < num_points; j++) {
            double z = gaussian_rng();
            S *= exp(drift + vol * z);               // Using precomputed values
            path_sum += S;
        }
        double average_price = path_sum / num_points;
        double payoff = fmax(average_price - K, 0);
        total_payoff += payoff;
    }

    double price = exp(-r * T) * total_payoff / num_paths;
    return price;
}

int main() {
    clock_t start, end;
    double cpu_time_used;

    // Seed the random number generator
    srand(time(NULL));

    // Start the clock
    start = clock();

    // Option parameters
    int num_paths = 10000000;
    int num_points = 4096;
    double S0 = 100.0;    // Initial stock price
    double K = 100.0;     // Strike price
    double T = 1.0;       // Time to maturity in years
    double r = 0.05;      // Risk-free interest rate
    double sigma = 0.2;   // Volatility

    double price = asian_option_monte_carlo(num_paths, S0, K, T, r, sigma, num_points);
    
    // End the clock
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Estimated Asian option price using Monte Carlo simulation: $%.2f\n", price);
    printf("Execution time: %f seconds\n", cpu_time_used);

    return 0;
}
