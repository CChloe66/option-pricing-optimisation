#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Utility function to return a standard normal random variable using the Box-Muller transform
double gaussian_rng() {
    double u = (double)rand() / RAND_MAX;
    double v = (double)rand() / RAND_MAX;
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

// Computes the Asian option price using Monte Carlo simulation
double asian_option_monte_carlo(int num_paths, double S0, double K, double T, double r, double sigma, int num_points) {
    double dt = T / num_points; // Time increment per simulation step
    double drift = (r - 0.5 * sigma * sigma) * dt; // Pre-computed drift term of the stock price path formula
    double vol = sigma * sqrt(dt); // Pre-computed volatility term of the stock price path formula
    double discount_factor = exp(-r * T); // Pre-computed discount factor applied to final payoff, computed outside the loop

    double total_payoff = 0.0; // Initialize total payoff accumulated over all paths

    for (int i = 0; i < num_paths; i++) {
        double path_sum = 0.0; // Sum of stock prices along the path to compute average
        double S = S0; // Initial stock price at the beginning of the path
        for (int j = 0; j < num_points; j++) {
            double z = gaussian_rng(); // Generate a standard normal random variable for stock price simulation
            S *= exp(drift + vol * z); // Update the stock price using pre-computed terms
            path_sum += S; // Accumulate the sum of stock prices for averaging
        }
        total_payoff += fmax(path_sum / num_points - K, 0); // Calculate the payoff for the path and add to total payoff
    }

    return discount_factor * total_payoff / num_paths; // Compute the average discounted payoff over all paths to estimate the option price
}

int main() {
    clock_t start, end;
    double cpu_time_used;

    // Seed the random number generator
    srand(time(NULL));

    // Start the clock to measure execution time
    start = clock();

    // Option parameters
    int num_paths = 10000000;
    int num_points = 4096;
    double S0 = 100.0; // Initial stock price
    double K = 100.0; // Strike price
    double T = 1.0; // Time to maturity in years
    double r = 0.05; // Risk-free interest rate
    double sigma = 0.2; // Volatility

    double price = asian_option_monte_carlo(num_paths, S0, K, T, r, sigma, num_points);
    
    // End the clock and calculate execution time
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Estimated Asian option price using Monte Carlo simulation: $%.2f\n", price);
    printf("Execution time: %f seconds\n", cpu_time_used);

    return 0;
}
