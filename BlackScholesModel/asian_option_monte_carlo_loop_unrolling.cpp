#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

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
        for (int j = 0; j < num_points; j += 2) { // Loop unrolling
            double z1 = gaussian_rng(); // Generate a standard normal random variable for stock price simulation
            S *= exp(drift + vol * z1); // Update the stock price using pre-computed terms
            path_sum += S;

            if (j + 1 < num_points) { // Check to prevent out-of-bounds access
                double z2 = gaussian_rng(); // Generate another random variable for the next step
                S *= exp(drift + vol * z2); // Update the stock price using pre-computed terms
                path_sum += S;
            }
        }
        total_payoff += fmax(path_sum / num_points - K, 0); // Calculate the payoff for the path and add to total payoff
    }

    return discount_factor * total_payoff / num_paths; // Compute the average discounted payoff over all paths to estimate the option price
}

int main() {
    srand(time(NULL)); // Seed the random number generator

    struct timeval start_wall, end_wall; // Variables to hold wall clock times
    gettimeofday(&start_wall, NULL); // Start the wall clock timing

    clock_t start_cpu = clock(); // Start the CPU timing

    int num_paths = 10000000;
    int num_points = 4096;
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;

    double price = asian_option_monte_carlo(num_paths, S0, K, T, r, sigma, num_points);

    clock_t end_cpu = clock(); // End the CPU timing
    gettimeofday(&end_wall, NULL); // End the wall clock timing

    double cpu_time_used = ((double) (end_cpu - start_cpu)) / CLOCKS_PER_SEC;
    double wall_time_used = (end_wall.tv_sec - start_wall.tv_sec) + (end_wall.tv_usec - start_wall.tv_usec) / 1000000.0;

    printf("Estimated Asian option price using Monte Carlo simulation: $%.2f\n", price);
    printf("CPU execution time: %f seconds\n", cpu_time_used);
    printf("Wall clock execution time: %f seconds\n", wall_time_used);

    return 0;
}
