#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>

// Utility function to return a standard normal random variable (Box-Muller transform)
double gaussian_rng(unsigned int *seed) {
    double u = (double)rand_r(seed) / RAND_MAX;
    double v = (double)rand_r(seed) / RAND_MAX;
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

// Computes the Asian option price using Monte Carlo simulation
double asian_option_monte_carlo(int num_paths, double S0, double K, double T, double r, double sigma, int num_points) {
    double dt = T / num_points;
    double drift = (r - 0.5 * sigma * sigma) * dt;
    double vol = sigma * sqrt(dt);
    double discount_factor = exp(-r * T);
    double total_payoff = 0.0;

    #pragma omp parallel reduction(+:total_payoff)
    {
        unsigned int seed = time(NULL) ^ omp_get_thread_num();  // Unique seed for each thread

        #pragma omp for
        for (int i = 0; i < num_paths; i++) {
            double path_sum = 0.0;
            double S = S0;
            for (int j = 0; j < num_points; j++) {
                double z = gaussian_rng(&seed);
                S *= exp(drift + vol * z);
                path_sum += S;
            }
            total_payoff += fmax(path_sum / num_points - K, 0);
        }
    }

    return discount_factor * total_payoff / num_paths;
}

int main() {
    struct timeval start_time, end_time;
    clock_t start_cpu, end_cpu;

    gettimeofday(&start_time, NULL);  // Start the wall clock timing
    start_cpu = clock();              // Start the CPU timing

    int num_paths = 10000000;
    int num_points = 4096;
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;

    double price = asian_option_monte_carlo(num_paths, S0, K, T, r, sigma, num_points);
    
    end_cpu = clock();                // End the CPU timing
    gettimeofday(&end_time, NULL);    // End the wall clock timing

    double cpu_time_used = ((double) (end_cpu - start_cpu)) / CLOCKS_PER_SEC;
    double elapsed_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    printf("Estimated Asian option price using Monte Carlo simulation: $%.2f\n", price);
    printf("CPU execution time: %f seconds\n", cpu_time_used);
    printf("Wall clock execution time: %f seconds\n", elapsed_time);

    return 0;
}

