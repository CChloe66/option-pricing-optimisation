#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

double gaussian_rng(unsigned int *seed) {
    double u = (double)rand_r(seed) / RAND_MAX;
    double v = (double)rand_r(seed) / RAND_MAX;
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

void monte_carlo_asian_option(int num_paths, double S0, double K, double T, double r, double sigma, int num_points) {
    double dt = T / num_points;
    double drift = (r - 0.5 * sigma * sigma) * dt;
    double vol = sigma * sqrt(dt);
    double discount_factor = exp(-r * T);
    double *results = (double*) malloc(num_paths * sizeof(double));

    #pragma omp parallel
    {
        unsigned int seed = time(NULL) ^ omp_get_thread_num();

        #pragma omp single
        {
            for (int i = 0; i < num_paths; i++) {
                #pragma omp task firstprivate(i, seed)
                {
                    // Stage 1: Initialize path
                    double S = S0;
                    double path_sum = 0.0;

                    // Stage 2: Compute path
                    for (int j = 0; j < num_points; j += 2) {
                        double z1 = gaussian_rng(&seed);
                        double z2 = gaussian_rng(&seed);
                        S *= exp(drift + vol * z1);
                        path_sum += S;
                        if (j + 1 < num_points) {
                            S *= exp(drift + vol * z2);
                            path_sum += S;
                        }
                    }

                    // Stage 3: Calculate payoff
                    results[i] = fmax(path_sum / num_points - K, 0) * discount_factor;
                }
            }
        }
    }

    double total_payoff = 0.0;
    for (int i = 0; i < num_paths; i++) {
        total_payoff += results[i];
    }
    double price = total_payoff / num_paths;
    free(results);

    printf("Estimated Asian option price using Monte Carlo simulation: $%.2f\n", price);
}

int main() {
    srand((unsigned)time(NULL));

    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL); // Start the wall clock timing

    clock_t start_cpu = clock(); // Start the CPU timing

    int num_paths = 10000000;
    int num_points = 4096;
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;

    monte_carlo_asian_option(num_paths, S0, K, T, r, sigma, num_points);

    clock_t end_cpu = clock(); // End the CPU timing
    gettimeofday(&end_time, NULL); // End the wall clock timing

    double cpu_time_used = ((double) (end_cpu - start_cpu)) / CLOCKS_PER_SEC;
    double wall_time_used = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    printf("CPU execution time: %f seconds\n", cpu_time_used);
    printf("Wall clock execution time: %f seconds\n", wall_time_used);

    return 0;
}

