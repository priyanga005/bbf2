#include <stdio.h>
#include <math.h>

// Butterworth IIR filter coefficients for a 2nd order low-pass filter
#define PI 3.14159265358979323846
#define CUTOFF_FREQ 1000.0  // Cutoff frequency in Hz
#define SAMPLING_FREQ 5000.0  // Sampling frequency in Hz

// Filter coefficients for 2nd order Butterworth low-pass
double b[3], a[3];

// Function to calculate filter coefficients for a Butterworth filter
void calculate_butterworth_coefficients() {
    double Wc = 2 * PI * CUTOFF_FREQ;  // Cutoff angular frequency
    double Ts = 1 / SAMPLING_FREQ;     // Sampling period
    double K = tan(Wc * Ts / 2);       // Bilinear transform

    // Coefficients for a second order low-pass Butterworth filter
    double norm = 1 / (1 + sqrt(2) * K + K * K);
    
    b[0] = K * K * norm;
    b[1] = 2 * b[0];
    b[2] = b[0];
    
    a[0] = 1;
    a[1] = 2 * (K * K - 1) * norm;
    a[2] = (1 - sqrt(2) * K + K * K) * norm;
}

// Apply the filter to a signal (difference equation)
double apply_filter(double x, double *y, double *x_hist) {
    // Shift history
    x_hist[0] = x_hist[1];
    x_hist[1] = x_hist[2];
    x_hist[2] = x;
    
    // Apply the filter difference equation
    y[0] = b[0] * x_hist[0] + b[
