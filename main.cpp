#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstring>
#include <fftw3.h> // FFTW for Fourier Transforms
#include <fstream> // fstream for generating CSV output file
#include <chrono>  // chrono for calculating elapsed time

using namespace std;
using namespace chrono;
using Complex = complex<double>;

//################### 00. Function Prototype Declarations ####################
    // 00.01 forcing function: F = compute_forcing()
    
void compute_AntiAlias( double Lx, const fftw_complex *u_hat);
    
    // 00.02 rhs function: du_hat_dt = compute_rhs(nu,Lx,k,u_hat)
void compute_rhs( double nu, double Lx, const vector <double> &kx, const fftw_complex *u_hat, vector<double> &du_hat_dt);

int main(){
    // 01.01.Parameters
    const double Lx{2.0*M_PI};         // Length of periodic domain
    const int N = pow(2,5);            // Number of spatial grid points, degree of freedom
    const double nu{0.001};            // Kinematic viscosity
    const double tf = 1;               // Final time
    const double dt = 1e-5;            // Time step
    const size_t Nt = static_cast<size_t>(round(tf / dt)); // Number of time steps

    // 01.02. Spatial grid
    vector <double> x (N,0);
    for(int i{0} ; i<N ; ++i)  x[i] = i*(Lx/N);     // i=0:N-1; (N-1 points + zero), where N's point is @ x[N-1]
                                                    // Lx divided to N-1 portion with N points
    // 01.03. Initial velocity field: u(x,0)= -U0 * sin(kx)
    const int n{1};                       // Mode of initial wavenumber
    const double k = (n* 2 * M_PI / Lx);  // Wave number of initial velocity field
    const double U0{1.0};                 // Amplitude of the initial velocity field
    vector <double> u (N,0);
    for(int i{0} ; i < N ; ++i){
        u[i] = -U0 * sin(k * x[i]);
    }
    // 3.1. Enforcing initial periodic condition
    u[N-1] = u[0];
    
    // 01.04. generating CSV file of initial velocity field
    ofstream outFile1("C:/Users/kiara/CppWorkspaces/Workspace002/Burgers_Forced/velocity_field_Initial.csv");  //ofstream is a class from <fstream> library
    // Write header
    outFile1 << "x,velocity\n";
    if (!outFile1) {
        cerr << "Error: Could not open initial data file for writing!" << endl;
        return -1;
    }
    // Write data
    for (int i{0}; i < N; ++i) outFile1 << x[i] << "," << u[i] << "\n";
    outFile1.close();
    
    // 01.05. Wave numbers
    vector <double> kx (N);
    for (int i = 0; i < N; ++i) kx[i] = (i < N / 2) ? i * 2.0 * M_PI / Lx : (i - N) * 2.0 * M_PI / Lx;

    // 01.06. Forcing term
    
//################### 02. Initialize FFTW PLAN ####################
    size_t fft_size = sizeof(fftw_complex) * N;
    // 02.02 Allocate memory for Fourier coefficients
    fftw_complex *u_hat{nullptr};
    u_hat = (fftw_complex*) fftw_malloc(fft_size);    // u_hat is a pointer of the type fftw_complex
    memset(u_hat, 0, fft_size);  //  Zero Initialize
        
//  02.04. Allocate intermediate arrays to improve performance 
    vector<double> R1(N,0);
    vector<double> R2(N,0);
    vector<double> R3(N,0);
    
    vector<double> u1(N,0);
            
    // 02.03 Create FFTW plans
    fftw_plan forward = fftw_plan_dft_r2c_1d(N, u.data(), u_hat, FFTW_ESTIMATE);
    // fftw_plan backward = fftw_plan_dft_c2r_1d(N, u_hat, u.data(), FFTW_ESTIMATE);
    
    // 02.04 Fourier transform
    
//################### 03. Processing Stage ####################  
auto start_time = chrono::high_resolution_clock::now();
  
//// 09. Main time-stepping loop
for (size_t t_cuntr{0}; t_cuntr < Nt+1; ++t_cuntr){
    
    fftw_execute(forward);
// 03.02.01. RK3 stage 1
    vector<double> &u1 = u;
    compute_rhs(nu, Lx, kx, u_hat, R1);
    
    // 03.02.02. RK3 stage 2
    for (int i = 0; i < N; ++i) {  
        u[i] = u1[i] + 0.5 * dt *  R1[i];  // Element-wise operation
    }
    fftw_execute(forward);
    compute_rhs(nu, Lx, kx, u_hat, R2);
    
    // 03.02.03. RK3 stage 3
    for (int i = 0; i < N; ++i) {  
        u[i] = u1[i] + 0.75 * dt * R2[i];  // Element-wise operation
    }
    fftw_execute(forward);
    compute_rhs(nu, Lx, kx, u_hat, R3);
    
// 03.02.04. Final stage
    for (int i = 0; i < N; ++i) {  
        u[i] = u1[i] + (dt / 9.0) * ( 2.0*R1[i] + 3.0*R2[i] + 4.0*R3[i] );
    }
    // 03.04. Enforcing final step periodic condition
    u[N-1] = u[0];
}

//########################## End of Time-stepping loop ##############################
    
// 03.08. End measuring time
auto end_time = chrono::high_resolution_clock::now();
auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
cout << "Elapsed time for the time-stepping loop: " << duration.count() << " ms" << endl; 
   
// 03.06. Clean up
    fftw_destroy_plan(forward);
    //fftw_destroy_plan(backward);
    fftw_free(u_hat);
    
    return 0;
}

//################### 04. Function Definitions #################### 
    
    // forcing function: F = compute_forcing()
    
    
//// rhs function: du_hat_dt = compute_rhs(nu,Lx,k,u_hat)
    void compute_rhs( double nu, double Lx, const vector <double> &kx, const fftw_complex *u_hat, vector<double> &du_dt){
    int N = kx.size();
 
    /// Diffusion term: (d²u/dx²)
        fftw_complex *ddu_hat_ddx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
        vector<double> ddu_ddx(N);
        // Create IFFT plan for diffusion term (d²u/dx²)
        fftw_plan ifft_plan_diffusion = fftw_plan_dft_c2r_1d(N, ddu_hat_ddx, ddu_ddx.data(), FFTW_ESTIMATE);
        // Compute second derivative in Fourier space (k^2 * u_hat)
        for (int i{0}; i < N; ++i) {
            ddu_hat_ddx[i][0] = pow(kx[i], 2) * u_hat[i][0];  // Real part
            ddu_hat_ddx[i][1] = pow(kx[i], 2) * u_hat[i][1];
        }
        // Execute inverse FFT for diffusion term (d²u/dx²)
        fftw_execute(ifft_plan_diffusion);
        // Normalize the IFFT results
        for (int i = 0; i < N; ++i) {
            ddu_ddx[i] /= N;
        }

    /// Nonlinear-advection term: (u.(du/dx))
        fftw_complex *nonlinear_advection_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
        vector<double> nonlinear_advection(N);
        // Create IFFT plan nonlinear-advection term (u.(du/dx))
        fftw_plan ifft_plan_nonlinear_advection = fftw_plan_dft_c2r_1d(N, nonlinear_advection_hat, nonlinear_advection.data(), FFTW_ESTIMATE);
        // Compute nonlinear-advection term using AntiAlias function
        for (int i = 0; i < N; ++i) {
            nonlinear_advection_hat[i] = compute_AntiAlias(Lx, u_hat );
        }
        // Execute inverse FFT for nonlinear-advection term (u.(du/dx))
        fftw_execute(ifft_plan_nonlinear_advection);
        // Normalize the IFFT results
        for (int i = 0; i < N; ++i) {
            nonlinear_advection[i] /= N;
        }
    /// Total RHS term
        for (int i = 0; i < N; ++i) {
            du_dt[i] = Nonlinear_advection[i] - nu*ddu_ddx[i];
        }
        
    fftw_destroy_plan(ifft_plan_diffusion);
    fftw_destroy_plan(ifft_plan_nonlinear_advection);
    fftw_free(ddu_hat_ddx);
    fftw_free(nonlinear_advection_hat);
}

//// Anti Aliasing: compute_AntiAlias(u_hat,Lx)     
    void AntiAlias_ZeroPadding(double Lx, vector<fftw_complex>& u_hat) {
    int N = u_hat.size();
    int M = (3 * N) / 2; // 3/2 Rule

    // Step 1: Zero-padding in Fourier space
    vector<fftw_complex> u_hat_pad(M, {0.0, 0.0});
    for (int i = 0; i < N / 2; ++i) {
        u_hat_pad[i][0] = u_hat[i][0];   
        u_hat_pad[i][1] = u_hat[i][1];   
        u_hat_pad[M - (N / 2) + i][0] = u_hat[N / 2 + i][0];  
        u_hat_pad[M - (N / 2) + i][1] = u_hat[N / 2 + i][1];
    }

    // Step 2: Perform IFFT to go to physical space
    vector<fftw_complex> w_pad(M);
    fftw_plan ifft_plan = fftw_plan_dft_1d(M, u_hat_pad.data(), w_pad.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft_plan);
    fftw_destroy_plan(ifft_plan);

    // Normalize IFFT
    for (int i = 0; i < M; ++i) {
        w_pad[i][0] /= M;
        w_pad[i][1] /= M;
    }

    // Step 3: Compute wave numbers in extended space
    vector<double> k_extended(M);
    for (int i = 0; i < M; ++i) 
        k_extended[i] = (i < M / 2) ? i * 2.0 * M_PI / Lx : (i - M) * 2.0 * M_PI / Lx;

    // Step 4: Compute derivative (du/dx) in Fourier space
    vector<fftw_complex> du_dx_extended(M);
    for (int i = 0; i < M; ++i) {
        du_dx_extended[i][0] = -k_extended[i] * w_pad[i][1];  
        du_dx_extended[i][1] = k_extended[i] * w_pad[i][0];   
    }

    // Step 5: Compute nonlinear term (-u * du/dx) in physical space
    vector<fftw_complex> nonlinear_term_extended(M);
    for (int i = 0; i < M; ++i) {
        nonlinear_term_extended[i][0] = -w_pad[i][0] * du_dx_extended[i][0];  
        nonlinear_term_extended[i][1] = -w_pad[i][1] * du_dx_extended[i][1];  
    }

    // Step 6: Transform back to Fourier space
    fftw_plan fft_plan = fftw_plan_dft_1d(M, nonlinear_term_extended.data(), nonlinear_term_extended.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);

    // Step 7: Modify `u_hat` directly (Project back to original grid)
    for (int i = 0; i < N / 2; ++i) {
        u_hat[i][0] = (3.0 / 2.0) * nonlinear_term_extended[i][0];
        u_hat[i][1] = (3.0 / 2.0) * nonlinear_term_extended[i][1];
        u_hat[N / 2 + i][0] = (3.0 / 2.0) * nonlinear_term_extended[M - (N / 2) + i][0];
        u_hat[N / 2 + i][1] = (3.0 / 2.0) * nonlinear_term_extended[M - (N / 2) + i][1];
    }
}
