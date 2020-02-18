// Part P0: Set the number of ghost cells, from NRPy+'s FD_CENTDERIVS_ORDER
#define NGHOSTS 3

// Step P1a: Import needed header files
#include "stdio.h"
#include "stdlib.h"
#include <stdint.h>
#include "math.h"
#include "time.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Step P1b: Import necessary gsl libraries for interpolating the initial data onto the grid

#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_interp.h"

// Step P2: Add needed #define's to set data type, the IDX4() macro, and the gridfunctions
// Step P2a: set REAL=double, so that all floating point numbers are stored to at least ~16 significant digits.
#define REAL double


// Step P3: Set free parameters
// Step P3a: Free parameters for the numerical grid

// Spherical coordinates parameter
const REAL RMAX    = 30.; /* Set to approximately the time you wish to evolve for, 
                            * so that at t=t_final data at the origin is not 
                            * affected by the boundary conditions */
// Time coordinate parameters
const REAL t_final =  50.;
const REAL CFL_FACTOR = 0.5; // Set the CFL Factor

// Step P3b: Free parameters for the spacetime evolution
const REAL eta = 2.; // Gamma-driving shift condition parameter.

// Step P4: Implement the algorithm for upwinding.
//          *NOTE*: This upwinding is backwards from
//          usual upwinding algorithms, because the
//          upwinding control vector in BSSN (the shift)
//          acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0

// Step P5: Set free parameters for Psi initial data
const REAL psi_posn_x = 0.0,psi_posn_y = 0.0,psi_posn_z = 0.0;

// Step P5b: Set free parameters for the scalar field
const REAL scalar_posn_x = 0.0;
const REAL scalar_posn_y = 0.0;
const REAL scalar_posn_z = 0.0;
const REAL br_on = 1.; // Turn on(1.)/off(0.) scalar field backreaction on the metric
const REAL pot1_on = 0.; // Turn on(1.)/off(0.) quadratic potential
const REAL pot2_on = 0.; // Turn on(1.)/off(0.) self-interactiong potential
// Make sure only one potential is on at a time
// Variables for the scalar field potential
const REAL scalarmass = 1.; // Scalar mass, \mu = c/\hbar m
const REAL fa = 0.05; // Decay constant, only relevant for the self-interacting potential

//Step P5c: Declare vars for initial data arrays
//          We use initial data profiles for the scalar 
//          and the conformal factor that is known to 
//          lead to stable scalar field evolution
REAL uu_in;
REAL vv_in;
REAL psi_in;
REAL alpha_in;
REAL r_scalar;
REAL r_psi;

// Step P6: Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//          data in a 1D array. In this case, consecutive values of "i" 
//          (all other indices held to a fixed value) are consecutive in memory, where 
//          consecutive values of "j" (fixing all other indices) are separated by 
//          Nxx_plus_2NGHOSTS[0] elements in memory. Similarly, consecutive values of
//          "k" are separated by Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1] in memory, etc.
#define IDX4(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * (k) ) )
// Assuming idx = IDX3(i,j,k). Much faster if idx can be reused over and over:
#define IDX4pt(g,idx)   ( (idx) + (Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2]) * (g) )

// Step P7: Set #define's for BSSN gridfunctions. C code generated above
#include "gridfunction_defines.h"

#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)

void xxCart(REAL *xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
#include "xxCart.h"
}

// Step P8: Include basic functions needed to impose curvilinear
//          parity and boundary conditions.
#include "curvilinear_parity_and_outer_boundary_conditions.h"

#include "enforce_detgammabar_constraint.h"

// Step P9: Find the CFL-constrained timestep
REAL find_timestep(const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3],REAL *xx[3], const REAL CFL_FACTOR) {
  const REAL dxx0 = dxx[0], dxx1 = dxx[1], dxx2 = dxx[2];
  REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.
  LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS[0]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[1]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[2]-NGHOSTS) {
    const REAL xx0 = xx[0][i0], xx1 = xx[1][i1], xx2 = xx[2][i2];
    REAL ds_dirn0, ds_dirn1, ds_dirn2;
#include "ds_dirn.h"
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
    // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2);
    dsmin = MIN(dsmin,MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2)));
  }
  return dsmin*CFL_FACTOR;
}

// Contains BSSN_ID() for arbitrary initial data array
#include "ID_array_psi.h" 

// Step P10: Declare the function for the exact solution. time==0 corresponds to the initial data.
void initial_data(const int Nxx_plus_2NGHOSTS[3],REAL *xx[3], REAL *in_gfs) {

// Step P11a: Declare initial data arrays
FILE *uu_file = fopen("BSSN_SF/InitialData/phiCC9.csv", "r");
FILE *vv_file = fopen("BSSN_SF/InitialData/PiCC9.csv", "r");
FILE *psi_file = fopen("BSSN_SF/InitialData/psiCC9.csv", "r");
FILE *alpha_file = fopen("BSSN_SF/InitialData/alphaCC9.csv", "r");

int temp;
int alen = 0;
while(fscanf(uu_file,"%lf\n",&temp)==1){
    alen++;
}
double r_arr[alen];
double uu_in_arr[alen];
double vv_in_arr[alen];
double psi_in_arr[alen];
double alpha_in_arr[alen];
    
rewind(uu_file);
for(int i=0;i<alen;i++){
    r_arr[i] = 0.01*i;
    fscanf(uu_file, "%lf\n", &uu_in_arr[i]);
    fscanf(vv_file, "%lf\n", &vv_in_arr[i]);
    fscanf(psi_file, "%lf\n", &psi_in_arr[i]);
    fscanf(alpha_file, "%lf\n", &alpha_in_arr[i]);
}

// Step P11b: Declare splines to interpolate onto the cartesian grid
    
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    
gsl_spline *spline_u = gsl_spline_alloc (gsl_interp_cspline, alen);
gsl_spline_init(spline_u, r_arr, uu_in_arr, alen);

gsl_spline *spline_v = gsl_spline_alloc (gsl_interp_cspline, alen);
gsl_spline_init(spline_v, r_arr, vv_in_arr, alen);

gsl_spline *spline_psi = gsl_spline_alloc (gsl_interp_cspline, alen);
gsl_spline_init(spline_psi, r_arr, psi_in_arr, alen);
    
gsl_spline *spline_alpha = gsl_spline_alloc (gsl_interp_cspline, alen);
gsl_spline_init(spline_alpha, r_arr, alpha_in_arr, alen);

#pragma omp parallel for
  LOOP_REGION(0,Nxx_plus_2NGHOSTS[0], 0,Nxx_plus_2NGHOSTS[1], 0,Nxx_plus_2NGHOSTS[2]) {
    const int idx = IDX3(i0,i1,i2);
    REAL xCart[3];
    xxCart(xx, i0,i1,i2, xCart);
    {
        r_psi = sqrt(pow(-psi_posn_x + xCart[0], 2) + pow(-psi_posn_y + xCart[1], 2) + pow(-psi_posn_z + xCart[2], 2));
        psi_in = gsl_spline_eval (spline_psi, r_psi, acc);
        alpha_in = gsl_spline_eval (spline_alpha, r_psi, acc);
    }
    BSSN_ID(xx[0][i0],xx[1][i1],xx[2][i2],xCart[0],xCart[1],xCart[2],
            &in_gfs[IDX4pt(HDD00GF,idx)],&in_gfs[IDX4pt(HDD01GF,idx)],&in_gfs[IDX4pt(HDD02GF,idx)],
            &in_gfs[IDX4pt(HDD11GF,idx)],&in_gfs[IDX4pt(HDD12GF,idx)],&in_gfs[IDX4pt(HDD22GF,idx)],
            &in_gfs[IDX4pt(TRKGF,idx)],
            &in_gfs[IDX4pt(ADD00GF,idx)],&in_gfs[IDX4pt(ADD01GF,idx)],&in_gfs[IDX4pt(ADD02GF,idx)],
            &in_gfs[IDX4pt(ADD11GF,idx)],&in_gfs[IDX4pt(ADD12GF,idx)],&in_gfs[IDX4pt(ADD22GF,idx)],
            &in_gfs[IDX4pt(LAMBDAU0GF,idx)],&in_gfs[IDX4pt(LAMBDAU1GF,idx)],&in_gfs[IDX4pt(LAMBDAU2GF,idx)],
            &in_gfs[IDX4pt(VETU0GF,idx)],&in_gfs[IDX4pt(VETU1GF,idx)],&in_gfs[IDX4pt(VETU2GF,idx)],
            &in_gfs[IDX4pt(BETU0GF,idx)],&in_gfs[IDX4pt(BETU1GF,idx)],&in_gfs[IDX4pt(BETU2GF,idx)],
            &in_gfs[IDX4pt(ALPHAGF,idx)],&in_gfs[IDX4pt(CFGF,idx)]);
    REAL xx0 = xCart[0];
    REAL xx1 = xCart[1];
    REAL xx2 = xCart[2];
    {
        r_scalar = sqrt(pow(-scalar_posn_x + xx0, 2) + pow(-scalar_posn_y + xx1, 2) + pow(-scalar_posn_z + xx2, 2));
        in_gfs[IDX4(UUGF, i0, i1, i2)] = gsl_spline_eval (spline_u, r_scalar, acc);
        in_gfs[IDX4(VVGF, i0, i1, i2)] = gsl_spline_eval (spline_v, r_scalar, acc);
    }
  }
}

// Step P12: Implement Hamiltonian constraint diagnostic
void Hamiltonian_constraint(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3], 
                            REAL *in_gfs, REAL *aux_gfs) {
#include "Hamiltonian.h"    
}

// Step P13: Declare the function to evaluate the BSSN RHSs
void rhs_eval(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3], const REAL *in_gfs,REAL *rhs_gfs) {
#include "BSSN_RHSs.h"
}

#include "ID_array_ADM.h"  

// main() function:
// Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Set up scalar wave initial data
// Step 2: Evolve scalar wave initial data forward in time using Method of Lines with RK4 algorithm,
//         applying quadratic extrapolation outer boundary conditions.
// Step 3: Output relative error between numerical and exact solution.
// Step 4: Free all allocated memory
int main(int argc, const char *argv[]) {

  // Step 0a: Read command-line input, error out if nonconformant
  if(argc != 4 || atoi(argv[1]) < NGHOSTS) {
      printf("Error: Expected one command-line argument: ./BSSNCurvilinear_Playground Nx0 Nx1 Nx2,\n");
      printf("where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
      printf("Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
      exit(1);
  }
  // Step 0b: Set up numerical grid structure, first in space...
  const int Nx0 = atoi(argv[1]);
  const int Nx1 = atoi(argv[2]);
  const int Nx2 = atoi(argv[3]);
  if(Nx0%2 != 0 || Nx1%2 != 0 || Nx2%2 != 0) {
    printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }
  const int Nxx[3] = { Nx0, Nx1, Nx2 };
  const int Nxx_plus_2NGHOSTS[3] = { Nxx[0]+2*NGHOSTS, Nxx[1]+2*NGHOSTS, Nxx[2]+2*NGHOSTS };
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2];
#include "xxminmax.h"

  // Step 0c: Allocate memory for gridfunctions
  REAL *evol_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *next_in_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *aux_gfs  = (REAL *)malloc(sizeof(REAL) * NUM_AUX_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k1_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k2_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k3_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k4_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);

  // Step 0d: Set up space and time coordinates
  // Step 0d.i: Set \Delta x^i on uniform grids.
  REAL dxx[3];
  for(int i=0;i<3;i++) dxx[i] = (xxmax[i] - xxmin[i]) / ((REAL)Nxx[i]);

  // Step 0d.ii: Set up uniform coordinate grids
  REAL *xx[3];
  for(int i=0;i<3;i++) {
    xx[i] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS[i]);
    for(int j=0;j<Nxx_plus_2NGHOSTS[i];j++) {
      xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxx[i]; // Cell-centered grid.
    }
  }

  // Step 0d.iii: Set timestep based on smallest proper distance between gridpoints and CFL factor 
  REAL dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);
  //printf("# Timestep set to = %e\n",(double)dt);
  int N_final = (int)(t_final / dt + 0.5); // The number of iterations in time.
                                           //Add 0.5 to account for C rounding down integers.

  // Step 0e: Find ghostzone mappings and parities:
  gz_map *bc_gz_map = (gz_map *)malloc(sizeof(gz_map)*Nxx_plus_2NGHOSTS_tot);
  parity_condition *bc_parity_conditions = (parity_condition *)malloc(sizeof(parity_condition)*Nxx_plus_2NGHOSTS_tot);
  set_up_bc_gz_map_and_parity_conditions(Nxx_plus_2NGHOSTS,xx,dxx,xxmin,xxmax,  bc_gz_map, bc_parity_conditions);

  // Step 1: Set up initial data to be exact solution at time=0:
  initial_data(Nxx_plus_2NGHOSTS, xx, evol_gfs);

  // Step 1b: Apply boundary conditions *FOR VALIDATION PURPOSES*
  apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, evol_gfs);
  enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);

  // Step 2: Evaluate Hamiltonian constraint violation
  //Hamiltonian_constraint(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);

  // Step 3: Start the timer, for keeping track of how fast the simulation is progressing.
  struct timespec start, end;
  //clock_gettime(CLOCK_REALTIME, &start);

  // Step 4: Integrate the initial data forward in time using the Method of Lines and RK4

  char filename2[100];
  sprintf(filename2,"BSSN_SF-evolution/quad_pot_uu_vv_cf.txt");
  FILE *evol = fopen(filename2, "w");
  for(int n=0;n<=N_final;n++) { // Main loop to progress forward in time.
    /***************************************************/
    /* Implement RK4 for Method of Lines timestepping: */
    /***************************************************/
    /* -= RK4: Step 1 of 4 =- */
    /* First evaluate k1 = RHSs expression             */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,evol_gfs, k1_gfs);
    /* Next k1 -> k1*dt, and then set the input for    */
    /*    the next RHS eval call to y_n+k1/2           */
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
      k1_gfs[i] *= dt;
      next_in_gfs[i] = evol_gfs[i] + k1_gfs[i]*0.5;
    }
    /* Finally, apply boundary conditions to           */
    /* next_in_gfs, so its data are set everywhere.    */
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, next_in_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, next_in_gfs);
    
    /* -= RK4: Step 2 of 4 =- */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,next_in_gfs, k2_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
      k2_gfs[i] *= dt;
      next_in_gfs[i] = evol_gfs[i] + k2_gfs[i]*0.5;
    }
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, next_in_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, next_in_gfs);

    /* -= RK4: Step 3 of 4 =- */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,next_in_gfs, k3_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
      k3_gfs[i] *= dt;
      next_in_gfs[i] = evol_gfs[i] + k3_gfs[i];
    }
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, next_in_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, next_in_gfs);

    /* -= RK4: Step 4 of 4 =- */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,next_in_gfs, k4_gfs);
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
      k4_gfs[i] *= dt;
      evol_gfs[i] += (1.0/6.0)*(k1_gfs[i] + 2.0*k2_gfs[i] + 2.0*k3_gfs[i] + k4_gfs[i]);
    }
    Hamiltonian_constraint(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, evol_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
    /* Output the solution of the scalar field and the conformal factor at diffrent time slices on a 2D grid */
    if(n%10 == 0) {
        char filename[100];
        sprintf(filename,"BSSN_SF-output2D/quad_pot_2d_t-%08d.txt",n);
        FILE *out2D = fopen(filename, "w");
        const int i0MIN=NGHOSTS; // In spherical, r=Delta r/2.
        const int i1mid=Nxx_plus_2NGHOSTS[1]/2;
        const int i2mid=Nxx_plus_2NGHOSTS[2]/2;
        LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS[0]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[1]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[2]-NGHOSTS) {
           REAL xx0 = xx[0][i0];
           REAL xx1 = xx[1][i1];
           REAL xx2 = xx[2][i2];
           REAL xCart[3];
#include "xxCart.h"
           int idx = IDX3(i0,i1,i2);
             ADMCart_ID(out2D, n*(double)dt , xx0, xx1, xx2, xCart[0], xCart[1], xCart[2],    
                        evol_gfs[IDX4pt(HDD00GF,idx)], evol_gfs[IDX4pt(HDD01GF,idx)], evol_gfs[IDX4pt(HDD02GF,idx)],
                        evol_gfs[IDX4pt(HDD11GF,idx)], evol_gfs[IDX4pt(HDD12GF,idx)], evol_gfs[IDX4pt(HDD22GF,idx)],
                        evol_gfs[IDX4pt(ADD00GF,idx)], evol_gfs[IDX4pt(ADD01GF,idx)], evol_gfs[IDX4pt(ADD02GF,idx)],
                        evol_gfs[IDX4pt(ADD11GF,idx)], evol_gfs[IDX4pt(ADD12GF,idx)], evol_gfs[IDX4pt(ADD22GF,idx)], 
                        evol_gfs[IDX4pt(TRKGF,idx)], 
                        evol_gfs[IDX4pt(LAMBDAU0GF,idx)],evol_gfs[IDX4pt(LAMBDAU1GF,idx)],evol_gfs[IDX4pt(LAMBDAU2GF,idx)],
                        evol_gfs[IDX4pt(VETU0GF,idx)],evol_gfs[IDX4pt(VETU1GF,idx)],evol_gfs[IDX4pt(VETU2GF,idx)], 
                        evol_gfs[IDX4pt(BETU0GF,idx)],evol_gfs[IDX4pt(BETU1GF,idx)],evol_gfs[IDX4pt(BETU2GF,idx)],  
                        evol_gfs[IDX4pt(ALPHAGF,idx)],evol_gfs[IDX4pt(CFGF,idx)]); 
                        
                        
         }
         fclose(out2D);
    }
       
    // Output time evolution at r=0
    int idx0 = IDX3(0,0,0);
    fprintf(evol,"%e %e %e %e\n", n*dt, evol_gfs[IDX4pt(UUGF,idx0)],evol_gfs[IDX4pt(VVGF,idx0)],evol_gfs[IDX4pt(CFGF,idx0)]);
        
    // Progress indicator printing to stdout
    // Measure average time per iteration
    //clock_gettime(CLOCK_REALTIME, &end);
    const long long unsigned int time_in_ns = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
    const REAL s_per_iteration_avg = ((REAL)time_in_ns / (REAL)n) / 1.0e9;

    const int iterations_remaining = N_final - n;
    const REAL time_remaining_in_mins = s_per_iteration_avg * (REAL)iterations_remaining / 60.0;

    const REAL num_RHS_pt_evals = (REAL)(Nxx[0]*Nxx[1]*Nxx[2]) * 4.0 * (REAL)n; // 4 RHS evals per gridpoint for RK4
    const REAL RHS_pt_evals_per_sec = num_RHS_pt_evals / ((REAL)time_in_ns / 1.0e9);

    // Progress indicator printing to stdout
    printf("%c[2K", 27); // Clear the line
    printf("It: %d t=%.2f | %.1f%%; ETA %.0f s | t/h %.2f | gp/s %.2e\r",  // \r is carriage return, move cursor to the beginning of the line
           n, n * (double)dt, (double)(100.0 * (REAL)n / (REAL)N_final),
           (double)time_remaining_in_mins*60, (double)(dt * 3600.0 / s_per_iteration_avg), (double)RHS_pt_evals_per_sec);
    fflush(stdout); // Flush the stdout buffer
  } // End main loop to progress forward in time.
  printf("\n"); // Clear the line.
  fclose(evol);
                               
  /* Step 4: Free all allocated memory */
  free(bc_parity_conditions);
  free(bc_gz_map);
  free(k4_gfs);
  free(k3_gfs);
  free(k2_gfs);
  free(k1_gfs);
  free(aux_gfs);
  free(next_in_gfs);
  free(evol_gfs);
  for(int i=0;i<3;i++) free(xx[i]);
  return 0;
}