/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// FUNCTIONS
////////////

double radiative_transfer_unpolarized(double *lightpath, int steps,
                                      double *frequency,
                                      double IQUV[num_frequencies][4],
                                      double tau[num_frequencies]) {

    int path_counter;
    double pitch_ang, nu_p;

    double X_u[4], k_d[4], k_u[4], k_u_s[4], dl_current, dl_current_s;
    double jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV;

    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm
    double g_dd[4][4], g_uu[4][4];

    double Icurrent;

    double dtau_old = 0;

    struct GRMHD modvar;
    modvar.B = 0;
    modvar.n_e = 0.;
    modvar.theta_e = 0;

    LOOP_i {
        modvar.B_u[i] = 0;
        modvar.U_u[i] = 0;
        modvar.B_d[i] = 0;
        modvar.U_d[i] = 0;
    }
    modvar.igrid_c = -1;

    for (path_counter = steps - 1; path_counter > 0; path_counter--) {
        // Current position, wave vector, and dlambda
        LOOP_i {
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter - 1) * 9 + 8]);

        metric_dd(X_u, g_dd);
        metric_uu(X_u, g_uu);

        if (get_fluid_params(X_u, &modvar)) {
            lower_index(X_u, k_u, k_d);
            pitch_ang = pitch_angle(X_u, k_u, modvar.B_u, modvar.U_u);

            if (pitch_ang < 1e-9)
                continue;

            for (int f = 0; f < num_frequencies; f++) {
                // Obtain pitch angle: still no units (geometric)

                Icurrent = IQUV[f][0];

                dl_current_s =
                    dl_current *
                    (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
                    (PLANCK_CONSTANT * frequency[f]);

                // CGS UNITS USED FROM HERE ON OUT
                //////////////////////////////////

                // Scale the wave vector to correct energy
                LOOP_i k_u_s[i] =
                    k_u[i] * PLANCK_CONSTANT * frequency[f] /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

                // lower the index of the wavevector
                lower_index(X_u, k_u_s, k_d);

                // Compute the photon frequency in the plasma frame:
                nu_p = freq_in_plasma_frame(modvar.U_u, k_d);
                // Obtain emission coefficient in current plasma conditions

#if (EMISUSER)
                evaluate_coeffs_user(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI,
                                     &aQ, &aU, &aV, nu_p, modvar, pitch_ang);
#else
                evaluate_coeffs_single(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI,
                                       &aQ, &aU, &aV, nu_p, modvar, pitch_ang);
#endif
                double C = Rg * PLANCK_CONSTANT /
                           (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

                double dtau = (aI * dl_current_s * C + dtau_old);
                double K_inv = aI;
                double j_inv = jI;

                tau[f] += dtau;
#if (DEBUG)
                if ((j_inv != j_inv || isnan(Icurrent))) {
                    fprintf(stderr, "NaN emissivity! I = %+.15e\n", Icurrent);
                    fprintf(stderr, "NaN emissivity! j_nu = %+.15e\n", j_inv);
                    fprintf(stderr, "NaN emissivity! nu_plasmaframe = %+.15e\n",
                            nu_p);
                    fprintf(stderr, "NaN emissivity! ne %e te %e B %e\n",
                            modvar.n_e, modvar.theta_e, modvar.B);
                    fprintf(stderr, "NaN emissivity! Unorm %e\n",
                            four_velocity_norm(X_u, modvar.U_u) + 1);
                    fprintf(stderr, "NaN emissivity! knorm %e\n",
                            four_velocity_norm(X_u, k_u));
                }
#endif

                if (jI == jI) {
                    double Ii = Icurrent;
                    double S = j_inv / K_inv;
                    if (K_inv == 0)
                        Icurrent = Ii;
                    else if (dtau < 1.e-5)
                        Icurrent = Ii - (Ii - S) * (0.166666667 * dtau *
                                                    (6. - dtau * (3. - dtau)));
                    else {
                        double efac = exp(-dtau);
                        Icurrent = Ii * efac + S * (1. - efac);
                    }
                }
                dtau_old = 0;

                IQUV[f][0] = Icurrent;
            }
        }
    }

    // IQUV[0] = Icurrent * pow(frequency, 3.);

    return 1;
}

double BB_spectrum(double nu, double Temperature){

    double Intensity, cons1, cons2;

    
    cons1 = 2 * PLANCK_CONSTANT / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    cons2 = PLANCK_CONSTANT / BOLTZMANN_CONSTANT;
    Intensity = cons1 * nu * nu * nu / (exp(cons2 * nu / Temperature) - 1.);

    return Intensity;
}

double Temperature1(double rho, double pp, double gammarel, double X_u[4], double U_u[4]) {

    double U_d[4];
    double TMA[4][4], kr[4][4];
    double TMA_BL[4][4];
    double TMA_rt,Temp;

    lower_index(X_u, U_u, U_d);
    double sigmaa = 5.6704e-5; //Stefan Boltzmann Constant in cgs units  

    LOOP_ij TMA_BL[i][j] = 0.;
    LOOP_ij kr[i][j] = 0.;
    
    for (int i = 0; i < DIM; i++){
        kr[i][i] = 1.;
    }
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            TMA[i][j] = (rho + pp * gammarel/(gammarel - 1.)) * U_u[i] * U_d[j] + pp * kr[i][j];
        }
    }
    stress_BL(X_u, TMA, TMA_BL);
    TMA_rt = fabs(TMA_BL[1][0]) * RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    Temp = pow(TMA_rt/sigmaa, 1./4.);
    return Temp;
}

void star_BB_emission(double *lightpath, int steps,
                        double *frequency, double IQUV[num_frequencies][4],
                        double alpha, double beta, 
                        int block, int pixel, double phi_global) {

    int path_counter;

    double X_u[4], k_d[4], k_u[4], k_u_s[4], photon_CSSuu[8], photon_BLuu[8], uBL_u[4], uBL_d[4];
    double TMArt, Temp;
    double nu_p;

    double sigmaa = 5.6704e-5; //Stefan Boltzmann Constant in cgs units  

    struct GRMHD modvar;

    #if (PPM)

    LOOP_i {
        modvar.B_u[i] = 0;
        modvar.B_d[i] = 0;
        modvar.U_u[i] = 0;
        modvar.U_d[i] = 0;
    }
    modvar.igrid_c = -1;
    modvar.B = 0;
    modvar.prim_rho = 0;
    modvar.prim_pp = 0.;

    path_counter = steps-1;

    LOOP_i {
        X_u[i] = lightpath[path_counter * 9 + i];
        k_u[i] = lightpath[path_counter * 9 + 4 + i];
    }

    double radii = get_r(X_u);

    if (get_fluid_params_star(X_u, &modvar)) {
        
        lower_index(X_u, k_u, k_d);

        Temp = Temperature1(modvar.prim_rho, modvar.prim_pp, modvar.gammarel, X_u, modvar.U_u);
        Temp = Temp * modvar.tr1;
        for (int f = 0; f < num_frequencies; f++) {

            IQUV[f][0] = 0.; //Initializing everything to 0
            // Scale the wave vector to correct energy
            LOOP_i k_u_s[i] =
                k_u[i] * PLANCK_CONSTANT * frequency[f] /
                (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // lower the index of the wavevector
            lower_index(X_u, k_u_s, k_d);

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(modvar.U_u, k_d);

            IQUV[f][0] = BB_spectrum(nu_p, Temp)/(nu_p * nu_p * nu_p);
            write_starBB_output(X_u, IQUV, Temp, block, pixel, alpha, beta, frequency, phi_global);

        }
    }
    #endif
}

