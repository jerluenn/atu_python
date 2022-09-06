/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_tetherunit_integrator.h"

#define NX     TETHERUNIT_INTEGRATOR_NX
#define NZ     TETHERUNIT_INTEGRATOR_NZ
#define NU     TETHERUNIT_INTEGRATOR_NU
#define NP     TETHERUNIT_INTEGRATOR_NP
#define NBX    TETHERUNIT_INTEGRATOR_NBX
#define NBX0   TETHERUNIT_INTEGRATOR_NBX0
#define NBU    TETHERUNIT_INTEGRATOR_NBU
#define NSBX   TETHERUNIT_INTEGRATOR_NSBX
#define NSBU   TETHERUNIT_INTEGRATOR_NSBU
#define NSH    TETHERUNIT_INTEGRATOR_NSH
#define NSG    TETHERUNIT_INTEGRATOR_NSG
#define NSPHI  TETHERUNIT_INTEGRATOR_NSPHI
#define NSHN   TETHERUNIT_INTEGRATOR_NSHN
#define NSGN   TETHERUNIT_INTEGRATOR_NSGN
#define NSPHIN TETHERUNIT_INTEGRATOR_NSPHIN
#define NSBXN  TETHERUNIT_INTEGRATOR_NSBXN
#define NS     TETHERUNIT_INTEGRATOR_NS
#define NSN    TETHERUNIT_INTEGRATOR_NSN
#define NG     TETHERUNIT_INTEGRATOR_NG
#define NBXN   TETHERUNIT_INTEGRATOR_NBXN
#define NGN    TETHERUNIT_INTEGRATOR_NGN
#define NY0    TETHERUNIT_INTEGRATOR_NY0
#define NY     TETHERUNIT_INTEGRATOR_NY
#define NYN    TETHERUNIT_INTEGRATOR_NYN
#define NH     TETHERUNIT_INTEGRATOR_NH
#define NPHI   TETHERUNIT_INTEGRATOR_NPHI
#define NHN    TETHERUNIT_INTEGRATOR_NHN
#define NPHIN  TETHERUNIT_INTEGRATOR_NPHIN
#define NR     TETHERUNIT_INTEGRATOR_NR


int main()
{

    tetherunit_integrator_solver_capsule *acados_ocp_capsule = tetherunit_integrator_acados_create_capsule();
    // there is an opportunity to change the number of shooting intervals in C without new code generation
    int N = TETHERUNIT_INTEGRATOR_N;
    // allocate the array and fill it accordingly
    double* new_time_steps = NULL;
    int status = tetherunit_integrator_acados_create_with_discretization(acados_ocp_capsule, N, new_time_steps);

    if (status)
    {
        printf("tetherunit_integrator_acados_create() returned status %d. Exiting.\n", status);
        exit(1);
    }

    ocp_nlp_config *nlp_config = tetherunit_integrator_acados_get_nlp_config(acados_ocp_capsule);
    ocp_nlp_dims *nlp_dims = tetherunit_integrator_acados_get_nlp_dims(acados_ocp_capsule);
    ocp_nlp_in *nlp_in = tetherunit_integrator_acados_get_nlp_in(acados_ocp_capsule);
    ocp_nlp_out *nlp_out = tetherunit_integrator_acados_get_nlp_out(acados_ocp_capsule);
    ocp_nlp_solver *nlp_solver = tetherunit_integrator_acados_get_nlp_solver(acados_ocp_capsule);
    void *nlp_opts = tetherunit_integrator_acados_get_nlp_opts(acados_ocp_capsule);

    // initial condition
    int idxbx0[NBX0];
    idxbx0[0] = 0;
    idxbx0[1] = 1;
    idxbx0[2] = 2;
    idxbx0[3] = 3;
    idxbx0[4] = 4;
    idxbx0[5] = 5;
    idxbx0[6] = 6;
    idxbx0[7] = 7;
    idxbx0[8] = 8;
    idxbx0[9] = 9;
    idxbx0[10] = 10;
    idxbx0[11] = 11;
    idxbx0[12] = 12;
    idxbx0[13] = 13;
    idxbx0[14] = 14;
    idxbx0[15] = 15;
    idxbx0[16] = 16;

    double lbx0[NBX0];
    double ubx0[NBX0];
    lbx0[0] = 0;
    ubx0[0] = 0;
    lbx0[1] = 0;
    ubx0[1] = 0;
    lbx0[2] = 0;
    ubx0[2] = 0;
    lbx0[3] = 1;
    ubx0[3] = 1;
    lbx0[4] = 0;
    ubx0[4] = 0;
    lbx0[5] = 0;
    ubx0[5] = 0;
    lbx0[6] = 0;
    ubx0[6] = 0;
    lbx0[7] = -5;
    ubx0[7] = 5;
    lbx0[8] = -5;
    ubx0[8] = 5;
    lbx0[9] = -5;
    ubx0[9] = 5;
    lbx0[10] = -5;
    ubx0[10] = 5;
    lbx0[11] = -5;
    ubx0[11] = 5;
    lbx0[12] = -5;
    ubx0[12] = 5;
    lbx0[13] = 0;
    ubx0[13] = 0;
    lbx0[14] = 0;
    ubx0[14] = 0;
    lbx0[15] = 0.05;
    ubx0[15] = 0.05;
    lbx0[16] = 0;
    ubx0[16] = 0;

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);

    // initialization for state values
    double x_init[NX];
    x_init[0] = 0.0;
    x_init[1] = 0.0;
    x_init[2] = 0.0;
    x_init[3] = 0.0;
    x_init[4] = 0.0;
    x_init[5] = 0.0;
    x_init[6] = 0.0;
    x_init[7] = 0.0;
    x_init[8] = 0.0;
    x_init[9] = 0.0;
    x_init[10] = 0.0;
    x_init[11] = 0.0;
    x_init[12] = 0.0;
    x_init[13] = 0.0;
    x_init[14] = 0.0;
    x_init[15] = 0.0;
    x_init[16] = 0.0;

    // initial value for control input
    double u0[NU];
    u0[0] = 0.0;

    // prepare evaluation
    int NTIMINGS = 1;
    double min_time = 1e12;
    double kkt_norm_inf;
    double elapsed_time;
    int sqp_iter;

    double xtraj[NX * (N+1)];
    double utraj[NU * N];


    // solve ocp in loop
    int rti_phase = 0;

    for (int ii = 0; ii < NTIMINGS; ii++)
    {
        // initialize solution
        for (int i = 0; i <= nlp_dims->N; i++)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x_init);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
        }
        ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_phase", &rti_phase);
        status = tetherunit_integrator_acados_solve(acados_ocp_capsule);
        ocp_nlp_get(nlp_config, nlp_solver, "time_tot", &elapsed_time);
        min_time = MIN(elapsed_time, min_time);
    }

    /* print solution and statistics */
    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", &xtraj[ii*NX]);
    for (int ii = 0; ii < nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", &utraj[ii*NU]);

    printf("\n--- xtraj ---\n");
    d_print_exp_tran_mat( NX, N+1, xtraj, NX);
    printf("\n--- utraj ---\n");
    d_print_exp_tran_mat( NU, N, utraj, NU );
    // ocp_nlp_out_print(nlp_solver->dims, nlp_out);

    printf("\nsolved ocp %d times, solution printed above\n\n", NTIMINGS);

    if (status == ACADOS_SUCCESS)
    {
        printf("tetherunit_integrator_acados_solve(): SUCCESS!\n");
    }
    else
    {
        printf("tetherunit_integrator_acados_solve() failed with status %d.\n", status);
    }

    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "kkt_norm_inf", &kkt_norm_inf);
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);

    tetherunit_integrator_acados_print_stats(acados_ocp_capsule);

    printf("\nSolver info:\n");
    printf(" SQP iterations %2d\n minimum time for %d solve %f [ms]\n KKT %e\n",
           sqp_iter, NTIMINGS, min_time*1000, kkt_norm_inf);

    // free solver
    status = tetherunit_integrator_acados_free(acados_ocp_capsule);
    if (status) {
        printf("tetherunit_integrator_acados_free() returned status %d. \n", status);
    }
    // free solver capsule
    status = tetherunit_integrator_acados_free_capsule(acados_ocp_capsule);
    if (status) {
        printf("tetherunit_integrator_acados_free_capsule() returned status %d. \n", status);
    }

    return status;
}
