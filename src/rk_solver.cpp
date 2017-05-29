//
// Created by arrouan on 18/01/17.
//

#include "rk_solver.h"

void rk_solver::solve() {
  for (current_step=0;current_step<nb_steps;current_step++) {

    if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
      h=x2-x;

    bool ok_or_failed_step;
    switch(solver_method_type) {
      case 0:
        ok_or_failed_step = rk4();
        break;
      default:
        printf("Unknown method\n");
        exit(-1);
    }

    if (ok_or_failed_step)
      ok_step++;
    else
      failed_step++;

    if ((x-x2)*(x2-x1) >= 0.0) {
      for (int i=0;i<nvar;i++) ystart[i]=y[i];
      if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS)
        out.save(x,y);
      return;
    }

    if (abs(next_step_size) <= min_step_size) throw("Step size too small in Odeint");
    h=s.hnext;
  }
}

void rk_solver::rk4(int current_time, int step_size) {
    Int n=y.size();
    VecDoub dym(n),dyt(n),yt(n);
    double step_1=step_size*0.5;
    double step_2=step_size/6.0;
    double step_3=current_time+step_1;

    // First step
    for (int i=0;i<n;i++)
      yt[i]=y[i]+hh*dydx[i];

    // Second step
    derivs(xh,yt,dyt);
    for (int i=0;i<n;i++)
      yt[i]=y[i]+hh*dyt[i];

    // Third step
    derivs(xh,yt,dym);
    for (int i=0;i<n;i++) {
      yt[i]=y[i]+h*dym[i];
      dym[i] += dyt[i];
    }

    // Fourth step
    derivs(x+h,yt,dyt);
    for (int i=0;i<n;i++)
      yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

void rk_solver::raevol_func(int indiv_id, int prot_id, int degradstep) {
  for (int j = 0;
       j < nb_rna_produce_protein[indiv_id][prot_id]; j++) {
    double enhancer_activity = 0;
    double operator_activity = 0;

    int rna_id = rna_produce_protein_array[indiv_id][prot_id][j];

    for (int i = 0; i <
                    nb_rna_influence_enhancing_coef[indiv_id][rna_id]; i++) {

      enhancer_activity +=
          rna_influence_enhancing_coef_array[indiv_id][rna_id][i]
          * protein_concentration_array[indiv_id][i];
      operator_activity +=
          rna_influence_operating_coef_array[indiv_id][rna_id][i]
          * protein_concentration_array[indiv_id][i];
    }

    double enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                     powf(enhancer_activity, hill_shape_n);
    double operator_activity_pow_n = operator_activity == 0 ? 0 :
                                     powf(operator_activity, hill_shape_n);
    delta += rna_basal_concentration_array[indiv_id][rna_id]
             * (hill_shape
                / (operator_activity_pow_n + hill_shape))
             * (1 +
                ((1 / rna_basal_concentration_array[indiv_id][rna_id]
                 ) -
                 1)
                * (enhancer_activity_pow_n /
                   (enhancer_activity_pow_n + hill_shape)));
  }

  delta -=
      degradrate *
      protein_concentration_array[indiv_id][prot_id];
  delta *= 1 / (double) degradstep;
}
