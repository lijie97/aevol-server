//
// Created by arrouan on 18/01/17.
//

#ifndef RAEVOL_CUDA_RK_SOLVER_H
#define RAEVOL_CUDA_RK_SOLVER_H

class rk_solver {
 public:
    void  solve();
    void  rk4();

 private:
    // Solver type 0=rk45
    int solver_method_type;

    // How many steps
    int nb_steps;
    // Current timestep
    int current_step;
    // Number of good steps (i.e. no need to recompute)
    int ok_step;
    // Number of bad steps (i.e. need to recompute)
    int failed_step;

    // Minimum step size
    double min_step_size;
    // Next step size
    double next_step_size;

    // Absolute tolerance
    double abs_tolerance;
    // Relative tolerance
    double rel_tolerance;
};


#endif //RAEVOL_CUDA_RK_SOLVER_H
