//
// Created by dparsons on 31/05/16.
//

#ifndef AEVOL_INDIVANALYSIS_H__
#define AEVOL_INDIVANALYSIS_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include "Individual.h"

namespace aevol {

/**
 *
 */
class IndivAnalysis : public Individual {
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  IndivAnalysis() = default; //< Default ctor
  IndivAnalysis(const IndivAnalysis&) = delete; //< Copy ctor
  IndivAnalysis(IndivAnalysis&&) = delete; //< Move ctor

  IndivAnalysis(const Individual&);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~IndivAnalysis() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  IndivAnalysis& operator=(const IndivAnalysis& other) = delete;
  /// Move assignment
  IndivAnalysis& operator=(IndivAnalysis&& other) = delete;

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  double compute_theoritical_f_nu();
  void compute_experimental_f_nu(int32_t nb_indiv,
                                 std::shared_ptr<JumpingMT> prng,
                                 FILE* output_summary = nullptr,
                                 FILE* output_detailed = nullptr,
                                 bool verbose = false);

  // ==========================================================================
  //                                Accessors
  // ==========================================================================

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
};

} // namespace aevol
#endif //AEVOL_INDIVANALYSIS_H__
