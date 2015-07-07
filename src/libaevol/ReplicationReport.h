// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
//*****************************************************************************


#ifndef AEVOL_REPLICATION_REPORT_H__
#define AEVOL_REPLICATION_REPORT_H__


// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>
#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "DnaReplicationReport.h"
#include "ae_enums.h"
#include "Observer.h"
#include "ObservableEvent.h"
#include "Individual.h"

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class Individual;





class ReplicationReport : public Observer {
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ReplicationReport() = default;
    ReplicationReport(Individual* indiv,
                      const Individual* parent,
                      Individual* donor = NULL);

    // Creates a completely independent copy of the original report
    ReplicationReport(const ReplicationReport& other);

    ReplicationReport(gzFile tree_file, Individual * indiv);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ReplicationReport(void) = default;

    // =================================================================
    //                              Accessors
    // =================================================================
    inline Individual * get_indiv(void) const;
    int32_t id() { return _id; };
    int32_t rank() { return _rank; };
    inline int32_t  get_genome_size(void) const;
    inline int32_t  get_parent_id(void) const;
    inline double   get_parent_metabolic_error(void) const;
    inline double   get_parent_secretion_error(void) const;
    inline int32_t  get_parent_genome_size(void) const;

    inline double   get_mean_align_score(void) const;
    inline int32_t	get_donor_id(void) const;
    inline double   get_donor_metabolic_error(void) const;
    inline int32_t  get_donor_genome_size(void) const;
    int32_t get_nb(MutationType t) const {
      return _dna_replic_report.get_nb(t);
    }

    // TODO <david.parsons@inria.fr> re-constify
    // => const DnaReplicationReport& dna_replic_report() const
    DnaReplicationReport& dna_replic_report() {
      return _dna_replic_report;
    }
     
    void            set_indiv(Individual * indiv);
    inline void     set_parent_id(int32_t parent_id);
    inline void     set_parent_metabolic_error(double parent_metabolic_error);
    inline void     set_parent_secretion_error(double parent_secretion_error);
    inline void     set_parent_genome_size(int32_t parent_genome_size);
    inline void     set_donor_id(int32_t donor_id);
    inline void     set_donor_metabolic_error(double donor_metabolic_error);
    inline void     set_donor_secretion_error(double donor_secretion_error);
    inline void     set_donor_genome_size(int32_t donor_genome_size);

    // =================================================================
    //                            Public Methods
    // =================================================================
    void init(Individual* offspring, Individual* parent);
    void signal_end_of_replication(Individual* indiv);
    void signal_end_of_generation();
    void write_to_tree_file(gzFile tree_file) const;


  void update(Observable& o, ObservableEvent e, void* arg) override;





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    Individual* _indiv = nullptr;
    int32_t _id = -1;
    int32_t _parent_id = -1;

    int32_t _rank = -1;

    int32_t _genome_size = -1;
    double _metabolic_error = -1;
    int16_t _nb_genes_activ = -1;
    int16_t _nb_genes_inhib = -1;
    int16_t _nb_non_fun_genes = -1;
    int16_t _nb_coding_RNAs = -1;
    int16_t _nb_non_coding_RNAs = -1;

    // List of each genetic unit's replication report
    DnaReplicationReport _dna_replic_report;
    
    double _parent_metabolic_error = -1;
    double _parent_secretion_error = -1;
    int32_t _parent_genome_size = -1;
    
    int32_t	_donor_id = -1;
    double _donor_metabolic_error = -1;
    double _donor_secretion_error = -1;
    int32_t _donor_genome_size = -1;
    
    // CK: I think that the attributes below are obsolete 
    // (HT events are now stored in the ae_dna_replic_reports)
    
    double _mean_align_score;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline Individual *ReplicationReport::get_indiv(void) const
{
  return _indiv;
}

inline int32_t ReplicationReport::get_genome_size(void) const
{
  return _genome_size;
}

int32_t ReplicationReport::get_parent_id(void) const
{
  return _parent_id;
}

double ReplicationReport::get_parent_metabolic_error(void) const
{
  return _parent_metabolic_error;
}

double ReplicationReport::get_parent_secretion_error(void) const
{
  return _parent_secretion_error;
}

int32_t ReplicationReport::get_parent_genome_size(void) const
{
  return _parent_genome_size;
}

inline int32_t	ReplicationReport::get_donor_id(void) const
{
  return _donor_id;
}

inline double   ReplicationReport::get_donor_metabolic_error(void) const
{
  return _donor_metabolic_error;
}

inline int32_t  ReplicationReport::get_donor_genome_size(void) const
{
  return _donor_genome_size;
}



inline double ReplicationReport::get_mean_align_score(void) const
{
  return _mean_align_score;
}

inline void ReplicationReport::set_indiv(Individual * indiv)
{
  _indiv = indiv;
}

void ReplicationReport::set_parent_id(int32_t parent_id)
{
  _parent_id = parent_id;
}

void ReplicationReport::set_parent_metabolic_error(double parent_metabolic_error)
{
  _parent_metabolic_error = parent_metabolic_error;
}

void ReplicationReport::set_parent_secretion_error(double parent_secretion_error)
{
  _parent_secretion_error = parent_secretion_error;
}

void ReplicationReport::set_parent_genome_size(int32_t parent_genome_size)
{
  _parent_genome_size = parent_genome_size;
}

inline void ReplicationReport::set_donor_id(int32_t donor_id)
{
  _donor_id = donor_id;
}

inline void  ReplicationReport::set_donor_metabolic_error(double donor_metabolic_error)
{
  _donor_metabolic_error = donor_metabolic_error;
}

inline void ReplicationReport::set_donor_secretion_error(double donor_secretion_error)
{
  _donor_secretion_error = donor_secretion_error;
}

inline void ReplicationReport::set_donor_genome_size(int32_t donor_genome_size)
{
  _donor_genome_size = donor_genome_size;
}

} // namespace aevol

#endif // AEVOL_REPLICATION_REPORT_H__
