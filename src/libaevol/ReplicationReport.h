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


 #ifndef __AE_REPLICATION_REPORT_H__
#define  __AE_REPLICATION_REPORT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "DnaReplicationReport.h"
#include "ae_enums.h"

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class Individual;





class ReplicationReport
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ReplicationReport(Individual * indiv,
                          const Individual * parent,
                          Individual * donor = NULL);

    // Creates a completely independent copy of the original report
    ReplicationReport(const ReplicationReport & model);

    ReplicationReport(gzFile tree_file, Individual * indiv);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ReplicationReport(void);

    // =================================================================
    //                              Accessors
    // =================================================================
    inline Individual * get_indiv(void) const;
    int32_t         get_id(void) const;
    int32_t         get_rank(void) const;
    inline int32_t  get_genome_size(void) const;

    inline int32_t  get_parent_id(void) const;
    inline double   get_parent_metabolic_error(void) const;
    inline double   get_parent_secretion_error(void) const;
    inline int32_t  get_parent_genome_size(void) const;

    inline double   get_mean_align_score(void) const;
    inline int32_t	get_donor_id(void) const;
    inline double   get_donor_metabolic_error(void) const;
    inline double   get_donor_secretion_error(void) const;
    inline int32_t  get_donor_genome_size(void) const;
     
    void            set_indiv(Individual * indiv);
    void            set_id(int32_t id);
    void            set_rank(int32_t rank);
    inline void     set_parent_id(int32_t parent_id);
    inline void     set_parent_metabolic_error(double parent_metabolic_error);
    inline void     set_parent_secretion_error(double parent_secretion_error);
    inline void     set_parent_genome_size(int32_t parent_genome_size);
    inline void     set_donor_id(int32_t donor_id);
    inline void     set_donor_metabolic_error(double donor_metabolic_error);
    inline void     set_donor_secretion_error(double donor_secretion_error);
    inline void     set_donor_genome_size(int32_t donor_genome_size);
  
    inline const std::list<DnaReplicReport*> get_dna_replic_reports() const;
    inline void add_dna_replic_report(DnaReplicReport* rep) { _dna_replic_reports.push_back(rep); }

    // =================================================================
    //                            Public Methods
    // =================================================================
    void signal_end_of_replication(void);
    void write_to_tree_file(gzFile tree_file) const;


    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ReplicationReport(void)
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };



    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    Individual *  _indiv;
    
    int32_t         _id;
    int32_t         _rank;
    int32_t         _parent_id;
    
    int32_t         _genome_size;
    double          _metabolic_error;
    int16_t         _nb_genes_activ;
    int16_t         _nb_genes_inhib;
    int16_t         _nb_non_fun_genes;
    int16_t         _nb_coding_RNAs;
    int16_t         _nb_non_coding_RNAs;
    
    std::list<DnaReplicReport*>  _dna_replic_reports; // List of each genetic unit's replication report
    
    double          _parent_metabolic_error;
    double          _parent_secretion_error;
    int32_t         _parent_genome_size;
    
    int32_t			    _donor_id;
    double          _donor_metabolic_error;
    double          _donor_secretion_error;
    int32_t         _donor_genome_size;
    
    // CK: I think that the attributes below are obsolete 
    // (HT events are now stored in the ae_dna_replic_reports)
    
    double          _mean_align_score;
    
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline Individual *ReplicationReport::get_indiv(void) const
{
  return _indiv;
}

inline int32_t ReplicationReport::get_id(void) const
{
  return _id;
}

inline int32_t ReplicationReport::get_rank(void) const
{
  return _rank;
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

inline const std::list<DnaReplicReport*> ReplicationReport::get_dna_replic_reports() const {
  return _dna_replic_reports;
}

inline void ReplicationReport::set_indiv(Individual * indiv)
{
  _indiv = indiv;
}

inline void ReplicationReport::set_id(int32_t id)
{
  _id = id;
}

inline void ReplicationReport::set_rank(int32_t rank)
{
  _rank = rank;
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

#endif // __AE_REPLICATION_REPORT_H__
