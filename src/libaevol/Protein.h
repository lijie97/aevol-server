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
// ****************************************************************************


#ifndef AEVOL_PROTEIN_H_
#define AEVOL_PROTEIN_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "macros.h"
#include "Dna.h"
#include "Codon.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class Individual;
class Rna;





class Protein
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Protein() = delete;
    Protein(const Protein &model) = delete;
    Protein(GeneticUnit* gen_unit, const Protein &model);
    Protein(GeneticUnit* gen_unit,
               const std::list<Codon*>& codon_list,
               Strand strand,
               int32_t shine_dal_pos,
               Rna * rna,
               double w_max);
    //Protein(Protein* parent);
    Protein(gzFile backup_file);
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Protein();

    // =================================================================
    //                              Accessors
    // =================================================================
    inline Strand strand()                const;
    inline const std::list<Rna *> rna_list()          const;
    inline int32_t            shine_dal_pos()         const;
    inline int32_t            first_translated_pos()  const;
    inline int32_t            last_translated_pos()   const;
           int32_t            last_STOP_base_pos()    const;
    inline double             mean()                  const;
    inline double             width()                 const; // returns the half-width
    inline double             height()                const;
    inline int32_t            length()                const; // Number of Amino-Acids (not including START and STOP)
    inline double             concentration()         const;
    inline  bool              is_functional()         const;

    Individual * indiv() const;

    // =================================================================
    //                            Public Methods
    // =================================================================
    void  add_RNA(Rna * rna);
    char* AA_sequence(char separator = ' ') const; // WARNING : creates a new char[...] (up to you to delete it!)
    virtual void  save(gzFile backup_file);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    GeneticUnit*  gen_unit_;
    Strand strand_;
    std::list<Rna *>  rna_list_;              // RNAs transcribing this protein
    int32_t           shine_dal_pos_;         // Index of the corresponding shine dalgarno sequence in the genome
    int32_t           first_translated_pos_;  // Index of the first base following the START codon
    int32_t           last_translated_pos_;   // Index of the last base before the STOP codon
    int32_t           length_;                // Number of Amino-Acids (START and STOP codon do NOT produce AAs)
    double            concentration_;
    bool              is_functional_;

    std::list<Codon *> AA_list_;

    // Phenotypic contribution (triangle) parameters
    double mean_;
    double width_;   // in fact, half-width
    double height_;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
inline Strand Protein::strand() const
{
  return strand_;
}

inline const std::list<Rna *> Protein::rna_list() const {
  return rna_list_;
}

int32_t Protein::shine_dal_pos() const
{
  return shine_dal_pos_;
}

int32_t Protein::first_translated_pos() const
{
  return first_translated_pos_;
}

int32_t Protein::last_translated_pos() const
{
  return last_translated_pos_;
}

double Protein::mean() const
{
  return mean_;
}

double Protein::width() const
{
  return width_;
}

double Protein::height() const
{
  return height_;
}

int32_t Protein::length() const
{
  return length_;
}

double Protein::concentration() const
{
  return concentration_;
}

bool Protein::is_functional() const
{
  return is_functional_;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
} // namespace aevol
#endif // AEVOL_PROTEIN_H_
