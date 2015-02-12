//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/** \class
 *  \brief
 */


// =================================================================
//                              Libraries
// =================================================================
#include <assert.h>



// =================================================================
//                            Project Files
// =================================================================

#include "ae_mutation.h"
#include "ae_gene_mutation.h"
#include "ae_rna.h"
#include "ae_macros.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                        Class ae_gene_mutation                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================



// =================================================================
//                             Constructors
// =================================================================




// Creates a copy of the mutation mut, but enriched with the generation when it occured
// and the position where it occurred in the RNA, relative to the first bp of the promoter
ae_gene_mutation::ae_gene_mutation(ae_mutation const & mut, int32_t gener, int32_t cdsPosBefore, ae_strand strandBefore, ae_gene_mutation_region region  ) : ae_mutation(mut)
{
  _generation = gener;
  _impact_on_metabolic_error = 0.0; /* should be set to its real value when known */
  _region = region;

  /* Compute _position_relative_to_shine_dal */

  switch ( _mut_type )
    {
    case SWITCH :
      _position_relative_to_shine_dal = new int32_t;
      if ( strandBefore == LEADING  ) {_position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;}
      else                              {_position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];}
      break;
    case S_INS :
      _position_relative_to_shine_dal = new int32_t;
      if ( strandBefore == LEADING  ) {_position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;}
      else                              {_position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];}
      break;
    case S_DEL :
      _position_relative_to_shine_dal = new int32_t;
      if ( strandBefore == LEADING  ) {_position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;}
      else                              {_position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];}
      break;
    case DUPL : 
      /* A duplication can affect a gene in two ways:
         1) The reinsertion point of the duplicated segment is located within the gene => stored
         2) The gene is partly or completely duplicated, but this does not change its sequence => nothing to store       
      */
      /* We should enter here only in case (1). Note that in this case, the relative positions for beginseg and endseg may be outside the gene */ 
      _position_relative_to_shine_dal = new int32_t[3];
      if ( strandBefore == LEADING  )
        {
          _position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;
          _position_relative_to_shine_dal[1] = _pos[1] - cdsPosBefore;
          _position_relative_to_shine_dal[2] = _pos[2] - cdsPosBefore;
        }
      else    
        {
          _position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];
          _position_relative_to_shine_dal[1] = cdsPosBefore - _pos[1];
          _position_relative_to_shine_dal[2] = cdsPosBefore - _pos[2];
        }
      break;
    case DEL : 
      _position_relative_to_shine_dal = new int32_t[2];
      if ( strandBefore == LEADING  )
        {
          _position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;
          _position_relative_to_shine_dal[1] = _pos[1] - cdsPosBefore;
        }
      else    
        {
          _position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];
          _position_relative_to_shine_dal[1] = cdsPosBefore - _pos[1];
        }
      break;
    case TRANS :
      _position_relative_to_shine_dal = new int32_t[4];
      if ( strandBefore == LEADING  )
        {
          _position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;
          _position_relative_to_shine_dal[1] = _pos[1] - cdsPosBefore;
          _position_relative_to_shine_dal[2] = _pos[2] - cdsPosBefore;
          _position_relative_to_shine_dal[3] = _pos[3] - cdsPosBefore;
        }
      else    
        {
          _position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];
          _position_relative_to_shine_dal[1] = cdsPosBefore - _pos[1];
          _position_relative_to_shine_dal[2] = cdsPosBefore - _pos[2];
          _position_relative_to_shine_dal[3] = cdsPosBefore - _pos[3];
        }
      break;
    case INV :
      _position_relative_to_shine_dal = new int32_t[2];
      if ( strandBefore == LEADING  )
        {
          _position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;
          _position_relative_to_shine_dal[1] = _pos[1] - cdsPosBefore;
        }
      else    
        {
          _position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];
          _position_relative_to_shine_dal[1] = cdsPosBefore - _pos[1];
        }
      break;
    case INSERT :
      _position_relative_to_shine_dal = new int32_t;
      if ( strandBefore == LEADING  ) _position_relative_to_shine_dal[0] = _pos[0] - cdsPosBefore;
      else                                 _position_relative_to_shine_dal[0] = cdsPosBefore - _pos[0];
      break;
    default :
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
}


// =================================================================
//                             Destructors
// =================================================================

ae_gene_mutation::~ae_gene_mutation()
{
  /* ae_mutation::~ae_mutation() will be called automatically by the compiler for the other attributes */
  switch ( _mut_type )
  {
  case SWITCH :
    delete _position_relative_to_shine_dal;
    break;
  case S_INS :
    delete _position_relative_to_shine_dal;
    break;
  case S_DEL :
    delete _position_relative_to_shine_dal;
    break;
  case DUPL :
    delete [] _position_relative_to_shine_dal;
    break;
  case DEL :
    delete [] _position_relative_to_shine_dal;
    break;
  case TRANS :
    delete [] _position_relative_to_shine_dal;
    break;
  case INV :
      delete [] _position_relative_to_shine_dal;
      break;
  case INSERT :
    delete _position_relative_to_shine_dal;
      break;
  default :
    fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n" , _mut_type, __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
    break;
  }
  
} ;




// =================================================================
//                        Public methods
// =================================================================

// 0 if local mut, 1 if rearrangement, 2 if transfer
int8_t ae_gene_mutation::type_of_event()
{
  if ((_mut_type == SWITCH) || (_mut_type == S_INS) || (_mut_type == S_DEL)) return 0;
  else if ((_mut_type == DUPL) || (_mut_type == DEL) || (_mut_type == TRANS) || (_mut_type == INV)) return 1;
  else return 2;
    

}



// str must be at least of size 60
void ae_gene_mutation::get_description_string_for_gene_mut(char * str)
{ 
   switch ( _mut_type )
  {
    case SWITCH :
    {
      sprintf( str, "%" PRId32 " SWITCH %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _impact_on_metabolic_error );
      break;
    }
    case S_INS :
    {
      sprintf( str, "%" PRId32 " SMALL_INS %" PRId32 " %" PRId32 " %s %.10f " , _generation, _position_relative_to_shine_dal[0], _length[0], _seq, _impact_on_metabolic_error );
      break;
    }
    case S_DEL :
    {
      sprintf( str, "%" PRId32 " SMALL_DEL %" PRId32 " %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _length[0], _impact_on_metabolic_error );
      break;
    }
    case DUPL :
    {
      sprintf( str, "%" PRId32 " INSERTION_OF_DUPLICATED_DNA %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _position_relative_to_shine_dal[1], _position_relative_to_shine_dal[2], _length[0], _impact_on_metabolic_error);
      break;
    }
    case DEL :
    {
      sprintf( str, "%" PRId32 " LARGE_DEL %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _position_relative_to_shine_dal[1], _length[0], _impact_on_metabolic_error );
      break;
    }
    case TRANS :
    {
      sprintf( str, "%" PRId32 " TRANSLOC %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _position_relative_to_shine_dal[1],  _position_relative_to_shine_dal[2],  _position_relative_to_shine_dal[3], _length[0], _impact_on_metabolic_error );    
      break;
    }
    case INV :
    {
      sprintf( str, "%" PRId32 " INV %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _position_relative_to_shine_dal[1], _length[0], _impact_on_metabolic_error );
      break;
    }
   case INSERT :
    {
      sprintf( str, "%" PRId32 " INSERTION_OF_FOREIGN_DNA %" PRId32 " %" PRId32 " %.10f " , _generation, _position_relative_to_shine_dal[0], _length[0], _impact_on_metabolic_error );
      break;
    }
    default :
    {
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n" , _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
  }

   if (_region == CDS)            strcat(str, "CDS ");
   else if (_region == UPSTREAM)  strcat(str, "UPSTREAM ");
   else if (_region == BOTH)      strcat(str, "BOTH ");


}

} // namespace aevol
