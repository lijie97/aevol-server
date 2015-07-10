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




// =================================================================
//                              Libraries
// =================================================================

// =================================================================
//                            Project Files
// =================================================================
#include "Individual_R.h"
#include "../ExpManager.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_individual_R                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

/*
 * Used at initialization
*/
Individual_R::Individual_R(ExpManager* exp_m,
                       std::shared_ptr<JumpingMT> mut_prng,
                       std::shared_ptr<JumpingMT> stoch_prng,
                       std::shared_ptr<MutationParams> param_mut,
                       double w_max,
                       int32_t min_genome_length,
                       int32_t max_genome_length,
                       bool allow_plasmids,
                       int32_t id,
                       const char* strain_name,
                       int32_t age) : Individual(exp_m,mut_prng,stoch_prng,param_mut,w_max,min_genome_length,
                                                 max_genome_length,allow_plasmids,id,strain_name,age) {

  _indiv_age = 0;
  _networked = false;
  _dist_sum = 0;
}

Individual_R::Individual_R(const Individual_R& other)
    : Individual( other )
{
  _indiv_age = 0;
  _networked = false;
  _dist_sum = 0;
  _inherited_protein_list = std::vector<Protein_R*>(other._inherited_protein_list);
}

Individual_R::Individual_R( Individual_R* parent, int32_t id,
                            std::shared_ptr<JumpingMT> mut_prng,
                            std::shared_ptr<JumpingMT> stoch_prng)
        : Individual( parent, id, mut_prng, stoch_prng )
{
  //~ printf( "ae_individual_R( parent ) : I have %d inherited proteins\n", parent->get_protein_list()->get_nb_elts() );
    _indiv_age = 0;
    _networked = false;
    _dist_sum = 0;

    for (const auto& prot : parent->_protein_list)
    {
    	if( prot->get_concentration() > parent->get_exp_m()->get_exp_s()->get_protein_presence_limit() )
    	{
    		Protein_R* inherited_prot = new Protein_R( prot->get_gen_unit(),
                                                   (Protein_R&)*prot);
    		inherited_prot->set_inherited( true );
    		_inherited_protein_list.push_back( inherited_prot );
    	}
    }
}

Individual_R::Individual_R(ExpManager* exp_m, gzFile backup_file) : Individual( exp_m, backup_file )
{
    _indiv_age = 0;
  if( get_exp_m()->get_exp_s()->get_with_heredity() )
  {
    // Retreive inherited proteins
//    _inherited_protein_list = new ae_list();
    int16_t nb_inherited_proteins = 0;
    gzread( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );
  
    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
	  _inherited_protein_list.push_back( new Protein_R( backup_file ) );
    }
  }  
}

// =================================================================
//                             Destructors
// =================================================================
Individual_R::~Individual_R( void )
{
  assert( !get_exp_m()->get_exp_s()->get_with_heredity()  );

  /*
  for (const auto& prot : parent->_protein_list)
	  Protein_R* dp = _inherited_protein_list[i];
	  delete dp;
  }*/

  _inherited_protein_list.clear();

  _rna_list_coding.clear();
}

// =================================================================
//                            Public Methods
// =================================================================
Individual_R* Individual_R::CreateClone(const Individual_R* dolly, int32_t id) {
  Individual_R* indiv = new Individual_R(*dolly);
  indiv->set_id(id);
  return indiv;
}

void Individual_R::Evaluate() {
		EvaluateInContext(_grid_cell->habitat());
}

void Individual_R::EvaluateInContext(const Habitat& habitat) {
	if (_evaluated == true) return; // Individual has already been evaluated, nothing to do.

    if (!_networked) {
      // ---------------------------------------------------------------------------
      // 1) Transcription - Translation - Folding - make_protein_list
      // ---------------------------------------------------------------------------

      _transcribed = false;
      _translated = false;
      _folded = false;

      do_transcription_translation_folding();

      if (_phenotype != NULL) {
        delete _phenotype;
        delete _phenotype_activ;
        delete _phenotype_inhib;

        _phenotype = NULL;
        _phenotype_activ = NULL;
        _phenotype_inhib = NULL;
      }
      _phenotype = new Phenotype();
      _phenotype_activ = new Phenotype();
      _phenotype_inhib = new Phenotype();

      //----------------------------------------------------------------------------
      // 2) Make a list of all the rna present in the individual
      //    and initialise the concentrations of the proteins
      //----------------------------------------------------------------------------
#ifdef __TRACING__
	    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	    #endif

      make_rna_list();

      for (const auto& prot : _protein_list) {
        ((Protein_R*) prot)->set_initial_concentration();
      }

#ifdef __TRACING__
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          ae_logger::addLog(MAKERNALIST,duration);
          t1 = t2;
        #endif
      //----------------------------------------------------------------------------
      // 3) Create influence graph (including the signals)
      //----------------------------------------------------------------------------
      set_influences();

      _networked = true;
#ifdef __TRACING__
	    t2 = high_resolution_clock::now();
	  duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	  ae_logger::addLog(SETINFLUS,duration);
	  t1 = t2;
	    #endif
    }
    //----------------------------------------------------------------------------
	// 4) Make the individual "live its life" and compute partial phenotypes and
	//    fitnesses
	//----------------------------------------------------------------------------

	#ifdef __TRACING__
		  t1 = high_resolution_clock::now();
	#endif

  update_concentrations();

	#ifdef __TRACING__
	    t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(UPDATECONCENT,duration);
	    t1 = t2;
	#endif
  if (std::fmod(((double)Time::get_time()), get_exp_m()->get_exp_s()->get_eval_step()) == 0.0)
	{
    update_phenotype();
	  _distance_to_target_computed = false;
	  _phenotype_computed = true;
	  compute_distance_to_target( habitat.phenotypic_target() );
    _dist_sum += _dist_to_target_by_feature[METABOLISM];
	}

	#ifdef __TRACING__
	    t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(UPDATEPHENO,duration);
	    t1 = t2;
	#endif


  //----------------------------------------------------------------------------
	// 5) Compute final fitness and final dist to target
	// ----------------------------------------------------------------------------

  if (std::fmod(((double)Time::get_time()), get_exp_m()->get_exp_s()->get_nb_indiv_age()) == 0.0)
  {
    // On devrait faire la somme du carré des erreurs afin d'éviter qu'elles puissent se compenser
    _dist_to_target_by_feature[METABOLISM] = _dist_sum / (double) (get_exp_m()->get_exp_s()->get_nb_indiv_age() / get_exp_m()->get_exp_s()->get_eval_step());
    _fitness_computed=false;
    // yoram attention il peut y avoir des soucis si on utilise des environnements segmentés ici
    compute_fitness(habitat.phenotypic_target());
    _phenotype_computed = true;
  }

	#ifdef __TRACING__
	  t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(FITNESS,duration);
	    t1 = t2;
	#endif
}

void Individual_R::set_influences( void )
// Compute the influence of each protein over each coding RNA
// As non-coding RNAs are completely inert, we don't care about their concentration
// so we don't care if proteins activate or inhibit their transcription.
{
	  for(auto& rna : _rna_list_coding) {
		  rna->set_influences( _protein_list );
	  }
}

void Individual_R::update_concentrations( void )
{
	// Compute all the changes that will be applied to the concentrations
	// Concentrations must not be changed at this stage
  for (auto& prot : _protein_list) {
		if (!((Protein_R*)prot)->is_signal()) ((Protein_R*)prot)->compute_delta_concentration();
	}

	// Apply the changes in concentrations we have just computed
  for (auto& prot : _protein_list) {
		if (!((Protein_R*)prot)->is_signal()) ((Protein_R*)prot)->update_concentration();
	}
}

// Multiply the concentration of each protein by <factor>
void Individual_R::multiply_concentrations( double factor )
{
  for (auto& prot : _protein_list) {
	 	  ((Protein_R*)prot)->multiply_concentration( factor );
	}
}

int8_t Individual_R::get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  const char* dna = gen_unit->get_dna()->get_data();
  int32_t  len    = gen_unit->get_dna()->get_length();
  int8_t quadon   = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[(pos+i) % len] == '1' )
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }
  else  // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[(pos-i) % len] != '1' ) // == and not != because we are on the complementary strand...
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }

  return quadon;
}

void Individual_R::save( gzFile backup_file )
{
  Individual::save( backup_file );
  // Test if there is heredity, and if the generation is the first one (no inherited protein list).
  if (this->get_exp_m()->get_exp_s()->get_with_heredity() && !_inherited_protein_list.empty() )
  {
    // Write inherited proteins
    int16_t nb_inherited_proteins = _inherited_protein_list.size();
    gzwrite( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );

    for (auto& prot : _inherited_protein_list) {
    	prot->save( backup_file );
    }
  }
}
// =================================================================
//                           Protected Methods
// =================================================================
void Individual_R::make_protein_list( void )
{
	  Individual::make_protein_list();

	  for (auto& prot : _inherited_protein_list)
		  _protein_list.push_back( prot );
}

void Individual_R::make_rna_list( void )
{
  Individual::make_rna_list();
  _rna_list_coding.clear();
  _rna_list_coding.reserve(_rna_list.size());

  // Parse the newly created RNA list and copy the coding RNAs in _rna_list_coding.
  for (const auto& gen_unit: _genetic_unit_list) {
    // Create proxies
    const auto& rna_list = gen_unit.get_rna_list();
    const auto& lead = rna_list[LEADING];
    const auto& lagg = rna_list[LAGGING];

    // append pointers to rna material to local _rna_list
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: rna_list[strand]) {
        //TODO Ugly fix, change it to avoid memory usage double
        _rna_list_coding.push_back(new Rna_R(const_cast<GeneticUnit*>(&gen_unit), rna));
    }
  }
}

void Individual_R::update_phenotype( void )
{
  // We will use two fuzzy sets :
  //   * _phenotype_activ for the proteins realising a set of functions
  //   * _phenotype_inhib for the proteins inhibitting a set of functions
  // The phenotype will then be given by the sum of these 2 fuzzy sets

  _phenotype_activ->reset();
  _phenotype_inhib->reset();
  _phenotype->reset();

  for (auto& prot : _protein_list) {
    if ( ((Protein_R*)prot)->get_is_functional() )
    {
      if ( ((Protein_R*)prot)->get_height() > 0 )
      {
//    	  added=true;
        _phenotype_activ->add_triangle(  ((Protein_R*)prot)->get_mean(),
                                         ((Protein_R*)prot)->get_width(),
                                         ((Protein_R*)prot)->get_height() * ((Protein_R*)prot)->get_concentration() );
      }
      else
      {
        _phenotype_inhib->add_triangle(  ((Protein_R*)prot)->get_mean(),
                                         ((Protein_R*)prot)->get_width(),
                                         ((Protein_R*)prot)->get_height() * ((Protein_R*)prot)->get_concentration() );
      }
    }
  }

  _phenotype_activ->clip(Fuzzy::max, Y_MAX );
  _phenotype_inhib->clip(Fuzzy::min,  -Y_MAX );

  _phenotype->add( *_phenotype_activ );
  _phenotype->add( *_phenotype_inhib );
  _phenotype->clip(Fuzzy::min,  Y_MIN );

//  if (added) {printf("PHENO: \n");_phenotype->print_points();}
//  _phenotype->simplify();
}
} // namespace aevol
