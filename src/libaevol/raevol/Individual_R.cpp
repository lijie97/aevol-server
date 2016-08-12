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

#include <unordered_map>
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

  _inherited_protein_list.resize(0);
}

Individual_R::Individual_R(const Individual_R& other)
    : Individual( other )
{
  _indiv_age = 0;
  _networked = false;
  _dist_sum = 0;


  if (exp_m_->exp_s()->get_with_heredity()) {
    _inherited_protein_list.clear();

    for (const auto& prot : other.protein_list_) {
      if (prot->concentration() >
          other.exp_m_->exp_s()->get_protein_presence_limit() && !((Protein_R*)prot)->is_signal()) {
        Protein_R* inherited_prot = new Protein_R(prot->get_gen_unit(),
                                                  (Protein_R&) *prot);
        inherited_prot->set_inherited(true);
        _inherited_protein_list.push_back(inherited_prot);
      }
    }
    //printf("Size of inherited is %d\n",_inherited_protein_list.size());
  }
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

  if (exp_m_->exp_s()->get_with_heredity()) {
    //_inherited_protein_list.resize(parent->protein_list_.size());
    //printf("%llu -- Parent protein size: %d\n",id, parent->protein_list_.size());

    for (const auto& prot : parent->protein_list_) {
      if (prot->concentration() >
          parent->exp_m_->exp_s()->get_protein_presence_limit() && ! ((Protein_R*)prot)->is_signal()) {
        Protein_R* inherited_prot = new Protein_R(prot->get_gen_unit(),
                                                  (Protein_R&) *prot, this->exp_m_);
        inherited_prot->set_inherited(true);
        _inherited_protein_list.push_back(inherited_prot);
        //printf("Add new herited protein\n");
      }
    }

    //printf("%llu -- Inherited protein size is %d\n",id,_inherited_protein_list.size());
  }
}

Individual_R::Individual_R(ExpManager* exp_m, gzFile backup_file) : Individual( exp_m, backup_file )
{
    _indiv_age = 0;
  _networked = false;

  if( exp_m_->exp_s()->get_with_heredity() )
  {
    // Retreive inherited proteins
    int32_t nb_inherited_proteins = 0;
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
Individual_R::~Individual_R( void ) noexcept
{

  if (exp_m_->exp_s()->get_with_heredity()) {
    for (unsigned int i = 0; i < _inherited_protein_list.size(); i++)
      delete _inherited_protein_list[i];

    _inherited_protein_list.clear();
  }

  for (unsigned int i = 0; i < _rna_list_coding.size(); i++) {
    _rna_list_coding[i] = NULL;
  }

  std::unordered_map<int,Protein_R * >::iterator iter = signal_list.begin();

  while (iter != signal_list.end())
  {
    Protein_R* clone = iter->second;
    delete clone;
    iter++;
  }

  signal_list.clear();

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
		EvaluateInContext(grid_cell_->habitat());
}

void Individual_R::EvaluateInContext(const Habitat_R& habitat, bool no_signal) {

    if (evaluated_ == true)
      return; // Individual has already been evaluated, nothing to do.
    for (const auto& prot : protein_list_) {
      ((Protein_R*) prot)->reset_concentration();
    }

    if (!_networked) {
      init_indiv(habitat);
    }
    /*for (const auto& prot : protein_list_) {
      printf("AT INIT %d ID %d Concentration of %d is %lf\n",AeTime::time(),id(),
             ((Protein_R*)prot)->get_id(),prot->concentration());
    }*/

    std::set<int>* eval = exp_m_->exp_s()->get_list_eval_step();
    // i is thus the age of the individual
    for (int8_t i = 1; i <= exp_m_->exp_s()->get_nb_indiv_age(); i++) {
      //Set the concentration of signals for this age
      for (auto prot1 : signal_list) {
        prot1.second->set_concentration(0.0);
      }
      if (no_signal)
        for (Protein_R* prot2 : habitat.phenotypic_target(i).signals()) {
          signal_list[prot2->get_id()]->set_concentration(0.9);
        }


      for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
        one_step();
      }


      /*for (const auto& prot : protein_list_) {
        printf("AT %d ID %d Concentration of %d is %lf\n",AeTime::time(),id(),
               ((Protein_R*)prot)->get_id(),prot->concentration());
      }*/

      // If we have to evaluate the individual at this age
      if (eval->find(i) != eval->end()) {
        //if (id_ % 1024 == 1) printf("Eval at %d\n",i);
        eval_step(habitat, i);
      }
    }


    final_step(habitat, exp_m_->exp_s()->get_nb_indiv_age());

//    protein_list_.clear();
//    protein_list_ = _initial_protein_list;



    //initialized = false;

}

void Individual_R::EvaluateInContext(const Habitat& habitat) {
  EvaluateInContext(dynamic_cast<const Habitat_R&> (habitat));
}

void Individual_R::init_indiv() {
  init_indiv(dynamic_cast<const Habitat_R&> (grid_cell_->habitat()));
}

void Individual_R::init_indiv(const Habitat_R& habitat)
{
  // ---------------------------------------------------------------------------
  // 1) Transcription - Translation - Folding - make_protein_list
  // ---------------------------------------------------------------------------
  transcribed_ = false;
  translated_ = false;
  folded_ = false;

  do_transcription_translation_folding();

  if (phenotype_ != NULL) {
    delete phenotype_;
    delete phenotype_activ_;
    delete phenotype_inhib_;

    phenotype_ = NULL;
    phenotype_activ_ = NULL;
    phenotype_inhib_ = NULL;
  }
  phenotype_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_activ_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_inhib_ = FuzzyFactory::fuzzyFactory->create_fuzzy();

  //----------------------------------------------------------------------------
  // 2) Make a list of all the rna present in the individual
  //    and initialise the concentrations of the proteins
  //----------------------------------------------------------------------------
  make_rna_list();

  int local_prot_id = 0;
  for (auto prot : protein_list_)
    ((Protein_R*)prot)->set_local_id(local_prot_id++);

  _initial_protein_list = protein_list_;

  //_protein_list.insert(_protein_list.end(), habitat.signals().begin(), habitat.signals().end());
  for(Protein_R* prot : habitat.signals()) {
    Protein_R* cloned = new Protein_R(prot);

    signal_list.insert({cloned->get_id(),cloned});
    protein_list_.push_back(cloned);
  }

  //----------------------------------------------------------------------------
  // 3) Create influence graph (including the signals)
  //----------------------------------------------------------------------------
  //printf("Protein %ld RNA %ld\n",protein_list_.size(),rna_list_.size());
  set_influences();

  _networked = true;
  //initialized = true;
}

void Individual_R::one_step( void )
{
  //----------------------------------------------------------------------------
  // 4) Make the individual "live its life" and compute partial phenotypes and
  //    fitnesses
  //----------------------------------------------------------------------------

  update_concentrations();
}

void Individual_R::eval_step( const Habitat_R& habitat, int8_t age ) {
  update_phenotype();
  distance_to_target_computed_ = false;
  phenotype_computed_ = true;

  for (int i=0; i<NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0;
  }

  compute_distance_to_target( habitat.phenotypic_target( age ) );

  /*printf("AT %d ID %d at %d dist is %lf on %d\n",AeTime::time(),id(),age,dist_to_target_by_feature_[METABOLISM],
         habitat.phenotypic_target( age ).get_id());*/

  _dist_sum += dist_to_target_by_feature_[METABOLISM];

  /*if (id_ % 1024 == 1)
    printf("Dist to target à l'age %d du nouveau clone : %f -- %f\n",
           age, dist_to_target_by_feature_[METABOLISM],_dist_sum);*/
}


void Individual_R::final_step( const Habitat_R& habitat, int8_t age ) {
  dist_to_target_by_feature_[METABOLISM] = _dist_sum / (double) (exp_m_->exp_s()->get_list_eval_step()->size());
  //printf("AT %d ID %d Meta error : %lf \n",AeTime::time(),id(),dist_to_target_by_feature_[METABOLISM]);


  fitness_computed_=false;
  // yoram attention il peut y avoir des soucis si on utilise des environnements segmentés ici
  compute_fitness(habitat.phenotypic_target( age ));

  /*if (id_ % 1024 == 1)
    printf("Nombre final d'évaluations : %d -- %f -- %e -- %e\n", exp_m_->exp_s()->get_list_eval_step()->size(),
           dist_to_target_by_feature_[METABOLISM], fitness_, fitness_by_feature_[METABOLISM]);*/

  phenotype_computed_ = true;
}

void Individual_R::set_influences()
// Compute the influence of each protein over each coding RNA
// As non-coding RNAs are completely inert, we don't care about their concentration
// so we don't care if proteins activate or inhibit their transcription.
{
  /*for (const auto& prot : protein_list_) {
   if (id() == 12608 && ((Protein_R*)prot)->is_signal())
      printf("ID 12608 set_influence Signal is %ld\n",((Protein_R*)prot)->get_id());
  }*/

	  for(auto& rna : _rna_list_coding) {
		  rna->set_influences( protein_list_ );
	  }
}

void Individual_R::update_concentrations( void )
{
	// Compute all the changes that will be applied to the concentrations
	// Concentrations must not be changed at this stage
  for (auto& prot : protein_list_) {
    if (!((Protein_R*)prot)->is_signal()) ((Protein_R*)prot)->compute_delta_concentration(exp_m_);
	}

	// Apply the changes in concentrations we have just computed
  for (auto& prot : protein_list_) {
		if (!((Protein_R*)prot)->is_signal()) ((Protein_R*)prot)->update_concentration();
	}
}

// Multiply the concentration of each protein by <factor>
void Individual_R::multiply_concentrations( double factor )
{
  for (auto& prot : protein_list_) {
	 	  ((Protein_R*)prot)->multiply_concentration( factor );
	}
}

/*int8_t Individual_R::get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  int8_t new_quadon = new_get_quadon(gen_unit,strand,pos);
  int8_t old_quadon = old_get_quadon(gen_unit,strand,pos);

  if (new_quadon != old_quadon) {
    printf("New/Old quadon are different %d %d\n",new_quadon,old_quadon);
    exit(-1);
  }
}*/

/*int8_t Individual_R::old_get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  const char* dna = gen_unit->dna()->data();
  int32_t  len    = gen_unit->dna()->length();
  int8_t quadon   = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[Utils::mod((pos+i),len)] == '1' )
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }
  else  // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {

      if ( dna[Utils::mod((pos-i),len)] != '1' ) // == and not != because we are on the complementary strand...
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }




  return quadon;
}*/

int8_t Individual_R::get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  const char* dna = gen_unit->dna()->data();
  int32_t  len    = gen_unit->dna()->length();
  int8_t quadon[4];

  if ( strand == LEADING )
  {
#ifdef __SIMD
#pragma omp simd
#endif
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      quadon[i] = (dna[((pos+i) % len < 0 ? (pos+i) % len + len : (pos+i) % len)] == '1') ? 1 << (QUADON_SIZE - i - 1) : 0;
    }
  }
  else  // ( strand == LAGGING )
  {
#ifdef __SIMD
#pragma omp simd
#endif
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      quadon[i] = (dna[((pos-i) % len < 0 ? (pos-i) % len + len : (pos-i) % len)] != '1') ? 1 << (QUADON_SIZE - i - 1) : 0;
    }
  }

  return quadon[0]+quadon[1]+quadon[2]+quadon[3];
}

void Individual_R::save( gzFile backup_file )
{
  Individual::save( backup_file );
  // Test if there is heredity, and if the generation is the first one (no inherited protein list).
  if (this->exp_m_->exp_s()->get_with_heredity() )
  {
    // Write inherited proteins
    int32_t nb_inherited_proteins = _inherited_protein_list.size();
    gzwrite( backup_file, &nb_inherited_proteins,  sizeof(int32_t) );

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

    if (this->exp_m_->exp_s()->get_with_heredity()) {
      for (auto& prot : _inherited_protein_list) {
        protein_list_.push_back(prot);
      }
    }
}

void Individual_R::make_rna_list( void )
{
  Individual::make_rna_list();
  _rna_list_coding = {};

  int local_rna_id = 0;

  // Parse the newly created RNA list and copy the coding RNAs in _rna_list_coding.
  for (const auto& gen_unit: genetic_unit_list_) {
    GeneticUnit* genu = const_cast<GeneticUnit*>(&gen_unit);
    // Create proxies
    const auto& rna_list = gen_unit.rna_list();
    const auto& lead = rna_list[LEADING];
    const auto& lagg = rna_list[LAGGING];

    // append pointers to rna material to local _rna_list
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: rna_list[strand]) {
        //TODO Ugly fix, change it to avoid memory usage double
        if (rna.is_coding()) {
          Rna_R* prna =  const_cast<Rna_R*>(&rna);
          //printf("COPY OR NOT : %ld == %ld",prna->get_id(),((Rna_R)rna).get_id());
          _rna_list_coding.push_back(
             prna);//new Rna_R(genu, rna));

          prna->set_local_id(local_rna_id++);
        }
    }

  }
}

void Individual_R::update_phenotype( void )
{
  // We will use two fuzzy sets :
  //   * _phenotype_activ for the proteins realising a set of functions
  //   * _phenotype_inhib for the proteins inhibitting a set of functions
  // The phenotype will then be given by the sum of these 2 fuzzy sets

  phenotype_activ_->reset();
  phenotype_inhib_->reset();
  phenotype_->reset();

  for (auto& prot : protein_list_) {
    if ( ((Protein_R*)prot)->is_functional() )
    {
      if ( ((Protein_R*)prot)->height() > 0 )
      {
//    	  added=true;
        phenotype_activ_->add_triangle(  ((Protein_R*)prot)->mean(),
                                         ((Protein_R*)prot)->width(),
                                         ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration() );

/*        printf("Add triangle ACTIV %f %f %f (%f %f)\n",((Protein_R*)prot)->mean(),
               ((Protein_R*)prot)->width(),
               ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(),
               ((Protein_R*)prot)->height(), ((Protein_R*)prot)->concentration() );*/
      }
      else
      {
        phenotype_inhib_->add_triangle(  ((Protein_R*)prot)->mean(),
                                         ((Protein_R*)prot)->width(),
                                         ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration() );

  /*      printf("Add triangle INHIB %f %f %f (%f %f)\n",((Protein_R*)prot)->mean(),
               ((Protein_R*)prot)->width(),
               ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(),
               ((Protein_R*)prot)->height(), ((Protein_R*)prot)->concentration() );*/
      }
    }
  }


    phenotype_activ_->clip(AbstractFuzzy::max,   Y_MAX);
    phenotype_inhib_->clip(AbstractFuzzy::min, - Y_MAX);

    phenotype_->add(*phenotype_activ_);
    phenotype_->add(*phenotype_inhib_);
    phenotype_->clip(AbstractFuzzy::min, Y_MIN);
    phenotype_->simplify();

//  _phenotype->simplify();
}

void    Individual_R::clear_everything_except_dna_and_promoters() {
    _networked = false;
    _rna_list_coding.clear();
    _dist_sum = 0.0;

    Individual::clear_everything_except_dna_and_promoters();
}

} // namespace aevol
