#ifndef Results_Container_h
#define Results_Container_h

#include <TObject.h>

#include <Results.h>

class Results_Container : public TObject
{
 public:
  Results_Container(){};
  ~Results_Container(){};
  
  void Reset()
  {
    status = false;

    best_neutrino_px = 99999;
    best_neutrino_py = 99999;
    best_neutrino_pz = 99999;
    
    best_chi2 = 99999;

    best_index_had_t_b = -1;
    best_index_w_u = -1;
    best_index_w_d = -1;
    best_index_lep_t_b = -1;
    
    best_initial_had_t_mass = 99999;
    best_initial_had_w_mass = 99999;
    best_initial_lep_t_mass = 99999;
    best_initial_lep_w_mass = 99999;
    
    best_fitted_had_t_mass = 99999;
    best_fitted_had_w_mass = 99999;
    best_fitted_lep_t_mass = 99999;
    best_fitted_lep_w_mass = 99999;

    vec_results.clear();
    vec_results.shrink_to_fit();
    
    return;
  }//void Reset
   
  bool status;

  float best_neutrino_px;
  float best_neutrino_py;
  float best_neutrino_pz;

  float best_chi2;
  
  int best_index_had_t_b;
  int best_index_w_u;
  int best_index_w_d;
  int best_index_lep_t_b;

  float best_initial_had_t_mass;
  float best_initial_had_w_mass;
  float best_initial_lep_t_mass;
  float best_initial_lep_w_mass;

  float best_fitted_had_t_mass;
  float best_fitted_had_w_mass;
  float best_fitted_lep_t_mass;
  float best_fitted_lep_w_mass;

  vector<Results> vec_results;
  
  ClassDef(Results_Container, 1);
};

#endif /* Results_Container_h */
