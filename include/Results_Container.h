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
    
    best_del_phi_had_t_lep_t = 99999;

    best_theta_w_u_w_d = 99999;
    best_theta_had_w_had_t_b = 99999;
    best_theta_lep_neu = 99999;
    best_theta_lep_w_lep_t_b = 99999;

    best_chi2 = 99999;
    
    best_chi2_jet_had_t_b = 99999;
    best_chi2_jet_w_u = 99999;
    best_chi2_jet_w_d = 99999;
    best_chi2_jet_lep_t_b = 99999;
    
    best_chi2_jet_extra = 99999;

    best_chi2_constraint_had_t = 99999;
    best_chi2_constraint_had_w = 99999;
    best_chi2_constraint_lep_t = 99999;
    best_chi2_constraint_lep_w = 99999;
    
    best_mva_score = 99999;

    best_index_had_t_b = -1;
    best_index_w_u = -1;
    best_index_w_d = -1;
    best_index_lep_t_b = -1;
    
    best_initial_had_t_mass = 99999;
    best_initial_had_w_mass = 99999;
    best_initial_lep_t_mass = 99999;
    best_initial_lep_t_partial_mass = 99999;
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

  float best_del_phi_had_t_lep_t;

  float best_theta_w_u_w_d;
  float best_theta_had_w_had_t_b;
  float best_theta_lep_neu;
  float best_theta_lep_w_lep_t_b;

  float best_chi2;
  
  float best_chi2_jet_had_t_b;
  float best_chi2_jet_w_u;
  float best_chi2_jet_w_d;
  float best_chi2_jet_lep_t_b;
  
  float best_chi2_jet_extra;

  float best_chi2_constraint_had_t;
  float best_chi2_constraint_had_w;
  float best_chi2_constraint_lep_t;
  float best_chi2_constraint_lep_w;

  float best_mva_score;
  
  int best_index_had_t_b;
  int best_index_w_u;
  int best_index_w_d;
  int best_index_lep_t_b;

  float best_initial_had_t_mass;
  float best_initial_had_w_mass;
  float best_initial_lep_t_mass;
  float best_initial_lep_t_partial_mass;
  float best_initial_lep_w_mass;

  float best_fitted_had_t_mass;
  float best_fitted_had_w_mass;
  float best_fitted_lep_t_mass;
  float best_fitted_lep_w_mass;

  vector<Results> vec_results;
  
  ClassDef(Results_Container, 1);
};

#endif /* Results_Container_h */
