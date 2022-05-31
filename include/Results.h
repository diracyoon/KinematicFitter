#ifndef Results_h
#define Results_h

#include <TObject.h>

class TKinFitterDriver;

class Results : public TObject
{
public:
  Results(){ Reset(); };
  ~Results(){}

  void Reset()
  {
      status = -1;
      cut = CUT_RESULT::INIT;

      index_neutrino_sol = 99999;
      chk_real_neu_pz = false;

      met_px = 99999.;
      met_py = 99999.;
      met_rebalance_px = 99999.;
      met_rebalance_py = 99999.;
      neutrino_pz_sol = 99999.;

      index_had_t_b = -1;
      index_w_u = -1;
      index_w_d = -1;
      index_lep_t_b = -1;
  
      tt_delta_phi = 99999.;

      chi2 = 99999.;
  
      chi2_jet_had_t_b = 99999.;
      chi2_jet_w_u = 99999.;
      chi2_jet_w_d = 99999.;
      chi2_jet_lep_t_b = 99999.;

      chi2_jet_extra = 99999.;

      chi2_constraint_had_t = 99999.;
      chi2_constraint_had_w = 99999.;
      chi2_constraint_lep_t = 99999.;
      chi2_constraint_lep_w = 99999.;
    
      initial_had_t_mass = 99999.;
      initial_had_w_mass = 99999.;
      initial_lep_t_mass = 99999.;
      initial_lep_w_mass = 99999.;

      fitted_had_t_mass = 99999.;
      fitted_had_w_mass = 99999.;
      fitted_lep_t_mass = 99999.;
      fitted_lep_w_mass = 99999.;
  }//void Reset()

  int status;
  CUT_RESULT cut;

  int index_neutrino_sol;
  bool chk_real_neu_pz;

  float met_px;
  float met_py;
  float met_rebalance_px;
  float met_rebalance_py;
  float neutrino_pz_sol;

  int index_had_t_b;
  int index_w_u;
  int index_w_d;
  int index_lep_t_b;

  vector<JET_ASSIGNMENT> vec_permutation;

  float tt_delta_phi;//initial only

  float chi2;

  float chi2_jet_had_t_b;
  float chi2_jet_w_u;
  float chi2_jet_w_d;
  float chi2_jet_lep_t_b;
  
  float chi2_jet_extra;
  
  float chi2_constraint_had_t;
  float chi2_constraint_had_w;
  float chi2_constraint_lep_t;
  float chi2_constraint_lep_w;

  float initial_had_t_mass;
  float initial_had_w_mass;
  float initial_lep_t_mass;
  float initial_lep_w_mass;
  
  float fitted_had_t_mass;
  float fitted_had_w_mass;
  float fitted_lep_t_mass;
  float fitted_lep_w_mass;

  ClassDef(Results, 1);
};

#endif /* Results_h */
