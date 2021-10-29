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

      index_neutrino_sol = -1;

      index_had_t_b = -1;
      index_w_u = -1;
      index_w_d = -1;
      index_lep_t_b = -1;
  
      chi2 = 99999.;
      
      initial_had_t_mass = 99999;
      initial_had_w_mass = 99999;
      initial_lep_t_mass = 99999;
      initial_lep_w_mass = 99999;

      fitted_had_t_mass = 99999;
      fitted_had_w_mass = 99999;
      fitted_lep_t_mass = 99999;
      fitted_lep_w_mass = 99999;
  }//void Reset()

  int status;
  CUT_RESULT cut;

  int index_neutrino_sol;

  int index_had_t_b;
  int index_w_u;
  int index_w_d;
  int index_lep_t_b;

  vector<JET_ASSIGNMENT> vec_permutation;

  float chi2;

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
