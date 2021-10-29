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
    neutrino_pz_sol[0] = -9999;
    neutrino_pz_sol[0] = -9999;
    chk_real_neu_pz = false;
   
    best_chi2 = 9999;
    best_had_t_mass = 9999;
    best_had_w_mass = 9999;
    best_lep_t_mass = 9999;
    best_lep_w_mass = 9999;


    vec_results.clear();
    vec_results.shrink_to_fit();
    
    return;
  }
  
  float neutrino_pz_sol[2];
  bool chk_real_neu_pz;
  
  float best_chi2;
  float best_had_t_mass;
  float best_had_w_mass;
  float best_lep_t_mass;
  float best_lep_w_mass;

  vector<Results> vec_results;
  
  ClassDef(Results_Container, 1);
};

#endif /* Results_Container_h */
