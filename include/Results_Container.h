#ifndef Results_Container_h
#define Results_Container_h

#include <TObject.h>

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
   
    return;
  }
  
  float neutrino_pz_sol[2];
  bool chk_real_neu_pz;
  
  //vector<Results> vec_results;
  
  ClassDef(Results_Container, 1);
};

#endif /* Results_Container_h */
