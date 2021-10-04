#ifndef Results_h
#define Results_h

#include <TObject.h>

class Results : public TObject
{
public:
  Results(){ Reset(); };
  ~Results(){}

  void Reset()
  {
      status = -1;
      cut = CUT_RESULT::PASS;

      index_had_t_b = -1;
      index_w_u = -1;
      index_w_d = -1;
      index_lep_t_b = -1;
  }

  enum CUT_RESULT {PASS, B_ASSIGNMENT, PREKINEMATIC, QUALITY};

  int status;
  CUT_RESULT cut;

  int index_had_t_b;
  int index_w_u;
  int index_w_d;
  int index_lep_t_b;
  
  ClassDef(Results, 1);
};

#endif /* Results_h */
