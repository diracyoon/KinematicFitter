#ifndef TKinFitterDriver_h
#define TKinFitterDriver_h

#include <TObject.h>

#include <Jet.h>
#include <Lepton.h>
#include <Particle.h>

//#include <TFitParticlePt.h>

#define W_MASS 80.379

using namespace std;

class TKinFitterDriver : public TObject
{
 public:
  TKinFitterDriver(){};
  TKinFitterDriver(int _data_year);

  bool End_Permutation();
  void Scan();
  void Set_Objects(vector<Jet>& _vec_jet, vector<bool>& _vec, Lepton& _lepton, Particle& _met);

  //name convention: mother_species
  enum JET_ASSIGNMENT {HAD_T_B, W_U, W_D, LEP_T_B, OTHERS};

 protected:
  int data_year;
  
  vector<Jet> vec_jet;
  vector<bool> vec_btag;
  Lepton lepton;
  Particle met;

  vector<TKinFitterDriver::JET_ASSIGNMENT> vec_permutation;
  int n_jet, n_btag;

  float neutrino_pz_sol[2];
  bool chk_real_neu_pz;

  TLorentzVector jet_had_t_b, jet_w_u, jet_w_d, jet_lep_t_b, neutrino;

  //TFitParticlePt* fit_lepton;

  bool BJet_Assignment_Cut();
  void Clear();
  void Fit();
  bool Pre_Kinematic_Cut();
  void Resol_Neutrino_Pz();//not used currently
  void Set_Current_Permutation();
  void Set_Lepton();
  void Set_Neutrino(const int& index);
  void Sol_Neutrino_Pz();

  ClassDef(TKinFitterDriver, 1);
};

#endif /* TKinFitterDriver_h */
