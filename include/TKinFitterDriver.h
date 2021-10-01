#ifndef TKinFitterDriver_h
#define TKinFitterDriver_h

#include <TObject.h>

#include <Jet.h>
#include <Lepton.h>
#include <Particle.h>

#include <PhysicsTools/KinFitter/interface/TFitConstraintMGaus.h>
#include <PhysicsTools/KinFitter/interface/TAbsFitParticle.h>
#include <TFitParticlePt.h>
#include <PhysicsTools/KinFitter/interface/TFitParticleMCCart.h>

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

  //for easy permutation handling
  TLorentzVector jet_had_t_b, jet_w_u, jet_w_d, jet_lep_t_b, neutrino;

  //fitting objects and handlers
  TFitParticlePt* fit_lepton;
  unique_ptr<TFitParticlePt> u_fit_lepton;

  TMatrixD error_lepton;

  TVector3 neutrino_vec3;
  TFitParticleMCCart* fit_neutrino;
  unique_ptr<TFitParticleMCCart> u_fit_neutrino;

  TFitParticlePt* fit_jet_had_t_b;
  TFitParticlePt* fit_jet_w_u;
  TFitParticlePt* fit_jet_w_d;
  TFitParticlePt* fit_jet_lep_t_b;
  
  unique_ptr<TFitParticlePt> u_fit_jet_had_t_b;
  unique_ptr<TFitParticlePt> u_fit_jet_w_u;
  unique_ptr<TFitParticlePt> u_fit_jet_w_d;
  unique_ptr<TFitParticlePt> u_fit_jet_lep_t_b;

  //vector<TAbsFitParticle*> 

  TFitConstraintMGaus* constrain_had_t_mgaus;
  TFitConstraintMGaus* constrain_lep_t_mgaus;
  TFitConstraintMGaus* constrain_had_w_mgaus;
  TFitConstraintMGaus* constrain_lep_w_mgaus;

  unique_ptr<TFitConstraintMGaus> u_ptr_constrain_had_t_mgaus;

  bool BJet_Assignment_Cut();
  void Clear();
  void Fit();
  bool Pre_Kinematic_Cut();
  void Resol_Neutrino_Pz();//not used currently
  void Set_Constraints();
  void Set_Jets();
  void Set_Lepton();
  void Set_Neutrino(const int& index);
  void Sol_Neutrino_Pz();

  ClassDef(TKinFitterDriver, 1);
};

#endif /* TKinFitterDriver_h */
