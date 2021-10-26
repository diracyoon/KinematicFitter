#ifndef TKinFitterDriver_h
#define TKinFitterDriver_h

#include <TObject.h>

#include "Jet.h"
#include "Lepton.h"
#include "Particle.h"
#include "Muon.h"
#include "Photon.h"
#include "Electron.h"
#include "Event.h"
#include "Gen.h"
#include "LHE.h"
#include <vector>

//in /cvmfs/cms.cern.ch/$(SCRAM_ARCH)/cms/$(cmsswrel)/src/PhysicsTools/KinFitter/interface
#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleMCCart.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintMGaus.h"

//Redefined by B. Oh
#include "TFitParticlePt.h"

#include "Results_Container.h"
#include "Results.h"

#define T_MASS 172.5
#define T_WIDTH 1.5
#define W_MASS 80.379
#define W_WIDTH 2.085

using namespace std;

class TKinFitterDriver : public TObject
{
 public:
  TKinFitterDriver(){};
  TKinFitterDriver(int _data_year); 
  
  void Scan();  
  void Set_Objects(vector<Jet>& _vec_jet, vector<float>& _vec_resolution_pt, vector<bool>& _vec_btag, Lepton& _lepton, Particle& _met);  
  
  //name convention: mother_species  
  enum JET_ASSIGNMENT {HAD_T_B, W_U, W_D, LEP_T_B, OTHERS};  
  
 protected:  
  int data_year;  
  
  vector<Jet> vec_jet;  
  vector<float> vec_resolution_pt;
  vector<bool> vec_btag;
  Lepton lepton;
  Particle met;

  vector<TKinFitterDriver::JET_ASSIGNMENT> vec_permutation;  
  int n_jet, n_btag;  
  
  //for easy permutation handling 
  TLorentzVector jet_had_t_b, jet_w_u, jet_w_d, jet_lep_t_b, neutrino; 
  vector<TLorentzVector> vec_jet_extra; 
  
  //fitting objects and handlers  
  TFitParticlePt* fit_lepton; 
  unique_ptr<TFitParticlePt> u_fit_lepton; 
  TMatrixD error_lepton; 

  TVector3 neutrino_vec3; 
  TFitParticleMCCart* fit_neutrino;
  unique_ptr<TFitParticleMCCart> u_fit_neutrino;
  //TMatrixD error_neutrino;

  TFitParticlePt* fit_jet_had_t_b;
  TFitParticlePt* fit_jet_w_u;
  TFitParticlePt* fit_jet_w_d;
  TFitParticlePt* fit_jet_lep_t_b;
  unique_ptr<TFitParticlePt> u_fit_jet_had_t_b;
  unique_ptr<TFitParticlePt> u_fit_jet_w_u;
  unique_ptr<TFitParticlePt> u_fit_jet_w_d;
  unique_ptr<TFitParticlePt> u_fit_jet_lep_t_b;
  TMatrixD error_jet_had_t_b;
  TMatrixD error_jet_w_u;
  TMatrixD error_jet_w_d;
  TMatrixD error_jet_lep_t_b;
  
  vector<TFitParticlePt*> vec_fit_jet_extra;
  vector<unique_ptr<TFitParticlePt>> vec_u_fit_jet_extra;
  vector<TMatrixD> vec_error_jet_extra;
  
  TFitConstraintMGaus* constraint_had_t_mgaus;
  TFitConstraintMGaus* constraint_lep_t_mgaus;
  TFitConstraintMGaus* constraint_had_w_mgaus;
  TFitConstraintMGaus* constraint_lep_w_mgaus;
  
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_had_t_mgaus;
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_lep_t_mgaus;
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_had_w_mgaus;
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_lep_w_mgaus;

  TFitConstraintEp* constraint_px;
  TFitConstraintEp* constraint_py;

  unique_ptr<TFitConstraintEp> u_ptr_constraint_px;
  unique_ptr<TFitConstraintEp> u_ptr_constraint_py;

  TKinFitter* fitter;
  unique_ptr<TKinFitter> u_ptr_fitter;

  Results results;
  Results_Container results_container;

  bool BJet_Assignment_Cut(); 
  float Calc_Chi2();
  float Calc_Each_Chi2(TAbsFitParticle* ptr);
  bool Check_Repetition();
  void Clear();
  bool End_Permutation(){ return next_permutation(vec_permutation.begin(), vec_permutation.end()); } 
  bool Pre_Kinematic_Cut(); 
  bool Quality_Cut();
  void Resol_Neutrino_Pz();//not used currently 
  void Save_Permutation(const bool& push=false);
  void Save_Results();
  void Set_Constraints();
  void Set_Fitter();
  void Set_Jets();  
  void Set_Lepton(); 
  void Set_Neutrino(const int& index); 
  void Sol_Neutrino_Pz(); 

  ClassDef(TKinFitterDriver, 1);
};

#endif /* TKinFitterDriver_h */
