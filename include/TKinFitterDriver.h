#ifndef TKinFitterDriver_h
#define TKinFitterDriver_h

#include <TObject.h>
#include <TMVA/Reader.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>

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

// in /cvmfs/cms.cern.ch/$(SCRAM_ARCH)/cms/$(cmsswrel)/src/PhysicsTools/KinFitter/interface
#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleMCCart.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintMGaus.h"

// #include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "TKinFitter_Mod.h"

// Redefined by B. Oh
#include "TFitParticlePt.h"

#include "Enum_Def.h"
#include "Results_Container.h"
#include "Results.h"
#include "../../../Analyzers/include/Vcb_Def.h"

using namespace std;

class TKinFitterDriver : public TObject
{
public:
  TKinFitterDriver(){};
  TKinFitterDriver(const TString &_data_era, const TString &_channel, bool _run_permutatation_tree = false, bool _run_chi2 = false, bool _rm_wm_constraint = false, bool _rm_bjet_energy_reg_nn = false);

  bool Check_Status() { return results_container.status; }
  Results_Container Get_Results() { return results_container; }
  void Scan();
  void Set_Objects(vector<Jet> &_vec_jet, vector<float> &_vec_resolution_pt, vector<bool> &_vec_btag, Lepton &_lepton, Particle &_met, bool _chk_matched = false, int *_index_matched_jet = NULL);

protected:
  TString data_era;
  TString channel;

  bool chk_bvsc_only = false;
  bool chk_two_step_reco = true;

  bool run_permutation_tree;
  bool bjet_energy_reg_nn;
  bool run_chi2;
  bool rm_wm_constraint;

  vector<Jet> vec_jet;
  vector<float> vec_resolution_pt;
  vector<bool> vec_btag;
  Lepton lepton;
  Particle met;
  Particle met_rebalance;

  bool chk_matched;
  int index_matched_jet[4];
  int permutation_to_index[4];

  vector<JET_ASSIGNMENT> vec_permutation;
  int n_jet;
  int n_btag;

  int correct_permutation[4];

  float neutrino_pz_sol[2];
  bool chk_real_neu_pz;

  int index_neutrino_sol;

  float neutrino_px_init;
  float neutrino_py_init;
  float neutrino_pz_init;

  // for easy permutation handling
  Jet jet_had_t_b, jet_w_u, jet_w_d, jet_lep_t_b, neutrino;
  vector<Jet> vec_jet_extra;

  // fitting objects and handlers
  TFitParticlePt *fit_lepton;
  unique_ptr<TFitParticlePt> u_fit_lepton;
  TMatrixD error_lepton;

  TVector3 neutrino_vec3;
  TFitParticleMCCart *fit_neutrino;
  unique_ptr<TFitParticleMCCart> u_fit_neutrino;
  // TMatrixD error_neutrino;

  TFitParticlePt *fit_jet_had_t_b;
  TFitParticlePt *fit_jet_w_u;
  TFitParticlePt *fit_jet_w_d;
  TFitParticlePt *fit_jet_lep_t_b;
  unique_ptr<TFitParticlePt> u_fit_jet_had_t_b;
  unique_ptr<TFitParticlePt> u_fit_jet_w_u;
  unique_ptr<TFitParticlePt> u_fit_jet_w_d;
  unique_ptr<TFitParticlePt> u_fit_jet_lep_t_b;
  TMatrixD error_jet_had_t_b;
  TMatrixD error_jet_w_u;
  TMatrixD error_jet_w_d;
  TMatrixD error_jet_lep_t_b;

  vector<TFitParticlePt *> vec_fit_jet_extra;
  vector<unique_ptr<TFitParticlePt>> vec_u_fit_jet_extra;
  vector<TMatrixD> vec_error_jet_extra;

  TFitConstraintMGaus *constraint_had_t_mgaus;
  TFitConstraintMGaus *constraint_lep_t_mgaus;
  TFitConstraintMGaus *constraint_had_w_mgaus;
  TFitConstraintMGaus *constraint_lep_w_mgaus;

  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_had_t_mgaus;
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_lep_t_mgaus;
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_had_w_mgaus;
  unique_ptr<TFitConstraintMGaus> u_ptr_constraint_lep_w_mgaus;

  TFitConstraintEp *constraint_px;
  TFitConstraintEp *constraint_py;

  unique_ptr<TFitConstraintEp> u_ptr_constraint_px;
  unique_ptr<TFitConstraintEp> u_ptr_constraint_py;

  // TKinFitter* fitter;
  // unique_ptr<TKinFitter> u_ptr_fitter;
  TKinFitter_Mod *fitter;
  unique_ptr<TKinFitter_Mod> u_ptr_fitter;

  Results results;
  Results_Container results_container;

  vector<float> vec_prekin_score;

  TMVA::Reader *reader_permutation_step_0[2];
  TMVA::Reader *reader_permutation_step_1[2];
  TMVA::Reader *reader_prekin_cut;

  float met_pt;
  float neutrino_p;
  float lepton_pt;
  float pt_ratio;

  float pt_had_t_b;
  float pt_w_u;
  float pt_w_d;
  float pt_lep_t_b;

  float bvsc_had_t_b;
  float cvsb_had_t_b;
  float cvsl_had_t_b;

  float bvsc_w_u;
  float cvsb_w_u;
  float cvsl_w_u;

  float bvsc_w_d;
  float cvsb_w_d;
  float cvsl_w_d;

  float bvsc_lep_t_b;
  float cvsb_lep_t_b;
  float cvsl_lep_t_b;

  float theta_w_u_w_d;
  float theta_had_w_had_t_b;
  float theta_lep_neu;
  float theta_lep_w_lep_t_b;
  float del_phi_had_t_lep_t;
  float had_t_mass;
  float had_w_mass;
  float lep_t_mass;
  float lep_t_partial_mass;
  float chi2_jet_had_t_b;
  float chi2_jet_w_u;
  float chi2_jet_w_d;
  float chi2_jet_lep_t_b;
  float chi2_jet_extra;
  float chi2_constraint_had_t;
  float chi2_constraint_had_w;
  float chi2_constraint_lep_t;
  float chi2_constraint_lep_w;
  float chi2;

  float prekin_cut_mva_score[4];

  void Remove_Ambiguity();
  bool BJet_Assignment_Cut();
  float Calc_Chi2();
  float Calc_Each_Chi2(TAbsFitConstraint *constraint, float mass, float width);
  float Calc_Each_Chi2(TAbsFitParticle *ptr);
  bool Check_Repetition();
  void Clear();
  float Define_PreKin_Cut();
  bool End_Permutation() { return next_permutation(vec_permutation.begin(), vec_permutation.end()); }
  float Get_PreKin_Cut_MVA_Score(const TString &fin_path);
  void Find_Best_Permutation();
  bool Included_Matched_Jet(const int &index);
  void Index_To_Permutation();
  int Permutation_To_Index(const int &permuation);
  bool Pre_Kinematic_Cut_Chi2();
  float Pre_Kinematic_Score();
  void Print_Permutation();
  bool Quality_Cut();
  void Rebalance_Met();
  void Resol_Neutrino_Pz(){}; // not used currently
  void Save_Permutation(const bool &push = false);
  void Save_Results();
  void Set_Constraints();
  void Set_Fitter();
  void Set_Jets();
  void Set_Lepton();
  void Set_Neutrino(const int &index);
  void Set_Variables_For_MVA(const Results &results);
  void Sol_Neutrino_Pz();
  void Update_Best(const Results &results);

  ClassDef(TKinFitterDriver, 1);
};

#endif /* TKinFitterDriver_h */
