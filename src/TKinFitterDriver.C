#include "TKinFitterDriver.h"

ClassImp(TKinFitterDriver);

//////////

TKinFitterDriver::TKinFitterDriver(int _data_year, bool _rm_wm_constraint)
{
  data_year = _data_year;
  rm_wm_constraint = _rm_wm_constraint;

  fitter = new TKinFitter("Fitter", "Fitter");
  u_ptr_fitter.reset(fitter);

  fitter->setMaxNbIter(5000);
  fitter->setMaxDeltaS(1e-2);//let 1e-2 be nominal, 1e-3 be high tolerence,1e-5 be very high tol.
  fitter->setMaxF(1e-2);//let 1e-2 be nominal, 1e-3 be high tolerence, 1e-4 be very high tol.
  fitter->setVerbosity(3);
  //fitter->setVerbosity(1);
  
  reader = new TMVA::Reader("!Color:!Silent");
  reader->AddSpectator("n_jets", &n_jet);
  reader->AddVariable("had_t_b_pt", &had_t_b_pt);
  reader->AddVariable("w_u_pt", &w_u_pt);
  reader->AddVariable("w_d_pt", &w_d_pt);
  reader->AddVariable("lep_t_b_pt", &lep_t_b_pt);
  reader->AddVariable("theta_w_u_w_d", &theta_w_u_w_d);
  reader->AddVariable("theta_had_w_had_t_b", &theta_had_w_had_t_b);
  reader->AddVariable("theta_lep_neu", &theta_lep_neu);
  reader->AddVariable("theta_lep_w_lep_t_b", &theta_lep_w_lep_t_b);
  reader->AddVariable("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  reader->AddVariable("had_t_mass", &had_t_mass);
  reader->AddVariable("had_w_mass", &had_w_mass);
  reader->AddVariable("lep_t_mass", &lep_t_mass);
  reader->AddVariable("lep_t_partial_mass", &lep_t_partial_mass);
  reader->AddVariable("chi2", &chi2);
  
  
  TString weight_file_base = getenv("SKFlat_WD");
  weight_file_base += "/external/KinematicFitter/data";
  
  TString weight_file = weight_file_base + "/4Jets/TMVAClassification_BDT.weights.xml"; 
  reader->BookMVA("BDT_4Jets", weight_file);

  weight_file = weight_file_base + "/5Jets/TMVAClassification_BDT.weights.xml";
  reader->BookMVA("BDT_5Jets", weight_file);

  weight_file = weight_file_base + "/6Jets/TMVAClassification_BDT.weights.xml";
  reader->BookMVA("BDT_6Jets", weight_file);
}//TKinFitterDriver::TKinFitterDriver(int _DataYear)

//////////

void TKinFitterDriver::Scan()
{
  if(n_jet<4 || n_btag<2) return;

  results_container.Reset();
  
  //loop over jet permutation
  do
    {
      Print_Permutation();
      
      results.Reset();
      
      //each t quark should contain at least one 
      // if(!BJet_Assignment_Cut())
      //   {
      //     results.cut = CUT_RESULT::B_ASSIGNMENT;
      //     Save_Permutation(true);
      
      //     cout << "Cut B Assignment" << endl;
      
      //     continue;
      //   }
      
      //check repetition
      //For KF, permutations within W cadidates are meaningless. So reject these permutation to save time
      // if(!Check_Repetition())
      //   {
      //     results.cut = CUT_RESULT::REPETITION;
      //     Save_Permutation(true);
      
      //     continue;
      //   }
      
      Set_Jets();
      
      Rebalance_Met();
      Sol_Neutrino_Pz();
      
      //loop over neutrino permuation
      int n_neutrino_sol;
      if(chk_real_neu_pz==true) n_neutrino_sol = 2;
      else if(chk_real_neu_pz==false) n_neutrino_sol = 1;
      for(int i=0; i<n_neutrino_sol; ++i)
	{
	  index_neutrino_sol = i;
	  cout << "test neutrino pz = " << neutrino_pz_sol[index_neutrino_sol] << endl;
	  Set_Neutrino(index_neutrino_sol);

	  //very useful to reduce computation time
	  if(!Pre_Kinematic_Cut())
	    {
	      //     results.cut = CUT_RESULT::PRE_KINEMATIC;
	      //     Save_Permutation(true);
	      
	      //     cout << "Cut Pre Kin" << endl;
	      
	      //     continue;
	    }
	  
	  //contruct TFitConstraint objects
	  Set_Constraints();
	  
	  //register things to fitter
	  Set_Fitter();
	  
	  //fit
	  results.status = fitter->fit();
	  
	  //not sure quality cut is useful at this stage
	  // if(!Quality_Cut())
          //   {
          //     results.cut = CUT_RESULT::QUALITY;
	  //     Save_Permutation(true);
	  
	  //     cout << "Cut Quality cut" << endl;
	  
	  //     float chi2_jet_had_t_b = Calc_Each_Chi2(fit_jet_had_t_b);
	  //     float chi2_jet_w_u = Calc_Each_Chi2(fit_jet_w_u);
	  //     float chi2_jet_w_d = Calc_Each_Chi2(fit_jet_w_d);
	  //     float chi2_jet_lep_t_b = Calc_Each_Chi2(fit_jet_lep_t_b);

	  //     cout << "chi2_jet_had_t_b = " << chi2_jet_had_t_b << ", chi2_jet_w_u = " << chi2_jet_w_u << ", chi2_jet_w_d = " << chi2_jet_w_d << ", chi2_jet_lep_t_b = " << chi2_jet_lep_t_b << endl;
	  
          //     continue;
          //   }
	  
	  results.cut = CUT_RESULT::PASS;
	  
	  //Save
 	  Save_Results();
	}//loop over neutrino pz solution
    }while(End_Permutation());//loop over jet permutation

  
  //find best permutation
  Find_Best_Permutation();
  
  return;
}//void TKinFitterDriver::Scan()

//////////

void TKinFitterDriver::Set_Objects(vector<Jet>& _vec_jet, vector<float>& _vec_resolution_pt, vector<bool>& _vec_btag, Lepton& _lepton, Particle& _met, bool _chk_matched, int* _index_matched_jet)
{
  Clear();
  
  vec_jet = _vec_jet;
  vec_resolution_pt = _vec_resolution_pt;
  vec_btag = _vec_btag;
  lepton = _lepton;
  met = _met;
  
  chk_matched = _chk_matched;
    
  n_jet = vec_jet.size();
  int n_vec_resolution_pt = vec_resolution_pt.size();
  int n_vec_btag = vec_btag.size();
  
  if(n_jet!=n_vec_resolution_pt || n_jet!=n_vec_btag)
    {
      cerr << "The size of two vectors are not same. Check it." << endl;
      cerr << "Size of jet vector = " << n_jet << ", Size of btagger vector = " << n_vec_btag << endl;

      exit(1);
    }
  
  //in case matched jets are provided
  if(chk_matched==true)
    {
      vec_jet.clear();
      vec_resolution_pt.clear();
      vec_btag.clear();
      
      for(int i=0; i<4; ++i) index_matched_jet[i] = _index_matched_jet[i];

      int index = 0;
      for(int i=0; i<n_jet; ++i)
        {
	  if(!Included_Matched_Jet(i))
	    {
	      ++index;
	      
	      continue;
	    }
	  permutation_to_index[vec_jet.size()] = index;
	  ++index;

	  vec_jet.push_back(_vec_jet[i]);
	  vec_resolution_pt.push_back(_vec_resolution_pt[i]);
	  vec_btag.push_back(_vec_btag[i]);
	  
	  cout << "(i = " << i << ", b-tag = " << _vec_btag[i] << "), ";
	}
      cout << endl;
      
      n_jet = 4;
    
      //correct permutation
      Index_To_Permutation();
    }

  for(int i=0; i<n_jet; ++i)
    { 
      if(i==0) vec_permutation.push_back(HAD_T_B);
      else if(i==1) vec_permutation.push_back(W_U);
      else if(i==2) vec_permutation.push_back(W_D);
      else if(i==3) vec_permutation.push_back(LEP_T_B);
      else  vec_permutation.push_back(OTHERS);
  
      if(vec_btag.at(i)==true) n_btag++;
    }
  
  Set_Lepton();
  
  return;
}//void TKinFitterDriver::Set_Objects(vector<Jet>& _vec_jet, vector<bool>& _vec_btag, Lepton& _lepton, Particle& _met)

//////////

bool TKinFitterDriver::BJet_Assignment_Cut()
{
  //current permutation
  int index_had_t_b=-1;
  int index_lep_t_b=-1;

  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);

      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B) index_had_t_b = i;
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B) index_lep_t_b = i;
    }

  if(vec_btag.at(index_had_t_b)==false || vec_btag.at(index_lep_t_b)==false) return false;

  return true;
}//bool BJet_Assignment_Cut()

//////////

float TKinFitterDriver::Calc_Chi2()
{
  float chi2 = 0;
  chi2 += Calc_Each_Chi2(fit_jet_had_t_b);
  chi2 += Calc_Each_Chi2(fit_jet_w_u);
  chi2 += Calc_Each_Chi2(fit_jet_w_d);
  chi2 += Calc_Each_Chi2(fit_jet_lep_t_b);

  for(unsigned int i=0; i<vec_fit_jet_extra.size(); i++){ chi2 += Calc_Each_Chi2(vec_fit_jet_extra.at(i)); }

  chi2 += Calc_Each_Chi2(fit_lepton);
  
  chi2 += Calc_Each_Chi2(constraint_had_t_mgaus, T_MASS, T_WIDTH);
  chi2 += Calc_Each_Chi2(constraint_lep_t_mgaus, T_MASS, T_WIDTH);
  chi2 += Calc_Each_Chi2(constraint_had_w_mgaus, W_MASS, W_WIDTH);
  chi2 += Calc_Each_Chi2(constraint_lep_w_mgaus, W_MASS, W_WIDTH);

  // cout << "test calc_chi2" << endl;
  // cout << "chi2_had_t_b = " << Calc_Each_Chi2(fit_jet_had_t_b) << endl;
  // cout << "chi2_w_u = " << Calc_Each_Chi2(fit_jet_w_u) << endl;
  // cout << "chi2_w_d = " << Calc_Each_Chi2(fit_jet_w_d) << endl;
  // cout << "chi2_lep_t_b = " << Calc_Each_Chi2(fit_jet_lep_t_b) << endl;
  // cout << "chi2_lepton = " << Calc_Each_Chi2(fit_lepton) << endl;
  // cout << "chi2_had_t = " << Calc_Each_Chi2(constraint_had_t_mgaus, T_MASS, T_WIDTH) << endl;
  // cout << "chi2_lep_t = " << Calc_Each_Chi2(constraint_lep_t_mgaus, T_MASS, T_WIDTH) << endl;
  // cout << "chi2_had_w = " << Calc_Each_Chi2(constraint_had_w_mgaus, W_MASS, W_WIDTH) << endl;
  // cout << "chi2_lep_w = " << Calc_Each_Chi2(constraint_lep_w_mgaus, W_MASS, W_WIDTH) << endl;
    
  return chi2;
}//float TKinFitterDriver::Calc_Chi2()

//////////

float TKinFitterDriver::Calc_Each_Chi2(TAbsFitConstraint* constraint, float mass, float width)
{
  const TMatrixD* currPar = constraint->getParCurr();
  float deltaY = 1 - (*currPar)(0, 0);
  float chi2 = (mass*mass)*(deltaY*deltaY)/(width*width); 

  return chi2;
}//float Calc_Each_Chi2(TAbsFitConstraint* ptr, float mass, float width)

//////////

float TKinFitterDriver::Calc_Each_Chi2(TAbsFitParticle* ptr)
{
  const TMatrixD* iniPar = ptr->getParIni();
  const TMatrixD* currPar = ptr->getParCurr();
  const TMatrixD* covMatrix = ptr->getCovMatrix();

  float chi2 = 0;
  for(Int_t i=0; i<iniPar->GetNcols(); i++)
    {
      float deltaY = (*iniPar)(i,i) - (*currPar)(i,i);
      chi2 += deltaY*deltaY/(*covMatrix)(i,i);
    }
  
  return chi2;
}//float TKinFitterDriver::Calc_Each_Chi2(TAbsFitParticle* ptr)
 
//////////

bool TKinFitterDriver::Check_Repetition()
{
  //Although vec_permutation distinguishes w_u, and w_d, fitter doesn't care this permutation. Therefore this permutation would be ignored to save computation time 
  //the above is obsolete because BJetRegression do distinguish

  //current permutation
  int index_had_t_b=-1;
  int index_w_u=-1;
  int index_w_d=-1;
  int index_lep_t_b=-1;
  
  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);
      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B) index_had_t_b = i;
      else if(jet_assignment==JET_ASSIGNMENT::W_U) index_w_u = i;
      else if(jet_assignment==JET_ASSIGNMENT::W_D) index_w_d = i;
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B) index_lep_t_b = i;
    }
  
  //check repetition 
  for(auto& results : results_container.vec_results)
    {
      if(results.index_neutrino_sol==index_neutrino_sol&&
	 results.index_had_t_b==index_had_t_b && 
	 results.index_lep_t_b==index_lep_t_b &&
	 results.index_w_u==index_w_d &&
	 results.index_w_d==index_w_u) return false;//repetition found
    }

  //no repetition
  return true;
}//bool TKinFitterDriver::Check_Repetition()

//////////

 void TKinFitterDriver::Clear()
 {
   vec_jet.clear();
   vec_jet.shrink_to_fit();

   vec_resolution_pt.clear();
   vec_resolution_pt.shrink_to_fit();
   
   vec_btag.clear();
   vec_btag.shrink_to_fit();
   
   vec_permutation.clear();
   vec_permutation.shrink_to_fit();
   
   chk_matched = false;
   
   for(int i=0; i<4; ++i)
     {
       index_matched_jet[i] = -999;
     }

   n_jet = 0;
   n_btag = 0;
   
   return;
   

 }//void TKinFitterDriver::Clear()

//////////

void TKinFitterDriver::Find_Best_Permutation()
{
  bool chk_mva = true;
  
  int i_best = -1;
  
  //using MVA
  if(chk_mva)
    {
      float mva_score_best = -999; 
      
      for(unsigned int i=0; i<results_container.vec_results.size(); ++i)
	{
	  Results results = results_container.vec_results.at(i);
	  
	  if(results.cut!=CUT_RESULT::PASS) continue;
	  
	  had_t_b_pt = vec_jet[results.index_had_t_b].Pt();
	  w_u_pt = vec_jet[results.index_w_u].Pt();
	  w_d_pt = vec_jet[results.index_w_d].Pt();
	  lep_t_b_pt = vec_jet[results.index_lep_t_b].Pt();
	  
	  theta_w_u_w_d = results.theta_w_u_w_d;
	  theta_had_w_had_t_b = results.theta_had_w_had_t_b;
	  theta_lep_neu = results.theta_lep_neu;
	  theta_lep_w_lep_t_b = results.theta_lep_w_lep_t_b;
	  del_phi_had_t_lep_t = results.del_phi_had_t_lep_t;
	  
	  had_t_mass = results.initial_had_t_mass;
	  had_w_mass = results.initial_had_w_mass;
	  lep_t_mass = results.initial_lep_t_mass;
	  lep_t_partial_mass = results.initial_lep_t_partial_mass;

	  chi2 = results.chi2;

	  cout << had_t_b_pt << " " << w_u_pt << " " << w_d_pt << " " << lep_t_b_pt << " " << theta_w_u_w_d << " " << theta_had_w_had_t_b << " " << theta_lep_neu << " " << theta_lep_w_lep_t_b << " " << del_phi_had_t_lep_t << " " << had_t_mass << " " << had_w_mass << " " << lep_t_mass << " " << lep_t_partial_mass << " " << chi2 << endl; 
	  float mva_score;
	  if(n_jet==4) mva_score = reader->EvaluateMVA("BDT_4Jets");
	  else if(n_jet==5) mva_score = reader->EvaluateMVA("BDT_5Jets");
	  else mva_score = reader->EvaluateMVA("BDT_6Jets");
	  
	  results.mva_score = mva_score;

	  if(mva_score_best<mva_score)
	    {
	      mva_score_best = mva_score; 
	      i_best = i;
	    }

	  cout << "test find_best_permutation " << i << " " << results.index_had_t_b << " " << results.index_w_u << " " << results.index_w_d << " " << results.index_lep_t_b << ") " << mva_score << " " << mva_score_best << " " << i_best << endl;
	}
    }//using MVA
  
  //using simple chi2 minimization
  else
    {
      float chi2_best = 99999;
      
      for(unsigned int i=0; i<results_container.vec_results.size(); i++)
	{
	  Results results = results_container.vec_results.at(i);
           
	  if(results.chi2<chi2_best && results.cut==CUT_RESULT::PASS)
	    {
	      i_best = i;
	      chi2_best = results.chi2;
	    }
	}
    }//using simple chi2 minimization
      
      
  if(i_best!=-1)
    {
      results_container.status = true;
            
      Results results = results_container.vec_results.at(i_best);
      
      results_container.best_chi2 = results.chi2;
      
      results_container.best_index_had_t_b = results.index_had_t_b;
      results_container.best_index_w_u = results.index_w_u;
      results_container.best_index_w_d = results.index_w_d;
      results_container.best_index_lep_t_b = results.index_lep_t_b;
      
      results_container.best_initial_had_t_mass = results.initial_had_t_mass;
      results_container.best_initial_had_w_mass = results.initial_had_w_mass;
      results_container.best_initial_lep_t_mass = results.initial_lep_t_mass;
      results_container.best_initial_lep_w_mass = results.initial_lep_w_mass;
      
      results_container.best_fitted_had_t_mass = results.fitted_had_t_mass;
      results_container.best_fitted_had_w_mass = results.fitted_had_w_mass;
      results_container.best_fitted_lep_t_mass = results.fitted_lep_t_mass;
      results_container.best_fitted_lep_w_mass = results.fitted_lep_w_mass;
    }
  
  cout << "Best Permutation = " << results_container.best_index_had_t_b << " " << results_container.best_index_w_u << " " << results_container.best_index_w_d << " " << results_container.best_index_lep_t_b << endl;
  //cout << results_container.best_chi2 << " " <<  results_container.best_fitted_had_t_mass << " " <<  results_container.best_fitted_had_w_mass << endl;
  
  return;
}//void TKinFitterDriver::Find_Best_Permutation()

//////////

bool TKinFitterDriver::Included_Matched_Jet(const int& index)
{
  for(int i=0; i<4; ++i)
    {
      if(index==index_matched_jet[i]) return true;
    }
  
  return false;
}//bool TKinFitterDriver::Included_Matched_Jet(const int& i)

//////////

void TKinFitterDriver::Index_To_Permutation()
{
  cout << "test Index_To_Permutation" << endl;

  for(int i=0; i<4; ++i) correct_permutation[i] = -999;

  for(int i=0; i<4; ++i)
    {
      int index = 0;
      int rank = 0;
      int smallest = 999;
      for(int j=0; j<4; ++j)
        {
          if(correct_permutation[j]!=-999)
            {
              ++rank;
              continue;
            }

          if(index_matched_jet[j]<smallest)
            {
              smallest = index_matched_jet[j];
              index = j;
            }
        }

      correct_permutation[index] = rank;
    }

  for(int i=0; i<4; ++i)
    {
      cout << correct_permutation[i] << " ";
    }
  cout << endl;
  
  return;
}//void TKinFitterDriver::Index_To_Permutation()

//////////

int Permutation_To_Index(const int& permutation)
{ 
  cout << "Permutation_To_Index" << endl;
  //vec_permutation
  int index = 0;

  return index;
}//void Permutation_To_Index()

//////////

bool TKinFitterDriver::Pre_Kinematic_Cut()
{
  //convention
  //true: event passes the cut
  //false: event is removed by the cut
  
  //Very useful to reduce computation time
  //Not for increase SNR 
  //Cut 5% of signal for each 
  
  //jet pt cut
  
  //tt delta phi
  TLorentzVector had_w = jet_w_u + jet_w_d;
  TLorentzVector had_t = jet_had_t_b + jet_w_u + jet_w_d;
  TLorentzVector lep_t = jet_lep_t_b + lepton + neutrino;
  
  float del_phi_had_t_lep_t = TMath::Abs(had_t.DeltaPhi(lep_t));
  if(n_jet==4 && del_phi_had_t_lep_t<1.9) return false; 
  if(n_jet==5 && del_phi_had_t_lep_t<0.82) return false;
  else if(6<=n_jet && del_phi_had_t_lep_t<0.42) return false;

  //had_t_mass
  float had_t_mass = had_t.M();
  if(n_jet==4 && (had_t_mass<132.5 || 270.5<had_t_mass)) return false;
  else if(n_jet==5 && (had_t_mass<129.5 || 280.5<had_t_mass)) return false;
  else if(6<=n_jet && (had_t_mass<127.5 || 321.5<had_t_mass)) return false;
  
  //lep_t_mass
  float lep_t_mass = lep_t.M();
  if(n_jet==4 && (lep_t_mass<138.5 || 329.5<lep_t_mass)) return false;
  else if(n_jet==5 && (lep_t_mass<135.5 || 345.5<lep_t_mass)) return false;
  else if(6<=n_jet && (lep_t_mass<134.5 || 363.5<lep_t_mass)) return false;
	     
  //had_w_mass
  float had_w_mass = had_w.M();
  if(n_jet==4 && (had_w_mass<51.5 || 153.5<had_w_mass)) return false;
  else if(n_jet==5 && (had_w_mass<50.5 || 162.5<had_w_mass)) return false;
  else if(n_jet<=6 && (had_w_mass<50.5 || 194.5<had_w_mass)) return false;
   
  return true;
}//bool TKinFitterDriver::Pre_Kinematic_Cut()

//////////

void TKinFitterDriver::Print_Permutation()
{
  //current permutation
  int index_had_t_b=-1;
  int index_w_u=-1;
  int index_w_d=-1;
  int index_lep_t_b=-1;

  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);
      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B) index_had_t_b = i;
      else if(jet_assignment==JET_ASSIGNMENT::W_U) index_w_u = i;
      else if(jet_assignment==JET_ASSIGNMENT::W_D) index_w_d = i;
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B) index_lep_t_b = i;
    }

  cout << "Index (had_t_b, w_u, w_d, lep_t_b) =  " << index_had_t_b << " "  << index_w_u << " " << index_w_d << " " << index_lep_t_b << endl;

  return;
}//void TKinFitterDriver::Print_Permutation()

//////////

bool TKinFitterDriver::Quality_Cut()
{
  float chi2_jet_had_t_b = Calc_Each_Chi2(fit_jet_had_t_b);
  float chi2_jet_w_u = Calc_Each_Chi2(fit_jet_w_u);
  float chi2_jet_w_d = Calc_Each_Chi2(fit_jet_w_d);
  float chi2_jet_lep_t_b = Calc_Each_Chi2(fit_jet_lep_t_b);

  //results.status=0 converged, results.status=1 not converged
  if(results.status==0 && chi2_jet_had_t_b<25 && chi2_jet_w_u<25 &&  chi2_jet_w_d<25 && chi2_jet_lep_t_b<25) return true;
  else return false;

  //dummy
  return true;
}//bool TKinFitterDriver::Quality_Cut()

//////////

void TKinFitterDriver::Rebalance_Met()
{
  float met_rebalance_px = 0;
  float met_rebalance_py = 0;

  for(unsigned int i=0; i<vec_jet.size(); ++i)
    {
      met_rebalance_px += vec_jet[i].Px();
      met_rebalance_py += vec_jet[i].Py();
    }

  //rebalacne met
  met_rebalance_px -= jet_had_t_b.Px();
  met_rebalance_py -= jet_had_t_b.Py();
  
  met_rebalance_px -= jet_w_u.Px();
  met_rebalance_py -= jet_w_u.Py();

  met_rebalance_px -= jet_w_d.Px();
  met_rebalance_py -= jet_w_d.Py();

  met_rebalance_px -= jet_lep_t_b.Px();
  met_rebalance_py -= jet_lep_t_b.Py();

  for(unsigned int i=0; i<vec_jet_extra.size(); ++i)
    {
      met_rebalance_px -= vec_jet_extra[i].Px();
      met_rebalance_py -= vec_jet_extra[i].Py();
    }

  met_rebalance_px = met.Px() + met_rebalance_px;
  met_rebalance_py = met.Py() + met_rebalance_py;
  
  met_rebalance.SetPxPyPzE(met_rebalance_px, met_rebalance_py, 0, TMath::Sqrt(met_rebalance_px*met_rebalance_px+met_rebalance_py*met_rebalance_py));
  
  cout << "test met_rebal " << met_rebalance_px << " " << met_rebalance_py << endl;

  return;
}//void TKineFitterDriver::Rebalance_Met()

//////////

void TKinFitterDriver::Save_Permutation(const bool& push)
{
  results.vec_permutation = vec_permutation;
  
  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);
  
      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B)
	{
	  if(chk_matched==true) results.index_had_t_b = permutation_to_index[i];
	  else results.index_had_t_b = i;
	}
      else if(jet_assignment==JET_ASSIGNMENT::W_U)
	{
	  if(chk_matched==true) results.index_w_u = permutation_to_index[i];
	  else results.index_w_u = i;
	}
      else if(jet_assignment==JET_ASSIGNMENT::W_D)
	{
	  if(chk_matched==true) results.index_w_d = permutation_to_index[i];
	  else results.index_w_d = i;
	}
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B)
	{
	  if(chk_matched==true) results.index_lep_t_b = permutation_to_index[i];
	  else results.index_lep_t_b = i;
	}
    }

  if(push) results_container.vec_results.push_back(results);

  results.index_neutrino_sol = index_neutrino_sol;

  return;
}//void TKinFitterDriver::Save_Permutation()

//////////

void TKinFitterDriver::Save_Results()
{
  Save_Permutation();
  
  results.chk_real_neu_pz = chk_real_neu_pz;

  //initial objects
  results.met_px = met.Px();
  results.met_py = met.Py();
  
  results.met_rebalance_px = neutrino_px_init;
  results.met_rebalance_py = neutrino_py_init;
  
  results.neutrino_pz_sol = neutrino_pz_init;
  
  results.chk_real_neu_pz = chk_real_neu_pz;

  TLorentzVector initial_had_t = jet_had_t_b + jet_w_u + jet_w_d;
  TLorentzVector initial_had_w = jet_w_u + jet_w_d;
  TLorentzVector initial_lep_t = jet_lep_t_b + lepton + neutrino;
  TLorentzVector initial_lep_w = lepton + neutrino;

  results.del_phi_w_u_w_d = TMath::Abs(jet_w_u.DeltaPhi(jet_w_d));
  results.del_phi_had_w_had_t_b = TMath::Abs(initial_had_w.DeltaPhi(jet_had_t_b));
  results.del_phi_lep_neu = TMath::Abs(lepton.DeltaPhi(neutrino));
  results.del_phi_lep_w_lep_t_b = TMath::Abs(initial_lep_w.DeltaPhi(jet_lep_t_b));
  results.del_phi_had_t_lep_t = TMath::Abs(initial_had_t.DeltaPhi(initial_lep_t));

  results.del_eta_w_u_w_d = TMath::Abs(jet_w_u.Eta()-jet_w_d.Eta());
  results.del_eta_had_w_had_t_b = TMath::Abs(initial_had_w.Eta()-jet_had_t_b.Eta());
  results.del_eta_lep_neu = TMath::Abs(lepton.Eta()-neutrino.Eta());
  results.del_eta_lep_w_lep_t_b = TMath::Abs(initial_lep_w.Eta()-jet_lep_t_b.Eta());
  results.del_eta_had_t_lep_t = TMath::Abs(initial_had_t.Eta()-initial_lep_t.Eta());

  results.del_r_w_u_w_d = jet_w_u.DeltaR(jet_w_d);
  results.del_r_had_w_had_t_b = initial_had_w.DeltaR(jet_had_t_b);
  results.del_r_lep_neu = lepton.DeltaR(neutrino);
  results.del_r_lep_w_lep_t_b = initial_lep_w.DeltaR(jet_lep_t_b);
  results.del_r_had_t_lep_t = initial_had_t.DeltaR(initial_lep_t);
  
  results.theta_w_u_w_d = jet_w_u.Angle(jet_w_d.Vect());
  results.theta_had_w_had_t_b = initial_had_w.Angle(jet_had_t_b.Vect());
  results.theta_lep_neu = lepton.Angle(neutrino.Vect());
  results.theta_lep_w_lep_t_b = initial_lep_w.Angle(jet_lep_t_b.Vect());
  results.theta_had_t_lep_t = initial_had_t.Angle(initial_lep_t.Vect());
  
  //construct fitted objects
  const TLorentzVector* tmp = fit_jet_had_t_b->getCurr4Vec();
  TLorentzVector fitted_jet_had_t_b(tmp->Px(), tmp->Py(), tmp->Pz(), tmp->E());

  tmp = fit_jet_w_u->getCurr4Vec();
  TLorentzVector fitted_jet_w_u(tmp->Px(), tmp->Py(), tmp->Pz(), tmp->E()); 

  tmp = fit_jet_w_d->getCurr4Vec();
  TLorentzVector fitted_jet_w_d(tmp->Px(), tmp->Py(), tmp->Pz(), tmp->E());

  tmp = fit_jet_lep_t_b->getCurr4Vec();
  TLorentzVector fitted_jet_lep_t_b(tmp->Px(), tmp->Py(), tmp->Pz(), tmp->E());
  
  tmp = fit_lepton->getCurr4Vec();
  TLorentzVector fitted_lepton(tmp->Px(), tmp->Py(), tmp->Pz(), tmp->E());
  
  tmp = fit_neutrino->getCurr4Vec();
  TLorentzVector fitted_neutrino(tmp->Px(), tmp->Py(), tmp->Pz(), tmp->E());
  
  TLorentzVector fitted_had_t = fitted_jet_had_t_b + fitted_jet_w_u + fitted_jet_w_d;
  TLorentzVector fitted_had_w = fitted_jet_w_u + fitted_jet_w_d;
  TLorentzVector fitted_lep_t = fitted_jet_lep_t_b + fitted_lepton + fitted_neutrino;
  TLorentzVector fitted_lep_w = fitted_lepton + fitted_neutrino;
  
  results.chi2_jet_had_t_b = Calc_Each_Chi2(fit_jet_had_t_b);
  results.chi2_jet_w_u = Calc_Each_Chi2(fit_jet_w_u);
  results.chi2_jet_w_d = Calc_Each_Chi2(fit_jet_w_d);;
  results.chi2_jet_lep_t_b = Calc_Each_Chi2(fit_jet_lep_t_b);
  
  results.chi2_jet_extra = 0;
  for(unsigned int i=0; i<vec_fit_jet_extra.size(); i++){ results.chi2_jet_extra += Calc_Each_Chi2(vec_fit_jet_extra.at(i)); }
  
  results.chi2_constraint_had_t = Calc_Each_Chi2(constraint_had_t_mgaus, T_MASS, T_WIDTH);
  results.chi2_constraint_had_w = Calc_Each_Chi2(constraint_had_w_mgaus, W_MASS, W_WIDTH);
  results.chi2_constraint_lep_t = Calc_Each_Chi2(constraint_lep_t_mgaus, T_MASS, T_WIDTH);
  results.chi2_constraint_lep_w = Calc_Each_Chi2(constraint_lep_w_mgaus, W_MASS, W_WIDTH);

  results.chi2 = Calc_Chi2();
  
  results.initial_had_t_mass = initial_had_t.M();
  results.initial_had_w_mass = initial_had_w.M();
  results.initial_lep_t_mass = initial_lep_t.M();
  results.initial_lep_t_partial_mass = (jet_lep_t_b + lepton).M();
  results.initial_lep_w_mass = initial_lep_w.M();

  results.fitted_had_t_mass = fitted_had_t.M();
  results.fitted_had_w_mass = fitted_had_w.M();
  results.fitted_lep_t_mass = fitted_lep_t.M();
  results.fitted_lep_w_mass = fitted_lep_w.M();

  results_container.vec_results.push_back(results);
  
  cout << "chi2 jet_lep_t_b = " << results.chi2_jet_lep_t_b << ", chi2_jet_extra = " << results.chi2_jet_extra <<  endl; 

  return;
}//void TKinFitterDriver::Save_Results()

//////////

void TKinFitterDriver::Set_Constraints()
{
  //hadronic top mass constraint
  constraint_had_t_mgaus = new TFitConstraintMGaus("Hadronic_Top_Mass", "Hadronic_Top_Mass", NULL, NULL, T_MASS, T_WIDTH);
  constraint_had_t_mgaus->addParticle1(fit_jet_had_t_b);
  constraint_had_t_mgaus->addParticle1(fit_jet_w_u);
  constraint_had_t_mgaus->addParticle1(fit_jet_w_d);
  
  u_ptr_constraint_had_t_mgaus.reset(constraint_had_t_mgaus);
  
  //leptonic top mass constraint
  constraint_lep_t_mgaus = new TFitConstraintMGaus("Leptonic_Top_Mass", "Leptonic_Top_Mass", NULL, NULL, T_MASS, T_WIDTH);
  constraint_lep_t_mgaus->addParticle1(fit_jet_lep_t_b);
  constraint_lep_t_mgaus->addParticle1(fit_lepton);
  constraint_lep_t_mgaus->addParticle1(fit_neutrino);
  
  u_ptr_constraint_lep_t_mgaus.reset(constraint_lep_t_mgaus);

  //hadronic W mass constraint
  constraint_had_w_mgaus = new TFitConstraintMGaus("Hadronic_W_Mass", "Hadronic_W_Mass", NULL, NULL, W_MASS, W_WIDTH);
  constraint_had_w_mgaus->addParticle1(fit_jet_w_u);
  constraint_had_w_mgaus->addParticle1(fit_jet_w_d);
  
  u_ptr_constraint_had_w_mgaus.reset(constraint_had_w_mgaus);

  //leptonic W mass constraint
  constraint_lep_w_mgaus = new TFitConstraintMGaus("Leptonic_W_Mass", "Leptonic_W_Mass", NULL, NULL, W_MASS, W_WIDTH);
  constraint_lep_w_mgaus->addParticle1(fit_lepton);
  constraint_lep_w_mgaus->addParticle1(fit_neutrino);

  u_ptr_constraint_lep_w_mgaus.reset(constraint_lep_w_mgaus);

  //px and py balance constraint
  constraint_px = new TFitConstraintEp("Px", "Px", TFitConstraintEp::component::pX, 0.);
  constraint_py = new TFitConstraintEp("Py", "Py", TFitConstraintEp::component::pY, 0.);

  constraint_px->addParticle(fit_lepton);
  constraint_py->addParticle(fit_lepton);

  constraint_px->addParticle(fit_neutrino);
  constraint_py->addParticle(fit_neutrino);

  constraint_px->addParticle(fit_jet_had_t_b);
  constraint_py->addParticle(fit_jet_had_t_b);
  
  constraint_px->addParticle(fit_jet_w_u);
  constraint_py->addParticle(fit_jet_w_u);

  constraint_px->addParticle(fit_jet_w_d);
  constraint_py->addParticle(fit_jet_w_d);

  constraint_px->addParticle(fit_jet_lep_t_b);
  constraint_py->addParticle(fit_jet_lep_t_b);

  for(auto& fit_jet : vec_fit_jet_extra)
    {
      constraint_px->addParticle(fit_jet);
      constraint_py->addParticle(fit_jet);
    }

  float px_initial = constraint_px->getInitValue();
  constraint_px->setConstraint(px_initial);

  float py_initial = constraint_py->getInitValue();
  constraint_py->setConstraint(py_initial);
  
  u_ptr_constraint_px.reset(constraint_px);
  u_ptr_constraint_py.reset(constraint_py);
  
  return;
}//void TKinFitterDriver::Set_Constraints()

//////////

void TKinFitterDriver::Set_Fitter()
{
  fitter->reset();

  //register measured particles
  fitter->addMeasParticle(fit_lepton);
  fitter->addMeasParticle(fit_jet_had_t_b);
  fitter->addMeasParticle(fit_jet_w_u);
  fitter->addMeasParticle(fit_jet_w_d);
  fitter->addMeasParticle(fit_jet_lep_t_b);
  for(auto& fit_jet_extra : vec_fit_jet_extra) fitter->addMeasParticle(fit_jet_extra);
  
  //register unmeasure particles
  fitter->addUnmeasParticle(fit_neutrino);

  //register constraints
  fitter->addConstraint(constraint_had_t_mgaus);
  fitter->addConstraint(constraint_lep_t_mgaus);
  if(!rm_wm_constraint) fitter->addConstraint(constraint_had_w_mgaus);
  fitter->addConstraint(constraint_lep_w_mgaus);
  fitter->addConstraint(constraint_px);
  fitter->addConstraint(constraint_py);

  return;
}//void TKinFitterDriver::Set_Fitter()

//////////

void TKinFitterDriver::Set_Jets()
{
  //clear
  error_jet_had_t_b.ResizeTo(1, 1);
  error_jet_w_u.ResizeTo(1, 1);
  error_jet_w_d.ResizeTo(1, 1);
  error_jet_lep_t_b.ResizeTo(1, 1);

  vec_jet_extra.clear();
  vec_jet_extra.shrink_to_fit();
  
  vec_u_fit_jet_extra.clear();
  vec_u_fit_jet_extra.shrink_to_fit();
  
  vec_fit_jet_extra.clear();

  vec_error_jet_extra.clear();
  vec_error_jet_extra.shrink_to_fit();
 
  //assign
  vector<float> vec_jer_extra;
  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);
      
      Jet jet = vec_jet.at(i);
      
      float pt = jet.Pt();
      float eta = jet.Eta();
      float phi = jet.Phi();
      float m = jet.M();

      //for easy handling
      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B)
	{
	  results.index_had_t_b = i;
	  
	  jet_had_t_b = jet;
	  
	  float corr = jet.GetBJetRegressionNN("BCorr"); 
	  jet_had_t_b.SetPtEtaPhiM(corr*pt, eta, phi, m);
	  
	  float jer = jet.GetBJetRegressionNN("BRes"); 
	  error_jet_had_t_b(0, 0) = jer*jer*pt*pt;
	} 
      else if(jet_assignment==JET_ASSIGNMENT::W_U)
	{
	  results.index_w_u = i;
	  
	  jet_w_u = jet;
	  
	  float corr = jet.GetBJetRegressionNN("CCorr");
	  jet_w_u.SetPtEtaPhiM(corr*pt, eta, phi, m);

	  float jer = jet.GetBJetRegressionNN("CRes");
	  error_jet_w_u(0, 0) = jer*jer*pt*pt;
	}
      else if(jet_assignment==JET_ASSIGNMENT::W_D) 
	{
	  results.index_w_d = i;
	  
	  jet_w_d = jet;

	  float corr = jet.GetBJetRegressionNN("BCorr");
          jet_w_d.SetPtEtaPhiM(corr*pt, eta, phi, m);
	  
	  float jer = jet.GetBJetRegressionNN("BRes");
	  error_jet_w_d(0, 0) = jer*jer*pt*pt;
	}
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B) 
	{
	  results.index_lep_t_b = i;
	  
	  jet_lep_t_b = jet;

	  float corr = jet.GetBJetRegressionNN("BCorr");
          jet_lep_t_b.SetPtEtaPhiM(corr*pt, eta, phi, m);
	  
	  float jer = jet.GetBJetRegressionNN("BRes");
	  error_jet_lep_t_b(0, 0) = jer*jer*pt*pt;
	}
      else if(jet_assignment==JET_ASSIGNMENT::OTHERS)
	{
	  vec_jet_extra.push_back(jet);
	  
	  float jer = 0;
	  //b tagged
	  if(vec_btag.at(i)==true) jer = jet.GetBJetRegressionNN("BRes");
	  else jer = vec_resolution_pt.at(i);
	  
	  vec_jer_extra.push_back(jer);
	}
    }
 
  //construct fitting objects
  //too verbose code. but easier to read
  //b jet from hadronic top
  fit_jet_had_t_b = new TFitParticlePt("B_From_Hadronic_Top", "B_From_Hadronic_Top", &jet_had_t_b, &error_jet_had_t_b);
  u_fit_jet_had_t_b.reset(fit_jet_had_t_b);

  //u type jet from W
  fit_jet_w_u = new TFitParticlePt("U_Type_From_W_From_Hadronic_Top", "U_Type_From_W_From_Hadronic_Top", &jet_w_u, &error_jet_w_u);
  u_fit_jet_w_u.reset(fit_jet_w_u);

  //d type jet from W
  fit_jet_w_d = new TFitParticlePt("D_Type_From_W_From_Hadronic_Top", "D_Type_From_W_From_Hadronic_Top", &jet_w_d, &error_jet_w_d);
  u_fit_jet_w_d.reset(fit_jet_w_d);

  //b jet from leptonic top  
  fit_jet_lep_t_b = new TFitParticlePt("B_From_Leptonic_Top", "B_From_Leptonic_Top", &jet_lep_t_b, &error_jet_lep_t_b);
  u_fit_jet_lep_t_b.reset(fit_jet_lep_t_b);
 
  //jet extra
  for(unsigned int i=0; i<vec_jet_extra.size(); i++)
    {
      TLorentzVector jet_extra = vec_jet_extra.at(i);
      float jer_extra = vec_jer_extra.at(i);

      TMatrixD error_jet_extra;
      error_jet_extra.ResizeTo(1, 1);
      
      float pt = jet_extra.Pt();
      error_jet_extra(0, 0) = jer_extra*jer_extra*pt*pt;
      vec_error_jet_extra.push_back(error_jet_extra);
      
      TString name_jet_extra = "Jet_Extra" + TString::Itoa(i, 10);
      TFitParticlePt* fit_jet_extra = new TFitParticlePt(name_jet_extra, name_jet_extra, &vec_jet_extra.at(i), &vec_error_jet_extra.at(i)); 
      vec_fit_jet_extra.push_back(fit_jet_extra);

      unique_ptr<TFitParticlePt> u_fit_jet_extra;
      u_fit_jet_extra.reset(fit_jet_extra);
      vec_u_fit_jet_extra.push_back(move(u_fit_jet_extra));
    }
  
  return;
}//void TKinFitterDriver::Set_Jets()

//////////

void TKinFitterDriver::Set_Lepton()
{
  //construct fitting object
  float pt = lepton.Pt();
  
  error_lepton.ResizeTo(1, 1);
  error_lepton(0, 0) = TMath::Power(0.0001*pt, 2);
  
  fit_lepton = new TFitParticlePt("Fit_Lepton", "Fit_Lepton", &lepton, &error_lepton);
  u_fit_lepton.reset(fit_lepton);
  
  return;
}//void TKinFitterDriver::Set_Lepton()

//////////

void TKinFitterDriver::Set_Neutrino(const int& index)
{
  //float px = met.Px();
  //float py = met.Py();
  neutrino_px_init = met_rebalance.Px();
  neutrino_py_init = met_rebalance.Py();
  neutrino_pz_init = neutrino_pz_sol[index];
  
  neutrino.SetPxPyPzE(neutrino_px_init, neutrino_py_init, neutrino_pz_init, TMath::Sqrt(neutrino_px_init*neutrino_px_init+neutrino_py_init*neutrino_py_init+neutrino_pz_init*neutrino_pz_init));
  
  //construct fitting object
  neutrino_vec3 = neutrino.Vect();
  
  fit_neutrino = new TFitParticleMCCart("Fit_Neutrino", "Fit_Neutrino", &neutrino_vec3, 0., NULL);
  u_fit_neutrino.reset(fit_neutrino);
  
  return;
}//void TKinFitterDriver::Set_Neurino(const int& index)
 
//////////
 
void TKinFitterDriver::Sol_Neutrino_Pz()
{
  double lepton_mass =  lepton.M();
  
  float met_px = met_rebalance.Px();
  float met_py = met_rebalance.Py();
  float met_pt = met_rebalance.Pt();

  double k = TMath::Power(W_MASS, 2.)/2.0 - lepton_mass*lepton_mass/2.0 + lepton.Px()*met_px + lepton.Py()*met_py;
  double a = TMath::Power(lepton.Px(), 2.0) + TMath::Power(lepton.Py(), 2.0);
  double b = -2*k*lepton.Pz();
  double c = TMath::Power(lepton.Pt(), 2.0)*TMath::Power(met_pt, 2.0) - TMath::Power(k, 2.0);

  double determinant = TMath::Power(b, 2.0) - 4*a*c; 

  //real solution
  if(determinant>=0)
    {
      neutrino_pz_sol[0] = (-b + TMath::Sqrt(determinant))/(2*a);
      neutrino_pz_sol[1] = (-b - TMath::Sqrt(determinant))/(2*a);
      
      chk_real_neu_pz = true;
      cout << "Neu Real sol = " << neutrino_pz_sol[0] << " " << neutrino_pz_sol[1] << endl;
    }
  //complex solution
  else
    {
      neutrino_pz_sol[0] = -b/(2*a);
      //Resol_Neutrino_Pt();
      chk_real_neu_pz = false;
    
      cout << "Neu Im sol = " << neutrino_pz_sol[0] << endl;
    }
  
  return;
}//void TKinFitterDriver::Sol_Neutrino_Pz()
 
//////////
