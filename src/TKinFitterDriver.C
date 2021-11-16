#include "TKinFitterDriver.h"

ClassImp(TKinFitterDriver);

//////////

TKinFitterDriver::TKinFitterDriver(int _data_year)
{
  data_year = _data_year;

  fitter = new TKinFitter("Fitter", "Fitter");
  u_ptr_fitter.reset(fitter);

  fitter->setMaxNbIter(5000);
  fitter->setMaxDeltaS(1e-2);//let 1e-2 be nominal, 1e-3 be high tolerence,1e-5 be very high tol.
  fitter->setMaxF(1e-2);//let 1e-2 be nominal, 1e-3 be high tolerence, 1e-4 be very high tol.
  fitter->setVerbosity(3);
  //fitter->setVerbosity(1);
  
  
}//TKinFitterDriver::TKinFitterDriver(int _DataYear)

//////////

void TKinFitterDriver::Scan()
{
  if(n_jet<4 || n_btag<2) return;

  results_container.Reset();
  
  Sol_Neutrino_Pz();
  
  //loop over neutrino pz solution
  index_neutrino_sol = 0;
  while(true)
    {
      if(index_neutrino_sol==2) break;
      
      //if neutrino pz solution gives complex, take only real part and consider only one possibility
      if(!results_container.chk_real_neu_pz && index_neutrino_sol!=0) break;
      
      Set_Neutrino();
      
      //loopsd over jet permutation
      do
	{
	  //Print_Permutation();
	  
	  results.Reset();
	
	  //each t quark should contain at least one 
	  if(!BJet_Assignment_Cut())
            {
              results.cut = CUT_RESULT::B_ASSIGNMENT;
              Save_Permutation(true);
          
              continue;
            }

	  //check repetition
	  //For KF, permutations within W cadidates are meaningless. So reject these permutation to save time
	  if(!Check_Repetition())
	    {
	      results.cut = CUT_RESULT::REPETITION;
	      Save_Permutation(true);
	      
	      continue;
	    }
	  
	  Set_Jets();

	  if(!Pre_Kinematic_Cut())
	    {
	      results.cut = CUT_RESULT::PRE_KINEMATIC;
	      Save_Permutation(true);
	  
	      continue;
	    }

	  //contruct TFitConstraint objects
	  Set_Constraints();

	  //register things to fitter
	  Set_Fitter();

	  //fit
	  results.status = fitter->fit();
	  
	  if(!Quality_Cut())
            {
              results.cut = CUT_RESULT::QUALITY;
	      Save_Permutation(true);

              continue;
            }

	  results.cut = CUT_RESULT::PASS;
	  
	  //Save
 	  Save_Results();

  	}while(End_Permutation());//loop over neutrino pz solution

      index_neutrino_sol++;
    }//loop over neutrino pz solution
  
  //find best permutation
  Find_Best_Permutation();
  
  return;
}//void TKinFitterDriver::Scan()

//////////

void TKinFitterDriver::Set_Objects(vector<Jet>& _vec_jet, vector<float>& _vec_resolution_pt, vector<bool>& _vec_btag, Lepton& _lepton, Particle& _met)
{
  Clear();
  
  vec_jet = _vec_jet;
  vec_resolution_pt = _vec_resolution_pt;
  vec_btag = _vec_btag;
  lepton = _lepton;
  met = _met;
  
  n_jet = vec_jet.size();
  int n_vec_resolution_pt = vec_resolution_pt.size();
  int n_vec_btag = vec_btag.size();
  if(n_jet!=n_vec_resolution_pt || n_jet!=n_vec_btag)
    {
      cerr << "The size of two vectors are not same. Check it." << endl;
      cerr << "Size of jet vector = " << n_jet << ", Size of btagger vector = " << n_vec_btag << endl;
      
      exit(1);
    }
  
  for(int i=0; i<n_jet; i++)
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
   
   n_jet = 0;
   n_btag = 0;
   
   return;
 }//void TKinFitterDriver::Clear()

//////////

void TKinFitterDriver::Find_Best_Permutation()
{
  const float chi2_max = 99999;
  float chi2_best = chi2_max;
  int i_best = 0;
  for(unsigned int i=0; i<results_container.vec_results.size(); i++)
    {
      Results results = results_container.vec_results.at(i);
           
      if(results.chi2<chi2_best && results.cut==CUT_RESULT::PASS)
	{
	  i_best = i;
	  chi2_best = results.chi2;
	}
    }
  
  if(chi2_best<0.99*chi2_max)//arbitary number to avoid ambiguity of float comparison
    {
      results_container.status = true;
      results_container.best_chi2 = chi2_best;
      
      Results results = results_container.vec_results.at(i_best);
      
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
  
  //cout << results_container.best_chi2 << " " <<  results_container.best_fitted_had_t_mass << " " <<  results_container.best_fitted_had_w_mass << endl;
  
  return;
}//void TKinFitterDriver::Find_Best_Permutation()

//////////

bool TKinFitterDriver::Pre_Kinematic_Cut()
{
  bool jet_pt_cut = false;
  if(jet_lep_t_b.Pt()>30) jet_pt_cut = true;//need to be tuned more
  if(!jet_pt_cut) return false;

  TLorentzVector had_t = jet_had_t_b + jet_w_u + jet_w_d;
  TLorentzVector lep_t = jet_lep_t_b + lepton + neutrino;

  bool tt_angular_cut = (TMath::Abs(had_t.DeltaPhi(lep_t) > 1.5)) ? true : false;
  if(!tt_angular_cut) return false;

  bool had_t_mass_cut = false;
  double had_t_mass = had_t.M();
  if(n_jet==4 && 100<had_t_mass && had_t_mass<240) had_t_mass_cut = true;
  else if(n_jet==5 && 120<had_t_mass && had_t_mass<220) had_t_mass_cut = true;
  else if(n_jet>=6 && 140<had_t_mass && had_t_mass<200) had_t_mass_cut = true;
  if(!had_t_mass_cut) return false;

  bool lep_t_mass_cut = false;
  TLorentzVector partial_lep_t = jet_lep_t_b + lepton;
  if(partial_lep_t.M()<150) lep_t_mass_cut = true;
  if(!lep_t_mass_cut) return false;

  return true;
}//bool TKinFitterDriver::Pre_Kinematic_Cut()

//////////

void TKinFitterDriver::Print_Permutation()
{
  cout << "Permutation: ";
  for(unsigned int i=0; i<vec_permutation.size(); i++){ cout << vec_permutation.at(i) << ", "; }
  cout << endl;

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
  if(results.status==0 && chi2_jet_had_t_b<9 && chi2_jet_w_u<9 &&  chi2_jet_w_d<9 && chi2_jet_lep_t_b<9) return true;
  else return false;

  //dummy
  return true;
}//bool TKinFitterDriver::Quality_Cut()

//////////

void TKinFitterDriver::Save_Permutation(const bool& push)
{
  results.index_neutrino_sol = index_neutrino_sol;
  results.vec_permutation = vec_permutation;
  
  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);
      
      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B) results.index_had_t_b = i;
      else if(jet_assignment==JET_ASSIGNMENT::W_U) results.index_w_u = i;
      else if(jet_assignment==JET_ASSIGNMENT::W_D) results.index_w_d = i;
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B) results.index_lep_t_b = i;
    }
  
  if(push) results_container.vec_results.push_back(results);

  return;
}//void TKinFitterDriver::Save_Permutation()

//////////

void TKinFitterDriver::Save_Results()
{
  Save_Permutation();
  
  //initial objects
  TLorentzVector initial_had_t = jet_had_t_b + jet_w_u + jet_w_d;
  TLorentzVector initial_had_w = jet_w_u + jet_w_d;
  TLorentzVector initial_lep_t = jet_lep_t_b + lepton + neutrino;
  TLorentzVector initial_lep_w = lepton + neutrino;

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
  
  results.chi2 = Calc_Chi2();
  
  results.initial_had_t_mass = initial_had_t.M();
  results.initial_had_w_mass = initial_had_w.M();
  results.initial_lep_t_mass = initial_lep_t.M();
  results.initial_lep_w_mass = initial_lep_w.M();

  results.fitted_had_t_mass = fitted_had_t.M();
  results.fitted_had_w_mass = fitted_had_w.M();
  results.fitted_lep_t_mass = fitted_lep_t.M();
  results.fitted_lep_w_mass = fitted_lep_w.M();
 
  results_container.vec_results.push_back(results);
  
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
  fitter->addMeasParticle(fit_jet_w_d);
  fitter->addMeasParticle(fit_jet_w_d);
  fitter->addMeasParticle(fit_jet_lep_t_b);
  for(auto& fit_jet_extra : vec_fit_jet_extra) fitter->addMeasParticle(fit_jet_extra);
  
  //register unmeasure particles
  fitter->addUnmeasParticle(fit_neutrino);

  //register constraints
  fitter->addConstraint(constraint_had_t_mgaus);
  fitter->addConstraint(constraint_lep_t_mgaus);
  fitter->addConstraint(constraint_had_w_mgaus);
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
      
      TLorentzVector jet = vec_jet.at(i);
      float pt = jet.Pt();
      float jer = vec_resolution_pt.at(i);
      
      //for easy handling
      if(jet_assignment==JET_ASSIGNMENT::HAD_T_B)
	{
	  results.index_had_t_b = i;
	  
	  jet_had_t_b = jet;
	  error_jet_had_t_b(0, 0) = jer*jer*pt*pt;
	} 
      else if(jet_assignment==JET_ASSIGNMENT::W_U)
	{
	  results.index_w_u = i;
	  
	  jet_w_u = jet;
	  error_jet_w_u(0, 0) = jer*jer*pt*pt;
	}
      else if(jet_assignment==JET_ASSIGNMENT::W_D) 
	{
	  results.index_w_d = i;
	  
	  jet_w_d = jet;
	  error_jet_w_d(0, 0) = jer*jer*pt*pt;
	}
      else if(jet_assignment==JET_ASSIGNMENT::LEP_T_B) 
	{
	  results.index_lep_t_b = i;
	  
	  jet_lep_t_b = jet;
	  error_jet_lep_t_b(0, 0) = jer*jer*pt*pt;
	}
      else if(jet_assignment==JET_ASSIGNMENT::OTHERS)
	{
	  vec_jet_extra.push_back(jet);
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

void TKinFitterDriver::Set_Neutrino()
{
  float px = met.Px();
  float py = met.Py();
  float pz = results_container.neutrino_pz_sol[index_neutrino_sol];
  
  neutrino.SetPxPyPzE(px, py, pz, TMath::Sqrt(px*px+py*py+pz*pz));
  
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
  
  double k = TMath::Power(W_MASS, 2.)/2.0 - lepton_mass*lepton_mass/2.0 + lepton.Px()*met.Px() + lepton.Py()*met.Py();
  double a = TMath::Power(lepton.Px(), 2.0) + TMath::Power(lepton.Py(), 2.0);
  double b = -2*k*lepton.Pz();
  double c = TMath::Power(lepton.Pt(), 2.0)*TMath::Power(met.Pt(), 2.0) - TMath::Power(k, 2.0);

  double determinant = TMath::Power(b, 2.0) - 4*a*c; 

  //real solution
  if(determinant>=0){
    results_container.neutrino_pz_sol[0] = (-b + TMath::Sqrt(determinant))/(2*a);
    results_container.neutrino_pz_sol[1] = (-b - TMath::Sqrt(determinant))/(2*a);
    
    results_container.chk_real_neu_pz = true;
  }
  //complex solution
  else{
    results_container.neutrino_pz_sol[0] = -b/(2*a);
    //Resol_Neutrino_Pt();
    results_container.chk_real_neu_pz = false;
  }
  
  return;
}//void TKinFitterDriver::Sol_Neutrino_Pz()
 
//////////
