#include "TKinFitterDriver.h"

ClassImp(TKinFitterDriver);

//////////

TKinFitterDriver::TKinFitterDriver(int _data_year)
{
  data_year = _data_year;
}//TKinFitterDriver::TKinFitterDriver(int _DataYear)

//////////

bool TKinFitterDriver::End_Permutation()
{
  return next_permutation(vec_permutation.begin(), vec_permutation.end());
}//void TKinFitterDriver::End_Permutation()

//////////

void TKinFitterDriver::Scan()
{
  if(n_jet<4 || n_btag<2) return;

  Sol_Neutrino_Pz();

  //loop over neutrino pz solution
  for(int i=0; i<2; i++)
    {
      //if neutrino pz solution gives complex, take only real part and consider only one possibility
      if(!chk_real_neu_pz && i!=0) continue;
      Set_Neutrino(i);

      //loop over jet permutation
      do
	{
	  Set_Current_Permutation();
	  
	  if(!BJet_Assignment_Cut()) continue;
	  if(!Pre_Kinematic_Cut()) continue;
	  
	  Fit();
	}while(End_Permutation());//loop over neutrino pz solution
    }//loop over neutrino pz solution

  return;
}//void TKinFitterDriver::Scan()

//////////

void TKinFitterDriver::Set_Objects(vector<Jet>& _vec_jet, vector<bool>& _vec_btag, Lepton& _lepton, Particle& _met)
{
  Clear();

  vec_jet = _vec_jet;
  vec_btag = _vec_btag;
  lepton = _lepton;
  met = _met;

  n_jet = vec_jet.size();
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

bool BJet_Assignment_Cut()
{
  //condition develeped for the CH analysis can't be used for the Vcb analysis
  //this part will be developed later...
 
  return true;
}//bool BJet_Assignment_Cut()

//////////

void TKinFitterDriver::Clear()
{
  vec_jet.clear();
  vec_btag.clear();
  vec_permutation.clear();
  
  n_jet = 0;
  n_btag = 0;

  return;
}//void TKinFitterDriver::Clear()

//////////

void TKinFitterDriver::Fit()
{
  return;
}//void TKinFitterDriver::Fit()

//////////

void TKinFitterDriver::Resol_Neutrino_Pz()
{
  return;
}//void TKinFitterDriver::Resol_Neutrino_Pz()

//////////

bool TKinFitterDriver::Pre_Kinematic_Cut()
{
  bool jet_pt_cut = false;
  if(jet_lep_t_b.Pt()>30) jet_pt_cut = true;//need to be tuned more
  if(!jet_pt_cut) return false;

  TLorentzVector had_t = jet_had_t_b + jet_w_u + jet_w_d;
  TLorentzVector lep_t = jet_lep_t_b + lepton + neutrino; 

  bool tt_angular_cut = (had_t.DeltaPhi(lep_t) > 1.5) ? true : false;
  if(!tt_angular_cut) return false;
  
  bool had_t_mass_cut = false;
  double had_t_mass = had_t.M();
  if(n_jet==4 && 100<had_t_mass && 240<had_t_mass) had_t_mass_cut = true;
  else if(n_jet==5 && 120<had_t_mass && 220<had_t_mass) had_t_mass_cut = true;
  else if(n_jet>=6 && 140<had_t_mass && 200<had_t_mass) had_t_mass_cut = true; 
  if(!had_t_mass_cut) return false;

  bool lep_t_mass_cut = false;
  TLorentzVector partial_lep_t = jet_lep_t_b + lepton;
  if(partial_lep_t.M()<150) lep_t_mass_cut = true;
  if(!lep_t_mass_cut) return false;


  return true;
}//bool TKinFitterDriver::Pre_Kinematic_Cut()

//////////

void TKinFitterDriver::Set_Current_Permutation()
{
  for(unsigned int i=0; i<vec_permutation.size(); i++)
    {
      float pt, eta, phi;
      JET_ASSIGNMENT jet_assignment = vec_permutation.at(i);
    
      //for easy handling
      if(jet_assignment==HAD_T_B) jet_had_t_b = vec_jet.at(jet_assignment);
      else if(jet_assignment==W_U) jet_w_u = vec_jet.at(jet_assignment);
      else if(jet_assignment==W_D) jet_w_d = vec_jet.at(jet_assignment);
      else if(jet_assignment==LEP_T_B) jet_lep_t_b = vec_jet.at(jet_assignment);
    }

  return;
}//void TKinFitterDriver::Set_Current_Permutation()

//////////

void TKinFitterDriver::Set_Lepton()
{
  //fit_lepton = new TFitParticlePt("Lepton", "Lepton", );

  return;
}//void TKinFitterDriver::Set_Lepton()

//////////

void TKinFitterDriver::Set_Neutrino(const int& index)
{
  float px = met.Px();
  float py = met.Py();
  float pz = neutrino_pz_sol[index];
  
  neutrino.SetPxPyPzE(px, py, pz, TMath::Sqrt(px*px+py*py+pz*pz));

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
    neutrino_pz_sol[0] = (-b + TMath::Sqrt(determinant))/(2*a);
    neutrino_pz_sol[1] = (-b - TMath::Sqrt(determinant))/(2*a);
    chk_real_neu_pz = true;
  }
  //complex solution
  else{
    neutrino_pz_sol[0] = -b/(2*a);
    //Resol_Neutrino_Pt();
    chk_real_neu_pz = false;
  }

  return;
}//void TKinFitterDriver::Sol_Neutrino_Pz()

//////////
