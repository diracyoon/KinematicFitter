#ifndef Jet_h
#define Jet_h

#include "Particle.h"

class Jet: public Particle {
public:

  Jet();
  ~Jet();

  void SetArea(double area);
  void SetGenFlavours(double pf, double hf);
  void SetTaggerResults(std::vector<double> ds);
  void SetEnergyFractions(double cH, double nH, double nEM, double cEM);
  void SetMultiplicities(double cM, double nM);
  void SetPileupJetId(double v);


  enum Tagger{
    CSVv2,
    DeepCSV, DeepCvsL, DeepCvsB,
    CvsL, CvsB,
    DeepFlavour_b, DeepFlavour_bb, DeepFlavour_lepb, DeepFlavour_c, DeepFlavour_uds, DeepFlavour_g,
  };
  enum WP{
    Loose, Medium, Tight
  };
  double GetTaggerResult(Tagger tg);
  bool IsTagged(Tagger tg, WP wp);

private:

  double  j_area;
  int j_partonFlavour;
  int j_hadronFlavour;
  double  j_CSVv2;
  double  j_DeepCSV;
  double  j_DeepFlavour_b;
  double  j_DeepFlavour_bb;
  double  j_DeepFlavour_lepb;
  double  j_DeepFlavour_c;
  double  j_DeepFlavour_uds;
  double  j_DeepFlavour_g;
  double  j_CvsL;
  double  j_CvsB;
  double  j_DeepCvsL;
  double  j_DeepCvsB;
  double  j_chargedHadronEnergyFraction;
  double  j_neutralHadronEnergyFraction;
  double  j_neutralEmEnergyFraction;
  double  j_chargedEmEnergyFraction;
  int j_chargedMultiplicity;
  int j_neutralMultiplicity;
  double j_PileupJetId;
};

Jet::Jet() : Particle() {
  j_area=-999.;
  j_partonFlavour=-999;
  j_hadronFlavour=-999;
  j_CSVv2=-999.;
  j_DeepCSV=-999.;
  j_DeepFlavour_b=-999;
  j_DeepFlavour_bb=-999;
  j_DeepFlavour_lepb=-999;
  j_DeepFlavour_c=-999;
  j_DeepFlavour_uds=-999;
  j_DeepFlavour_g=-999;
  j_CvsL=-999.;
  j_CvsB=-999.;
  j_DeepCvsL=-999.;
  j_DeepCvsB=-999.;
  j_chargedHadronEnergyFraction=-999.;
  j_neutralHadronEnergyFraction=-999.;
  j_neutralEmEnergyFraction=-999.;
  j_chargedEmEnergyFraction=-999.;
  j_chargedMultiplicity=-999;
  j_neutralMultiplicity=-999;
  j_PileupJetId=-999.;
}

Jet::~Jet(){

}

void Jet::SetArea(double area){
  j_area = area;
}
void Jet::SetGenFlavours(double pf, double hf){
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void Jet::SetTaggerResults(std::vector<double> ds){
  j_CSVv2             = ds.at(0);
  j_DeepCSV           = ds.at(1);
  j_DeepCvsL          = ds.at(2);
  j_DeepCvsB          = ds.at(3);
  j_DeepFlavour_b     = ds.at(4);
  j_DeepFlavour_bb    = ds.at(5);
  j_DeepFlavour_lepb  = ds.at(6);
  j_DeepFlavour_c     = ds.at(7);
  j_DeepFlavour_uds   = ds.at(8);
  j_DeepFlavour_g     = ds.at(9);
  j_CvsL              = ds.at(10);
  j_CvsB              = ds.at(11);
}
void Jet::SetEnergyFractions(double cH, double nH, double nEM, double cEM){
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
}
void Jet::SetMultiplicities(double cM, double nM){
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}
void Jet::SetPileupJetId(double v){
  j_PileupJetId = v;
}

double Jet::GetTaggerResult(Tagger tg){
  if(tg==CSVv2) return j_CSVv2;
  else if(tg==DeepCSV) return j_DeepCSV;
  else if(tg==DeepFlavour_b) return j_DeepFlavour_b;
  else if(tg==DeepFlavour_bb) return j_DeepFlavour_bb;
  else if(tg==DeepFlavour_lepb) return j_DeepFlavour_lepb;
  else if(tg==DeepFlavour_c) return j_DeepFlavour_c;
  else if(tg==DeepFlavour_uds) return j_DeepFlavour_uds;
  else if(tg==DeepFlavour_g) return j_DeepFlavour_g;
  else if(tg==CvsL) return j_CvsL;
  else if(tg==CvsB) return j_CvsB;
  else if(tg==DeepCvsL) return j_DeepCvsL;
  else if(tg==DeepCvsB) return j_DeepCvsB;
  else{
    cout << "[Jet::GetTaggerResult] ERROR; Wrong tagger" << endl;
    return -999;
  }
}


bool Jet::IsTagged(Tagger tg, WP wp){

  double tr = GetTaggerResult(tg);

  if(tg==CSVv2){
    if(wp==Loose) return (tr>0.5803);
    else if(wp==Medium) return (tr>0.8838);
    else if(wp==Tight) return (tr>0.9693);
    else return false;
  }
  else if(tg==DeepCSV){
    if(wp==Loose) return (tr>0.1522);
    else if(wp==Medium) return (tr>0.4941);
    else if(wp==Tight) return (tr>0.8001);
    else return false;
  }
/* FIXME
  else if(tg==DeepFlavour){
    if(wp==Loose) return (tr>0.0574);
    else if(wp==Medium) return (tr>0.4318);
    else if(wp==Tight) return (tr>0.9068);
    else return false;
  }  
*/
  else if(tg==CvsL){
    if(wp==Loose) return (tr>-0.53);
    else if(wp==Medium) return (tr>0.07);
    else if(wp==Tight) return (tr>0.87);
    else return false;
  }  
  else if(tg==CvsB){
    if(wp==Loose) return (tr>-0.26);
    else if(wp==Medium) return (tr>-0.10);
    else if(wp==Tight) return (tr>-0.3);
    else return false;
  }  
  else if(tg==DeepCvsL){
    if(wp==Loose) return (tr>0.05);
    else if(wp==Medium) return (tr>0.15);
    else if(wp==Tight) return (tr>0.8);
    else return false;
  }  
  else if(tg==DeepCvsB){
    if(wp==Loose) return (tr>0.33);
    else if(wp==Medium) return (tr>0.28);
    else if(wp==Tight) return (tr>0.1);
    else return false;
  }  
  else{
    cout << "[Jet::IsTagged] ERROR; Wrong tagger" << endl;
    return false;
  }


}

#endif
