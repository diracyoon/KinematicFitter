// Classname: TFitConstraintMGaus
// Author: Jan E. Sundermann, Verena Klose (TU Dresden)      


//________________________________________________________________
// 
// TFitConstraintMGaus::
// --------------------
//
// Fit constraint: mass conservation ( m_i - m_j - alpha * MassConstraint == 0 )
//

#include <iostream>
#include <iomanip>
#include "TFitConstraintMGaus_mod.h"

ClassImp(TFitConstraintMGaus_mod)

//#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TClass.h"


//----------------
// Constructor --
//----------------
TFitConstraintMGaus_mod::TFitConstraintMGaus_mod()
  : TFitConstraintM() 
{

  init();

}

TFitConstraintMGaus_mod::TFitConstraintMGaus_mod(std::vector<TAbsFitParticle*>* ParList1,
					 std::vector<TAbsFitParticle*>* ParList2, 
					 Double_t Mass,
					 Double_t Width)
  : TFitConstraintM(ParList1, ParList2, Mass ) 
{

  init();
  setMassConstraint( Mass, Width );

}

TFitConstraintMGaus_mod::TFitConstraintMGaus_mod(const TString &name, const TString &title,
					 std::vector<TAbsFitParticle*>* ParList1,
					 std::vector<TAbsFitParticle*>* ParList2, 
					 Double_t Mass,
					 Double_t Width)
  : TFitConstraintM( name, title, ParList1, ParList2, Mass )
{

  init();
  setMassConstraint( Mass, Width );

}

void
TFitConstraintMGaus_mod::init() {

  _nPar = 1;
  _iniparameters.ResizeTo(1,1);
  _iniparameters(0,0) = 1.;
  _parameters.ResizeTo(1,1);
  _parameters = _iniparameters;

}

//--------------
// Destructor --
//--------------
TFitConstraintMGaus_mod::~TFitConstraintMGaus_mod() {

}

//--------------
// Operations --
//--------------

void TFitConstraintMGaus_mod::setMassConstraint(Double_t Mass, Double_t Width) { 
  
  _TheMassConstraint = Mass;
  _width = Width;
  setCovMatrix( 0 );
  if(!Mass) std::cout
    << "Error occured!\n"
    << "Object type : TFitConstraintMGaus_mod\n"
    << "Object name : " << GetName() << "\n"
    << "Object title: " << GetTitle() << "\n"
    << "Mass of 0 GeV not supported, please choose a larger mass!\n";
  _covMatrix(0,0) = (Width*Width) / (Mass * Mass);

}

Double_t TFitConstraintMGaus_mod::getInitValue() {
  // Get initial value of constraint (before the fit)

  Double_t InitValue = 
    CalcMass( &_ParList1, true ) - 
    CalcMass( &_ParList2, true ) - 
    _iniparameters(0,0) * _TheMassConstraint;

  return InitValue;

}

Double_t TFitConstraintMGaus_mod::getCurrentValue() {
  // Get value of constraint after the fit

  Double_t CurrentValue =
    CalcMass(&_ParList1,false) - 
    CalcMass(&_ParList2,false) - 
    _parameters(0,0)*_TheMassConstraint;

  return CurrentValue;

}

//BHO
Double_t TFitConstraintMGaus_mod::getChi2() {

  Double_t diff =
    CalcMass(&_ParList1,false) - 
    CalcMass(&_ParList2,false) - 
    _TheMassConstraint;
  Double_t chi2 = (diff*diff)/(_width*_width);

  return chi2;

}


TMatrixD* TFitConstraintMGaus_mod::getDerivativeAlpha() { 
  // Calculate dF/dAlpha = -1 * M

  TMatrixD* DerivativeMatrix = new TMatrixD(1,1);
  DerivativeMatrix->Zero();

  (*DerivativeMatrix)(0,0) = -1. * _TheMassConstraint;

  return DerivativeMatrix;

}

TString TFitConstraintMGaus_mod::getInfoString() {
  // Collect information to be used for printout

  std::stringstream info;
  info << std::scientific << std::setprecision(6);

  info << "__________________________" << std::endl
       << std::endl;
  info << "OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << std::endl;

  info << "initial value: " << getInitValue() << std::endl;
  info << "current value: " << getCurrentValue() << std::endl;
  info << "mean mass: " << _TheMassConstraint << std::endl;
  info << "width: " << _width << std::endl;
  info << "initial mass: " << _iniparameters(0,0)*_TheMassConstraint  << std::endl;
  info << "current mass: " << _parameters(0,0)*_TheMassConstraint  << std::endl;

  return info.str();

}

void TFitConstraintMGaus_mod::print() {
  // Print constraint contents

  //edm::LogVerbatim("KinFitter") << this->getInfoString();

}

