#include <cassert>
#include <iostream>
#include <math.h>
#include <complex>

#include <Eigen/Eigenvalues>

#include "PMNS_OQS.h"

using namespace std;

using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_OQS::PMNS_OQS()
  : PMNS_Fast(), fPhi(), fR(), fRt(),
    fRho(3, vectorC(3, 0)), fD(9, vectorC(9, 0)),
    fM(9, vectorC(9, 0)), fMd(9, 9),
    fHeff(3, vectorC(3, 0)), fHGM(9, vectorC(9, 0)),
    fMEvec(9, 9)
{
  SetStdPath();
  InitializeVectors();
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_OQS::~PMNS_OQS() {}

void PMNS_OQS::InitializeVectors()
{
  SetPhi(1, 0);
  SetPhi(2, 0);
  fEval = vectorC(9, 0);
  fEvec = matrixC(9, vectorC(9, (0, 0)));
}

// set Heff in vacuum-mass basis
void PMNS_OQS::SetHeff(NuPath p){

  cout << "setHeff\n";
  
  double Ve = 0;
  
  // Set the neutrino path                                                            
  SetCurPath(p);
  
  double rho = fPath.density;
  double zoa = fPath.zoa;
  
  double kr2GNe = kK2 * M_SQRT2 * kGf * rho * zoa; // Electron matter potential in eV
  
  if(fIsNuBar == 0){
    Ve = kr2GNe*2.;
  } else {
    Ve = -kr2GNe*2.;
  }
  
  double s12 = sin(fTheta[0][1]);
  double s13 = sin(fTheta[0][2]);
  double s23 = sin(fTheta[1][2]);
  
  double c12 = cos(fTheta[0][1]);
  double c13 = cos(fTheta[0][2]);
  double c23 = cos(fTheta[1][2]);

  complexD idelta(0.0, fDelta[0][2]);
  if (fIsNuBar){
    idelta = conj(idelta);
  }
  complexD iphi1(0.0, fPhi[0]);
  complexD iphi2(0.0, fPhi[1]);

  cout << "setHeff2\n";

  fHeff[0][0] = c12*c12 * c13*c13 * Ve;
  fHeff[0][1] = c12 * c13*c13 * exp(iphi1) * s12 * Ve;
  fHeff[0][2] = c12 * c13 * exp(-idelta + iphi2) * s13 * Ve;

  fHeff[1][0] = c12 * c13*c13 * exp(-iphi1) * s12 * Ve;
  fHeff[1][1] = c13*c13 * s12*s12 * Ve;
  fHeff[1][2] = c13 * exp(-idelta - iphi1 + iphi2) * s12 * s13 * Ve;

  fHeff[2][0] = c12 * c13 * exp(idelta - iphi2) * s13 * Ve;
  fHeff[2][1] = c13 * exp(idelta + iphi1 - iphi2) * s12 * s13 * Ve;
  fHeff[2][2] = s13*s13 * Ve;  

  // add mass terms
  fHeff[1][1] += fDm[1] / (2 * kGeV2eV * fEnergy);
  fHeff[2][2] += fDm[2] / (2 * kGeV2eV * fEnergy);

  cout << "Heff:\n";
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      cout << fHeff[i][j] << " ";
    }
    cout << endl;
  }

}


// Set Heff in Gell-Mann basis: set only right-triangle as the left part is -right
void PMNS_OQS::SetHGM()
{

  //double k = 1./2.;
  double k = 1.;

  cout << "test Heff in setHGM " << fHeff[0][1] << " real " << real(fHeff[0][1]) << " imag " << imag(fHeff[0][1]) << endl;  
  
  fHGM[1][2] =  k * (fHeff[0][0] - fHeff[1][1]);
  fHGM[1][3] =  2. * k * imag(fHeff[0][1]);
  fHGM[1][4] = -k * imag(fHeff[1][2]);
  fHGM[1][5] = -k * real(fHeff[1][2]);
  fHGM[1][6] = -k * imag(fHeff[0][2]);
  fHGM[1][7] = -k * real(fHeff[0][2]);

  fHGM[2][3] =  2. * k * real(fHeff[0][1]);
  fHGM[2][4] =  k * real(fHeff[1][2]);
  fHGM[2][5] = -k * imag(fHeff[1][2]);
  fHGM[2][6] = -k * real(fHeff[0][2]);
  fHGM[2][7] =  k * imag(fHeff[0][2]);

  fHGM[3][4] = -k * imag(fHeff[0][2]);
  fHGM[3][5] = -k * real(fHeff[0][2]);
  fHGM[3][6] =  k * imag(fHeff[1][2]);
  fHGM[3][7] =  k * real(fHeff[1][2]);

  fHGM[4][5] =  k * (fHeff[0][0] - fHeff[2][2]);
  fHGM[4][6] = -k * imag(fHeff[0][1]);
  fHGM[4][7] =  k * real(fHeff[0][1]);
  fHGM[4][8] =  sqrt(3.) * k * imag(fHeff[0][2]);

  fHGM[5][6] = -k * real(fHeff[0][1]);
  fHGM[5][7] = -k * imag(fHeff[0][1]);
  fHGM[5][8] =  sqrt(3.) * k * real(fHeff[0][2]);
  
  fHGM[6][7] =  k * (fHeff[1][1] - fHeff[2][2]);
  fHGM[6][8] =  sqrt(3.) * k * imag(fHeff[1][2]);
  
  fHGM[7][8] = sqrt(3.) * k * real(fHeff[1][2]);

  cout << "HGM: \n";
  for(int i = 1; i < 9; ++i){
    for(int j = 1; j < 9; ++j){

      if(fHGM[j][i] == complexD(0, 0)) fHGM[j][i] = -fHGM[i][j];
      
      cout << fHGM[i][j] << " ";
    }
    cout << endl;
  }
		     
}


void PMNS_OQS::SetDissipatorElement(int i, int j, double val)
{
  fD[i][j] = val;
  fD[j][i] = conj(fD[i][j]);
}


void PMNS_OQS::SetM()
{

  for(int k = 0; k < 9; ++k){
    for(int j = 0; j < 9; ++j){
      fM[k][j] = fHGM[k][j] + fD[k][j];
    }
  }  
  
}


void PMNS_OQS::SetPhi(int i, double val)
{
  fPhi[i - 1] = val;  
}


template <typename T> void PMNS_OQS::SolveEigenSystem()
{

  cout << "M to diagonalize: \n";

  for(int i = 0; i < 9; ++i){
    for(int j = 0; j < 9; ++j){
      fMd(i, j) = fM[i][j];

      cout << fMd(i,j) << " "; 
    }
    cout << "\n";
  }
  
  Eigen::Ref<T> M(fMd);
  
  Eigen::ComplexEigenSolver<T> eigensolver(M);

  // Fill fEval and fEvec vectors from GSL objects
  for (int i = 0; i < 9; i++) {
    fEval[i] = eigensolver.eigenvalues()(i);
    for (int j = 0; j < 9; j++) {
      fEvec[i][j] = eigensolver.eigenvectors()(i, j);
    }
  }

}


void PMNS_OQS::Diagonalise()
{

  SolveEigenSystem<Eigen::MatrixXcd>();
  
  // Mark eigensystem as solved
  fGotES = true;
  
  // Fill cache if activated
  FillCache();
}


// rotate fRho into vacuum mass basis or flavour basis
void PMNS_OQS::RotateState(bool to_mass)
{

  matrixC UM(3, vectorC(3, 0));

  double s12 = sin(fTheta[0][1]);
  double s13 = sin(fTheta[0][2]);
  double s23 = sin(fTheta[1][2]);
  
  double c12 = cos(fTheta[0][1]);
  double c13 = cos(fTheta[0][2]);
  double c23 = cos(fTheta[1][2]);

  complexD idelta(0.0, fDelta[0][2]);
  if (fIsNuBar){
    idelta = conj(idelta);
  }
  complexD iphi1(0.0, fPhi[0]);
  complexD iphi2(0.0, fPhi[1]);

  UM[0][0] =  c12*c13;
  UM[0][1] =  s12*c13*exp(iphi1);
  UM[0][2] =  s13*exp(-idelta + iphi2);

  UM[1][0] = -exp(-iphi1) * (s12*c23 + c12*s13*s23*exp(idelta));
  UM[1][1] =  c12 * c23 - s12 * s13 * s23 * exp(idelta);
  UM[1][2] =  s23 * c13 * exp(iphi2 - iphi1);

  UM[2][0] =  exp(-iphi2) * (s12*s23 - c12*c23*s13*exp(idelta));
  UM[2][1] = -exp(iphi1 - iphi2) * (c12*s23 + s12*s13*c23*exp(idelta));
  UM[2][2] =  c13 * c23;

  matrixC mult(3, vectorC(3, 0));
  
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      mult[i][j] = 0;
      for (int k = 0; k < 3; k++) {
	if (to_mass)
	  mult[i][j] += fRho[i][k] * UM[k][j];
	else
	  mult[i][j] += fRho[i][k] * conj(UM[j][k]);
      }
    }
   }
  
   for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      fRho[i][j] = 0;
      for (int k = 0; k < 3; k++) {
	if (to_mass)
	  fRho[i][j] += conj(UM[k][i]) * mult[k][j];
	else
          fRho[i][j] += UM[i][k] * mult[k][j];
      }
      if (j > i) fRho[j][i] = conj(fRho[i][j]);
    }
   }

   cout << "fRho in SU3 basis: " << endl;

   for (int i = 0; i < 3; i++) {
     for (int j = 0; j < 3; j++) {
       cout << fRho[i][j] << " ";
     }
     cout << endl;
   }
}


void PMNS_OQS::ChangeBaseToGM()
{
  
  cout << "BaseToGM: \n" << fRho[0][0] << " " << fRho[0][1] << " " << fRho[0][2]
       << " " << fRho[1][0] << " " << fRho[1][1] << " " << fRho[1][2]
       << " " << fRho[2][0] << " " << fRho[2][1] << " " << fRho[2][2] << endl;

  double k = 2.;

  fR[0] =  k * (fRho[0][0] + fRho[1][1] + fRho[2][2]) / sqrt(6.);
  fR[1] =  k * real(fRho[0][1]);
  fR[2] = -k * imag(fRho[0][1]);
  fR[3] =  k * (fRho[0][0] - fRho[1][1]) / 2.;
  fR[4] =  k * real(fRho[0][2]);
  fR[5] = -k * imag(fRho[0][2]);
  fR[6] =  k * real(fRho[1][2]);
  fR[7] = -k * imag(fRho[1][2]);
  fR[8] =  k * (fRho[0][0] + fRho[1][1] - 2.*fRho[2][2]) / (2.*sqrt(3.));

  cout << "Rmu " << fR[0] << " " << fR[1] << " " << fR[2] << endl;  
}


void PMNS_OQS::ChangeBaseToSU3()
{

  double k = 1./2.;
  
  fRho[0][0] = k * (sqrt(2./3.)*fRt[0] + fRt[3] + fRt[8]/sqrt(3.));
  fRho[0][1] = k * (fRt[1] - complexD(0.0,1.0) * fRt[2]);
  fRho[0][2] = k * (fRt[4] - complexD(0.0,1.0) * fRt[5]);

  fRho[1][0] = k * (fRt[1] + complexD(0.0,1.0) * fRt[2]);
  fRho[1][1] = k * (sqrt(2./3.)*fRt[0] - fRt[3] + fRt[8]/sqrt(3.));
  fRho[1][2] = k * (fRt[6] - complexD(0.0,1.0) * fRt[7]);

  fRho[2][0] = k * (fRt[4] + complexD(0.0,1.0) * fRt[5]);
  fRho[2][1] = k * (fRt[6] + complexD(0.0,1.0) * fRt[7]);
  fRho[2][2] = k * (1. / sqrt(3.) * (sqrt(2.)*fRt[0] - 2.*fRt[8]));

  cout << "BaseToSU3: \n" << fRho[0][0] << " " << fRho[0][1] << " " << fRho[0][2] 
       << "\n" << fRho[1][0] << " " << fRho[1][1] << " " << fRho[1][2]
       << "\n" << fRho[2][0] << " " << fRho[2][1] << " " << fRho[2][2] << endl;	

  cout << "ABS: " << abs(fRho[0][0]) << " " << abs(fRho[1][1]) << " " << abs(fRho[2][2]) << endl;  
}


//.............................................................................
///
/// Propagate the current neutrino state through a given path
///
/// @param p - A neutrino path segment
///
void PMNS_OQS::PropagatePath(NuPath p)
{

  cout << "start of propagatepath\n";

  // define Heff in vacuum mass basis
  SetHeff(p);

  // define Heff in Gell-Mann basis
  SetHGM();

  // Define M matrix as Heff_GM + Dissipator
  SetM();

  // Solve for eigensystem
  Diagonalise();

  cout << "MATRIX WITH EIGENVECTORS: \n";

  for(int i = 0; i < 9; ++i){
    for(int j = 0; j < 9; ++j){
      fMEvec(i, j) = fEvec[i][j];

      cout << fEvec[i][j] << " ";
    }
    cout << endl;
  }

  cout << "EIGENVALUES: \n";
  //  for(int i = 0; i < 9; ++i){
  for(int i = 0; i < 9; ++i){
    cout << fEval[i] << " ";
  }
  cout << endl;
  
  Eigen::MatrixXcd fMEvecInv = fMEvec.inverse();

  matrixC mult(9, vectorC(9, 0));
  matrixC diag(9, vectorC(9, 0));

  // Multiplying matrix a and b and storing in array mult.
  for(int i = 0; i < 9; ++i){
    for(int j = 0; j < 9; ++j){
      for(int k = 0; k < 9; ++k){
	mult[i][j] += fM[i][k] * fMEvec(k, j);
      }
    }
  }

  for(int i = 0; i < 9; ++i){
    for(int j = 0; j < 9; ++j){
      for(int k = 0; k < 9; ++k){
	diag[i][j] += fMEvecInv(i, k) * mult[k][j];
      }
    }
  }


  for(int i = 0; i < 9; ++i){
    for(int j = 0; j < 9; ++j){
      diag[i][j] *= 1e10;
    }
  }
  
  cout << "TEST D-1 x M x D: MUST BE DIAGONAL \n";
  for(int i = 0; i < 9; ++i){
    cout << diag[i][i] << " ";
  }
  cout << endl;

  // rotate density matrix to vacuum mass basis
  RotateState(true);
  
  ChangeBaseToGM();

  // Convert path length to eV
  double lengthIneV = kKm2eV * p.length;

  for(int i = 0; i < 9; ++i){

    fRt[i] = 0;
    
    for(int k = 0; k < 9; ++k){
      for(int j = 0; j < 9; ++j){
  
       	fRt[i] += exp(fEval[k] * lengthIneV) * fEvec[i][k] * fMEvecInv(k, j) * fR[j];  
      }
    }

    cout << "fRt[" << i << "] = " << fRt[i] << endl;

  }

  ChangeBaseToSU3();

  RotateState(false);
  
  cout << "\n end of propagatepath() \n";
}


//.............................................................................
///
/// Reset the neutrino state back to a pure flavour where it starts
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flv - The neutrino starting flavour.
///
void PMNS_OQS::ResetToFlavour(int flv)
{
  PMNS_Base::ResetToFlavour(flv);

  assert(flv >= 0 && flv < fNumNus);

  for (int i = 0; i < fNumNus; ++i) {
    for (int j = 0; j < fNumNus; ++j) {
      if (i == flv && i == j)
	fRho[i][j] = one;
      else
        fRho[i][j] = zero;
    }
  }

}

//.............................................................................
///
/// Compute oscillation probability of flavour flv from current state
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flv - The neutrino final flavour.
///
/// @return Neutrino oscillation probability
///
double PMNS_OQS::P(int flv)
{
  assert(flv >= 0 && flv < fNumNus);

  return sqrt(fRho[flv][flv].real()*fRho[flv][flv].real() +
	      fRho[flv][flv].imag()*fRho[flv][flv].imag());
}

//.............................................................................
///
/// Set the density matrix from a pure state: redefinition from PMNS_Base.
///
/// @param nu_in - The neutrino initial state in flavour basis.
///
void PMNS_OQS::SetPureState(vectorC nu_in)
{
  assert(nu_in.size() == fNumNus);

  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      fRho[i][j] = conj(nu_in[i]) * nu_in[j];
    }
  }
}


////////////////////////////////////////////////////////////////////////
