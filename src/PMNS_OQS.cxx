#include <cassert>
#include <complex>
#include <iostream>
#include <math.h>

#include <Eigen/Eigenvalues>

#include "PMNS_OQS.h"

using namespace std;

using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_DensityMatrix::PMNS_DensityMatrix
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_OQS::PMNS_OQS()
    : PMNS_DensityMatrix(), fPhi(), fR(), fRt(), fD(9, vectorC(9, 0)),
      fM(9, vectorC(9, 0)), fMd(9, 9), fHGM(3, vectorC(3, 0)),
      fHeff(3, vectorC(3, 0)), fMEvec(9, 9)
{
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

  fEvalC = vectorC(9, 0);
  fEvec  = matrixC(9, vectorC(9, (0, 0)));
}

// set Heff in vacuum-mass basis
void PMNS_OQS::SetHeff(NuPath p)
{
  double Ve = 0;

  // Set the neutrino path
  SetCurPath(p);

  double rho = fPath.density;
  double zoa = fPath.zoa;

  double kr2GNe =
      kK2 * M_SQRT2 * kGf * rho * zoa; // Electron matter potential in eV

  if (fIsNuBar == 0) { Ve = kr2GNe * 2.; }
  else {
    Ve = -kr2GNe * 2.;
  }

  double s12 = sin(fTheta[0][1]);
  double s13 = sin(fTheta[0][2]);
  double s23 = sin(fTheta[1][2]);

  double c12 = cos(fTheta[0][1]);
  double c13 = cos(fTheta[0][2]);
  double c23 = cos(fTheta[1][2]);

  complexD idelta(0.0, fDelta[0][2]);
  if (fIsNuBar) { idelta = conj(idelta); }
  complexD iphi1(0.0, fPhi[0]);
  complexD iphi2(0.0, fPhi[1]);

  fHeff[0][0] = c12 * c12 * c13 * c13 * Ve;
  fHeff[0][1] = c12 * c13 * c13 * exp(iphi1) * s12 * Ve;
  fHeff[0][2] = c12 * c13 * exp(-idelta + iphi2) * s13 * Ve;

  fHeff[1][0] = c12 * c13 * c13 * exp(-iphi1) * s12 * Ve;
  fHeff[1][1] = c13 * c13 * s12 * s12 * Ve;
  fHeff[1][2] = c13 * exp(-idelta - iphi1 + iphi2) * s12 * s13 * Ve;

  fHeff[2][0] = c12 * c13 * exp(idelta - iphi2) * s13 * Ve;
  fHeff[2][1] = c13 * exp(idelta + iphi1 - iphi2) * s12 * s13 * Ve;
  fHeff[2][2] = s13 * s13 * Ve;

  // add mass terms
  fHeff[1][1] += fDm[1] / (2 * kGeV2eV * fEnergy);
  fHeff[2][2] += fDm[2] / (2 * kGeV2eV * fEnergy);
}

// Set Heff in Gell-Mann basis: set only right-triangle as the left part is
// -right
void PMNS_OQS::SetHGM()
{
  double k = 1. / 2.;

  fHGM[1][2] = k * (real(fHeff[0][0]) - real(fHeff[1][1]));
  fHGM[1][3] = 2. * k * imag(fHeff[0][1]);
  fHGM[1][4] = -k * imag(fHeff[1][2]);
  fHGM[1][5] = -k * real(fHeff[1][2]);
  fHGM[1][6] = -k * imag(fHeff[0][2]);
  fHGM[1][7] = -k * real(fHeff[0][2]);

  fHGM[2][3] = 2. * k * real(fHeff[0][1]);
  fHGM[2][4] = k * real(fHeff[1][2]);
  fHGM[2][5] = -k * imag(fHeff[1][2]);
  fHGM[2][6] = -k * real(fHeff[0][2]);
  fHGM[2][7] = k * imag(fHeff[0][2]);

  fHGM[3][4] = -k * imag(fHeff[0][2]);
  fHGM[3][5] = -k * real(fHeff[0][2]);
  fHGM[3][6] = k * imag(fHeff[1][2]);
  fHGM[3][7] = k * real(fHeff[1][2]);

  fHGM[4][5] = k * (real(fHeff[0][0]) - real(fHeff[2][2]));
  fHGM[4][6] = -k * imag(fHeff[0][1]);
  fHGM[4][7] = k * real(fHeff[0][1]);
  fHGM[4][8] = sqrt(3.) * k * imag(fHeff[0][2]);

  fHGM[5][6] = -k * real(fHeff[0][1]);
  fHGM[5][7] = -k * imag(fHeff[0][1]);
  fHGM[5][8] = sqrt(3.) * k * real(fHeff[0][2]);

  fHGM[6][7] = k * (real(fHeff[1][1]) - real(fHeff[2][2]));
  fHGM[6][8] = sqrt(3.) * k * imag(fHeff[1][2]);

  fHGM[7][8] = sqrt(3.) * k * real(fHeff[1][2]);

  for (int i = 1; i < 9; ++i) {
    for (int j = 1; j < 9; ++j) { fHGM[j][i] = -fHGM[i][j]; }
  }
}

void PMNS_OQS::SetDissipatorElement(int i, int j, double val)
{
  fD[i][j] = val;
  fD[j][i] = conj(fD[i][j]);
}

void PMNS_OQS::SetM()
{
  for (int k = 0; k < 9; ++k) {
    for (int j = 0; j < 9; ++j) { fM[k][j] = fHGM[j][j] + fD[k][j]; }
  }
}

void PMNS_OQS::SetPhi(int i, double val) { fPhi[i - 1] = val; }

template <typename T> void PMNS_OQS::SolveEigenSystem()
{
  cout << "M to diagonalize: \n";

  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      fMd(i, j) = fM[i][j];

      cout << fMd(i, j) << " ";
    }
    cout << "\n";
  }

  Eigen::Ref<T> M(fMd);

  Eigen::ComplexEigenSolver<T> eigensolver(M);

  // Fill fEvalC and fEvec vectors from GSL objects
  for (int i = 0; i < 9; i++) {
    fEvalC[i] = eigensolver.eigenvalues()(i);
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

void PMNS_OQS::ChangeBaseToGM()
{
  cout << "BaseToGM: \n"
       << fRho[0][0] << " " << fRho[0][1] << " " << fRho[0][2] << " "
       << fRho[1][0] << " " << fRho[1][1] << " " << fRho[1][2] << " "
       << fRho[2][0] << " " << fRho[2][1] << " " << fRho[2][2] << endl;

  fR[0] = (fRho[0][0] + fRho[1][1] + fRho[2][2]) / sqrt(6.);
  fR[1] = (fRho[0][1] + fRho[1][0]) / 2.;
  fR[2] = (fRho[0][1] - fRho[1][0]) * complexD(0.0, 1.0) / 2.;
  fR[3] = (fRho[0][0] - fRho[1][1]) / 2.;
  fR[4] = fRho[0][2];
  fR[5] = 0;
  fR[6] = fRho[1][2];
  fR[7] = 0;
  fR[8] = (fRho[0][0] + fRho[1][1] - 2. * fRho[2][2]) / (2. * sqrt(3.));

  cout << "Rmu " << fR[0] << " " << fR[1] << " " << fR[2] << endl;
}

// pensa se fare rho o nuova variabile rho(t)
void PMNS_OQS::ChangeBaseToSU3()
{
  fRho[0][0] = sqrt(2. / 3.) * fRt[0] + fRt[3] + fRt[8] / sqrt(3.);
  fRho[0][1] = fRt[1] - complexD(0.0, 1.0) * fRt[2];
  fRho[0][2] = fRt[4] - complexD(0.0, 1.0) * fRt[5];
  fRho[1][0] = fRt[1] + complexD(0.0, 1.0) * fRt[2];
  fRho[1][1] = sqrt(2. / 3.) * fRt[0] - fRt[3] + fRt[8] / sqrt(3.);
  fRho[1][2] = fRt[6] - complexD(0.0, 1.0) * fRt[7];
  fRho[2][0] = fRt[4] + complexD(0.0, 1.0) * fRt[5];
  fRho[2][1] = fRt[6] + complexD(0.0, 1.0) * fRt[7];
  fRho[2][2] = 1. / sqrt(3.) * (sqrt(2.) * fRt[0] - 2. * fRt[8]);

  cout << "BaseToSU3: \n"
       << fRho[0][0] << " " << fRho[0][1] << " " << fRho[0][2] << "\n"
       << fRho[1][0] << " " << fRho[1][1] << " " << fRho[1][2] << "\n"
       << fRho[2][0] << " " << fRho[2][1] << " " << fRho[2][2] << endl;

  cout << "ABS: " << abs(fRho[0][0]) << " " << abs(fRho[1][1]) << " "
       << abs(fRho[2][2]) << endl;
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
///
/// @param p - A neutrino path segment
///
void PMNS_OQS::PropagatePath(NuPath p)
{
  SetHeff(p);

  SetHGM();

  SetM();

  // Solve for eigensystem
  Diagonalise();

  cout << "MATRIX WITH EIGENVECTORS: \n";

  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      fMEvec(i, j) = fEvec[i][j];

      cout << fEvec[i][j] << " ";
    }
    cout << endl;
  }

  cout << "EIGENVALUES: \n";
  //  for(int i = 0; i < 9; ++i){
  for (int i = 0; i < 9; ++i) { cout << fEvalC[i] << " "; }
  cout << endl;

  Eigen::MatrixXcd fMEvecInv = fMEvec.inverse();

  matrixC mult(9, vectorC(9, 0));
  matrixC diag(9, vectorC(9, 0));

  // Multiplying matrix a and b and storing in array mult.
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      for (int k = 0; k < 9; ++k) { mult[i][j] += fMd(i, k) * fMEvec(k, j); }
    }
  }

  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      for (int k = 0; k < 9; ++k) {
        diag[i][j] += fMEvecInv(i, k) * mult[k][j];
      }
    }
  }

  cout << "TESTARE D-1 x M x D: DEVE VENIRE DIAGONALE \n";
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) { cout << diag[i][j] << " "; }
    cout << endl;
  }

  ChangeBaseToGM();

  // Convert path length to eV
  double lengthIneV = kKm2eV * p.length;

  for (int i = 0; i < 9; ++i) {
    fRt[i] = 0;

    for (int k = 0; k < 9; ++k) {
      for (int j = 0; j < 9; ++j) {
        fRt[i] +=
            exp(fEvalC[k] * lengthIneV) * fEvec[i][k] * fMEvecInv(k, j) * fR[j];
      }
    }

    cout << "fRt[" << i << "] = " << fRt[i] << endl;
  }

  ChangeBaseToSU3();

  cout << "\n fine di propagatepath() \n";
}

////////////////////////////////////////////////////////////////////////
