#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_Fast {
    public:

    PMNS_OQS();          ///< Constructor
    virtual ~PMNS_OQS(); ///< Destructor
    
    /// Compute the probability matrix
    virtual void SetPhi(int i, double val);
    virtual void SetDissipatorElement(int i, int j, double val);
    
    protected:
    
    virtual void InitializeVectors();
    virtual void SetHeff(NuPath p);
    virtual void SetHGM();
    virtual void SetM();
    virtual void Diagonalise();
    virtual void ChangeBaseToGM();
    virtual void ChangeBaseToSU3();
    virtual void RotateState();
    
    /// Specialized solver for NxN matrices
    template <typename T> void SolveEigenSystem();
    
    // Resetting and propagating
    virtual void ResetToFlavour(int flv); ///< Reset neutrino state to pure flavour flv
    virtual void SetPureState(vectorC nu_in); ///< Set the density matrix from a pure state
    virtual void PropagatePath(NuPath p); ///< Propagate neutrino through a single path
    virtual double P(int flv); ///< Return the probability of final state in flavour flv

    double fPhi[2];  ///< Majorana phases
    complexD fR[9];
    complexD fRt[9];

    matrixC fRho; ///< The neutrino density matrix state
    matrixC fHeff;
    matrixC fHGM; ///< Heff in Gell-Mann basis
    matrixC fD;  ///< Off-diagonal, 9x9 dissipator
    matrixC fM;  ///< M

    Eigen::MatrixXcd fMd;
    Eigen::MatrixXcd fMd8;
    Eigen::MatrixXcd fMEvec;

    vectorC fEval; ///< Eigenvalues of the Hamiltonian
    matrixC fEvec; ///< Eigenvectors of the Hamiltonian
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
