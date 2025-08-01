#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_DensityMatrix.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_DensityMatrix {
    public:
      PMNS_OQS();          ///< Constructor
      virtual ~PMNS_OQS(); ///< Destructor

      virtual void SetParameterisation(int par);
      virtual void SetPhi(int i, double val);
      virtual void SetDecoElement(int i, double val);
      virtual void SetDecoAngle(int i, int j, double th);
      virtual void SetPower(int n);
      virtual void SetIsNuBar(bool isNuBar);

      virtual matrixD ProbMatrix(int nflvi, int nflvf);

    protected:
      virtual void SetDissipator();
      virtual void SetHeff(NuPath p);
      virtual void SetHGM();
      virtual void SetM();
      virtual void BuildHms();
      virtual void BuildUM();
      virtual void RotateState(bool to_mass); ///< Rotate rho to/from mass basis
      virtual void ChangeBaseToGM();
      virtual void ChangeBaseToSU3();

      virtual void Propagate();
      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p);

      int     fPower;
      int     fParameterisation;
      double  fPhi[2]; ///< Majorana phases
      vectorD fR;
      vectorD fRt;

      matrixC fHeff;
      matrixD fHGM;

      matrixD fD;   ///< Off-diagonal, 9x9 dissipator
      matrixC fUM;  ///< PMNS Matrix
      vectorD fa;   ///< a vector
      matrixD fcos; ///< cosines ai . aj

      Eigen::MatrixXd fM; ///< Buffer matrix for exponential

      bool fBuiltDissipator; ///< Flag to rebuilt D
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
