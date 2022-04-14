#include "calculateRR.H"
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


fvScalarMatrix YEqn
(
    fvm::ddt(rho, Y)
  + C1*mvConvection->fvmDiv(phi, Y)
  - C2/Le*fvm::laplacian(Y)
 ==
    -omega
  + fvOptions(rho, Y)
);

YEqn.relax();

fvOptions.constrain(YEqn);

YEqn.solve("Y");

fvOptions.correct(Y);
