#include "calculateRR.H"
dimensionedScalar tmp3(dimensionSet(0,0,0,1,0),1.0);
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho,T)
        +C1*mvConvection->fvmDiv(phi,T)
        -C2*fvm::laplacian(T)
        -C3*dpdt
        ==
        omega*tmp3
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

#include "updaeRhoByEOS.H"
    
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
