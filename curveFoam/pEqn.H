#include "calculateRR.H"
/*
tmp<volScalarField> Hby
(
    volScalarField::New
    (
        "Hby",
        p.mesh(),
        dimensionSet(1,-1,-2,0,0)
    )
);
*/
dimensionedScalar tmp1(dimensionSet(0,2,-2,-1,0),1.0);
dimensionedScalar tmp2(dimensionSet(0,2,-2,0,0),1.0);
dimensionedScalar tmp5(dimensionSet(0,-2,2,0,0),1.0);

const scalarField& pCells=p.primitiveFieldRef();

forAll (pCells, iCell)
{
    Hby.ref()[iCell]=1+Gamma*Lambda*pCells[iCell];
}
const volScalarField::Boundary& pBf = p.boundaryFieldRef();
forAll (pBf, patchi)
{
    forAll(pBf[patchi], facei)
    {
        Hby.boundaryFieldRef()[patchi][facei]=1+Gamma*Lambda*pBf[patchi][facei];
    } 
}

surfaceScalarField rhof("rhof", fvc::interpolate(Hby));
/*
fvScalarMatrix pEqn
(
    Lambda*fvm::ddt(p)
  + C/Pr*C1*fvm::laplacian(rhof, p)
  ==q*(C2*fvc::laplacian(T)*tmp1+omega*tmp2)
);
*/
fvScalarMatrix pEqn
(
    Lambda*fvm::ddt(p)
  + C/Pr*C1*fvm::laplacian(rhof, p)
  ==q*(C2*fvc::laplacian(T)*tmp1+omega*tmp2)
);

pEqn.solve();

U = C/Pr*fvc::grad(p);
//phi = -pEqn.flux()*tmp5;
#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"


phi = fvc::interpolate(rho*U) & mesh.Sf();

U.correctBoundaryConditions();
fvOptions.correct(U);

 
dpdt = fvc::ddt(p);
