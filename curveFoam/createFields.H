#include "createRDeltaT.H"

#include "readAditionalParameters.H"

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimDensity
);

volScalarField rhoEOS
(
    IOobject
    (
        "rhoEOS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimDensity
);
volScalarField Hby
(
    IOobject
    (
        "Hby",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(1,-1,-2,0,0)
);

#include "updaeRhoByEOS.H"
rho=rhoEOS;

volScalarField omega
(
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1,-3,-1,0,0)
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

fields.add(Y);
fields.add(T);

#include "createMRF.H"
#include "createFvOptions.H"
