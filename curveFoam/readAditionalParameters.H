Info<<"Reading parameters\n"<<endl;

IOdictionary additionalParaters
(
    IOobject
    (
        "additionalParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
/*
dimensionedScalar qulv
(
    "qulv",
    dimensionSet(0,-1,0,0,0),
    1
);
qulv=dimensionedScalar(additionalParaters.lookup("qulv"));
*/
const scalar qulv(additionalParaters.lookupOrDefault<scalar>("qulv", 100));
const scalar Le(additionalParaters.lookupOrDefault<scalar>("Le", 1));
const scalar Pr(additionalParaters.lookupOrDefault<scalar>("Pr", 0.707));
const scalar Lambda(additionalParaters.lookupOrDefault<scalar>("Lambda", 0.25));
const scalar Gamma(additionalParaters.lookupOrDefault<scalar>("Gamma",1.4));
const scalar q(additionalParaters.lookupOrDefault<scalar>("q",5));
const scalar beta(additionalParaters.lookupOrDefault<scalar>("beta", 10));
const scalar zeta(additionalParaters.lookupOrDefault<scalar>("zeta",1.0548));



dimensionedScalar C("C",dimensionSet(-1,3,1,0,0),(4*(qulv+1)*(qulv+1)*std::log(qulv+1)*std::log(qulv+1)-qulv*qulv*(qulv+2)*(qulv+2))/(8*qulv*qulv*qulv*qulv*(qulv+2)));
const scalar C1=2/(qulv+2);
dimensionedScalar C2("C2",dimensionSet(1,-1,-1,0,0),2*std::log(qulv+1)/qulv/(qulv+2));
dimensionedScalar C3("C3",dimensionSet(0,-2,2,1,0),Lambda*(Gamma-1)/q);

Info<< "C="<<C<<endl;
Info<< "C1="<<C1<<endl;
Info<< "C2="<<C2<<endl;
Info<< "C3="<<C3<<endl;


/*
dimensionedScalar C("C",dimensionSet(-1,3,1,0,0),-1/12.0);
scalar C1=1.0;
dimensionedScalar C2("C2",dimensionSet(1,-1,-1,0,0),1.0);
dimensionedScalar C3("C3",dimensionSet(0,-2,2,1,0),Lambda*(Gamma-1)/q);
*/
