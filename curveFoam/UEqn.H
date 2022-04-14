// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);


U=C/Pr*fvc::grad(p);

fvOptions.correct(U);

