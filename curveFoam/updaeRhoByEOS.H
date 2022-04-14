{
    const scalarField& TCells=T.primitiveFieldRef();
    const scalarField& pCells=p.primitiveFieldRef();

    forAll (TCells, iCell)
    {
        rhoEOS.primitiveFieldRef()[iCell]=(1+Gamma*Lambda*pCells[iCell])/(1+q*TCells[iCell]);
    }

    const volScalarField::Boundary& TBf = T.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p.boundaryFieldRef();

    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pT = TBf[patchi];

        forAll(pT, facei)
        {
            rhoEOS.boundaryFieldRef()[patchi][facei]=(1+Gamma*Lambda*pp[facei])/(1+q*pT[facei]);
        }
    }
}
