{
    const scalarField& TCells=T.primitiveFieldRef();
    const scalarField& rhoCells=rho.primitiveFieldRef();
    const scalarField& YCells=Y.primitiveFieldRef();

    forAll (TCells, iCell)
    {
        omega.primitiveFieldRef()[iCell]=beta*beta*(1+q)*(1+q)/2/Le/zeta/zeta*rhoCells[iCell]*rhoCells[iCell]*YCells[iCell]*std::exp(beta*(TCells[iCell]-1)/(1+q*(TCells[iCell]-1)/(1+q)));
    }

    const volScalarField::Boundary& TBf = T.boundaryFieldRef();
    const volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    const volScalarField::Boundary& YBf = Y.boundaryFieldRef();

    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& pY = YBf[patchi];

        forAll(pT, facei)
        {
            omega.boundaryFieldRef()[patchi][facei]=beta*beta*(1+q)*(1+q)/2/Le/zeta/zeta*prho[facei]*prho[facei]*pY[facei]*std::exp(beta*(pT[facei]-1)/(1+q*(pT[facei]-1)/(1+q)));
        }
    }
}
