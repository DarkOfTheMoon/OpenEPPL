/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Harbin Engineering University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::deferred

    Author: Zhang HuiJie <yitianbuji@gmail.com，huijiezhang@hrbeu.edu.cn>
    Institution: Harbin Engineering University

\*---------------------------------------------------------------------------*/

#include "deferred.H"
#include "upwind.H"
#include "convectionScheme.H"
#include "gaussConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::deferred::div
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fv::convectionScheme<Type>> tHOScheme =
        fv::convectionScheme<Type>::New
        (
            mesh,
            faceFlux,
            mesh.divScheme("div("+faceFlux.name()+","+vf.name()+")")
        );
    const fv::convectionScheme<Type>& HOScheme = tHOScheme();

    // TODO 仅考虑了迎风格式，可能还要考虑别的例外情形。
    if (isA<fv::gaussConvectionScheme<Type>>(HOScheme))
    {
        const fv::gaussConvectionScheme<Type>& ghoscheme =
            dynamic_cast<const fv::gaussConvectionScheme<Type>&>(HOScheme);

        if(isA<const upwind<Type>>(ghoscheme.interpScheme()))
        {
            return HOScheme.fvmDiv(faceFlux, vf);
        }
    }

    const upwind<Type> upwindScheme(mesh, faceFlux);
    const surfaceScalarField weights(upwindScheme.weights());

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.lower() = -weights.primitiveField()*faceFlux.primitiveField();
    fvm.upper() = fvm.lower() + faceFlux.primitiveField();
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    GeometricField<Type, fvsPatchField, surfaceMesh> sfCorr
    (
        HOScheme.interpolate(faceFlux, vf)
      - upwindScheme.interpolate(vf)
    );
    fvm += fvc::surfaceIntegrate(faceFlux*sfCorr);

    return tfvm;
}


// ************************************************************************* //
