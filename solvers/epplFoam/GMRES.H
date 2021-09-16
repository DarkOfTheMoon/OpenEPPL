/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::GMRES

Description
    Preconditioned Generalised Minimal Residual solver with
    run-time selectable preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    GMRES.C

\*---------------------------------------------------------------------------*/

#ifndef GMRES_H
#define GMRES_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class GMRES Declaration
\*---------------------------------------------------------------------------*/

class GMRES
:
    public lduMatrix::solver
{
    // Private Data

        //- Krylov space dimension
        label nDirs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        GMRES(const GMRES&);

        //- Disallow default bitwise assignment
        void operator=(const GMRES&);


        //- Givens rotation
        void givensRotation
        (
            const scalar& H,
            const scalar& beta,
            scalar& c,
            scalar& s
        ) const;


public:

    //- Runtime type information
    TypeName("GMRES");


    // Constructors

        //- Construct from matrix components and solver data stream
        GMRES
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~GMRES()
    {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
