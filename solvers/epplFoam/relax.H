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

    Author: Zhang HuiJie <yitianbuji@gmail.com，huijiezhang@hrbeu.edu.cn>
    Institution: Harbin Engineering University

\*---------------------------------------------------------------------------*/

#ifndef relax_H
#define relax_H

#include "fvMatrix.H"

// ************************************************************************* //

namespace Foam
{
template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > relax(const fvMatrix<Type>& m)
{
    tmp<fvMatrix<Type> > tmatrix(new fvMatrix<Type>(m));
    tmatrix.ref().relax();

    return tmatrix;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > relax(tmp<fvMatrix<Type>> tm)
{
    tm.ref().relax();
    return tm;
}

}// End of namespace Foam

#endif

// ************************************************************************* //
