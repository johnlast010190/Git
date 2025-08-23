/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2017 Engys Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/makeFvOption.H"
#include "SemiImplicitSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvOption(SemiImplicitSource, scalar);
makeFvOption(SemiImplicitSource, vector);
makeFvOption(SemiImplicitSource, sphericalTensor);
makeFvOption(SemiImplicitSource, symmTensor);
makeFvOption(SemiImplicitSource, tensor);


#define addFvOptionAndRename(Option, Name, Type)                              \
                                                                              \
    Foam::fv::option::adddictionaryConstructorToTable                         \
        <Foam::fv::Option<Foam::Type>>                                        \
        add_##Option##Type##Name(#Type#Name)

addFvOptionAndRename(SemiImplicitSource, TimeDependentSemiImplicitSource, scalar);
addFvOptionAndRename(SemiImplicitSource, TimeDependentSemiImplicitSource, vector);
addFvOptionAndRename(SemiImplicitSource, TimeDependentSemiImplicitSource, sphericalTensor);
addFvOptionAndRename(SemiImplicitSource, TimeDependentSemiImplicitSource, symmTensor);
addFvOptionAndRename(SemiImplicitSource, TimeDependentSemiImplicitSource, tensor);


// ************************************************************************* //
