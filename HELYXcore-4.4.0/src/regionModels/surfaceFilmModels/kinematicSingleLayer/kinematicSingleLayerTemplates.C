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
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/kinematicSingleLayer/kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
tmp<Type> kinematicSingleLayer::constrainFilmField
(
    const tmp<Type>& tfield,
    const typename Type::cmptType& value
)
{
    tmp<Type> tresult(tfield);
    Type& result = tresult.ref();

    typename Type::Boundary& fieldBf = result.boundaryFieldRef();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        fieldBf[patchi] = value;

        DebugInFunction
            << "Constraining " << tfield().name()
            << " boundary " << tfield().boundaryField()[patchi].patch().name()
            << " to " << value << endl;
    }

    forAll(passivePatchIDs(), i)
    {
        const label patchi = passivePatchIDs()[i];
        fieldBf[patchi] = value;
        DebugInFunction
            << "Constraining " << tfield().name()
            << " boundary " << tfield().boundaryField()[patchi].patch().name()
            << " to " << value << endl;
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
