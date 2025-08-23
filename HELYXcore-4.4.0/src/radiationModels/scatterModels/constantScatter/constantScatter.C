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
    (c) 2015 Engys Ltd.
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "scatterModels/constantScatter/constantScatter.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace scatterModels
{
    defineTypeNameAndDebug(constantScatter, 0);
    addToRunTimeSelectionTable
    (
        scatterModel,
        constantScatter,
        dictionary
    );
    // Backward compatibility
    addNamedToRunTimeSelectionTable
    (
        scatterModel,
        constantScatter,
        dictionary,
        constantScatter
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::scatterModels::constantScatter::constantScatter
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& typeNameDerived
)
:
    scatterModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeNameDerived + "Coeffs")),
    sigma_("sigma", dimless/dimLength, coeffsDict_),
    C_("C", dimless, coeffsDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::scatterModels::constantScatter::~constantScatter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::scatterModels::constantScatter::sigmaEff() const
{
    return volScalarField::New("sigma", mesh_, sigma_*(3.0 - C_));
}


// ************************************************************************* //
