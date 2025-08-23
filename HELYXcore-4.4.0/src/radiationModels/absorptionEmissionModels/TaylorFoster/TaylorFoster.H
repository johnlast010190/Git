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
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::radiationModels::TaylorFoster

Description
    TaylorFoster radiation absorption and emission coefficients for continuous
    phase. The gas absorption coeff. is increased by a fraction related to
    smoke/soot according to:
    ksoot = b1*rhom[1+bT(T-2000)]
    with rhom soot mass concentration in [kg/m3] and b1 and bT constants.

SourceFiles
    TaylorFoster.C

\*---------------------------------------------------------------------------*/

#ifndef radiationTaylorFoster_H
#define radiationTaylorFoster_H

#include "absorptionEmissionModels/constantAbsorptionEmission/constantAbsorptionEmission.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                  Class constantAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class TaylorFoster
:
    public constantAbsorptionEmission
{
private:
    word rhoName_;

    word sootName_;

    word TName_;

public:

    //- Runtime type information
    TypeName("TaylorFoster");


    // Constructors

    //- Construct from components
    TaylorFoster
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const word& typeNameDerived = typeName
    );


    //- Destructor
    virtual ~TaylorFoster();


    // Member Functions

    //- Absorption coefficient for continuous phase
    virtual tmp<volScalarField> aCont(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
