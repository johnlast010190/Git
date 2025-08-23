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
    (c) 2017 OpenFOAM Foundation

Class
    Foam::surfaceTensionModels::liquidPropertiesSurfaceTension

Description
    Temperature-dependent surface tension model in which the surface tension
    function provided by the phase Foam::liquidPropertiesSurfaceTension class is used.

Usage
    \table
        Property     | Description               | Required    | Default value
        phase        | Phase name                | yes         |
    \endtable

    Example of the surface tension specification:
    \verbatim
        sigma
        {
            type    liquidProperties;
            phase   water;
        }
    \endverbatim

See also
    Foam::surfaceTensionModel

SourceFiles
    liquidPropertiesSurfaceTension.C

\*---------------------------------------------------------------------------*/

#ifndef liquidPropertiesSurfaceTension_H
#define liquidPropertiesSurfaceTension_H

#include "surfaceTensionModels/surfaceTensionModel/surfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace surfaceTensionModels
{

/*---------------------------------------------------------------------------*\
              Class liquidPropertiesSurfaceTension Declaration
\*---------------------------------------------------------------------------*/

class liquidPropertiesSurfaceTension
:
    public surfaceTensionModel
{
    // Private data

        //- Name of the liquid phase
        word phaseName_;


public:

    //- Runtime type information
    TypeName("liquidProperties");


    // Constructors

        //- Construct from dictionary and mesh
        liquidPropertiesSurfaceTension
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~liquidPropertiesSurfaceTension();


    // Member Functions

        //- Surface tension coefficient
        virtual tmp<volScalarField> sigma() const;

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict);

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
