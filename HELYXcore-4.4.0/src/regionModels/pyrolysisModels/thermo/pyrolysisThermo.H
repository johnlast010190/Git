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
    (c) 2016 OpenCFD Ltd
    (c) 2024-2025 Engys Ltd.

Class
    Foam::regionModels::pyrolysisThermo

Description
    Pyrolysis model which solves only the energy equation in the region.

SourceFiles
    pyrolysisThermo.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisThermo_H
#define pyrolysisThermo_H

#include "pyrolysisModels/pyrolysisModel/pyrolysisModel.H"
#include "fields/volFields/volFieldsFwd.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "solidThermo/solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                     Class pyrolysisThermo Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisThermo
:
    public pyrolysisModel
{
public:

    //- Runtime type information
    TypeName("thermo");


    // Constructors

        //- Construct from type name and mesh
        pyrolysisThermo
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name and mesh and dict
        pyrolysisThermo
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );

        //- Disallow default bitwise copy construction
        pyrolysisThermo(const pyrolysisThermo&) = delete;


    //- Destructor
    virtual ~pyrolysisThermo();


    // Member Functions

        //- Evolve the pyrolysis equations
        virtual void evolveRegion() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pyrolysisThermo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
