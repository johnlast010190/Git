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
    (c) 2020-2022 OpenFOAM Foundation

Class
    Foam::fluidMulticomponentThermo

Description
    Base-class for multi-component fluid thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    fluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMulticomponentThermo_H
#define fluidMulticomponentThermo_H

#include "fluidThermo/fluidThermo.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("fluidMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidMulticomponentThermo,
        objectRegistry,
        (const objectRegistry& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<fluidMulticomponentThermo> New
        (
            const objectRegistry& obr,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~fluidMulticomponentThermo();


    // Member Functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;
};


/*---------------------------------------------------------------------------*\
                     Class fluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentThermo::implementation
:
    virtual public fluidMulticomponentThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
