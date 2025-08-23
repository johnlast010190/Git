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
    (c) 2011-2023 OpenFOAM Foundation

Class
    Foam::solidMulticomponentThermo

Description
    Foam::solidMulticomponentThermo

SourceFiles
    solidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidMulticomponentThermo_H
#define solidMulticomponentThermo_H

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "solidThermo/solidThermo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/


class solidMulticomponentThermo
:
    virtual public solidThermo
{

public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;

    //- Runtime type information
    TypeName("solidMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidMulticomponentThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );

     // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidMulticomponentThermo,
        dictionary,
        (const objectRegistry& obr, const dictionary& dict, const word& phaseName),
        (obr, dict, phaseName)
    );


    // Selectors

        //- Standard selection based onobjectRegistry
        static autoPtr<solidMulticomponentThermo> New
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

         //- Standard selection based on objectRegistry amd dictionary
        static autoPtr<solidMulticomponentThermo> New
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidMulticomponentThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;

};


/*---------------------------------------------------------------------------*\
              Class solidMulticomponentThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class solidMulticomponentThermo::implementation
:
    virtual public solidMulticomponentThermo
{

public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry&, const word& phaseName);

        //- Construct from mesh and phase name
        implementation
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
                 Class solidMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public solidThermo::implementation,
    public solidMulticomponentThermo::implementation
{
public:

    // Constructors

        //- Construct from object registry and phase name
        composite
        (
            const objectRegistry& obr,
            const word& phaseName
        )
        :
            basicThermo::implementation(obr, phaseName),
            solidThermo::implementation(obr, phaseName),
            solidMulticomponentThermo::implementation(obr, phaseName)
        {}

        //- Construct from object registry and phase name
        composite
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName = word::null
        )
        :
            basicThermo::implementation(obr, phaseName),
            solidThermo::implementation(obr, dict, phaseName),
            solidMulticomponentThermo::implementation(obr, dict, phaseName)
        {}


    // Member Functions

        //- The dpdt term should not be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
