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
    (c) 2011-2022 OpenFOAM Foundation

Class
    Foam::rhoMulticomponentThermo

Group
    grpRhoThermo

Description
    Base-class for multi-component fluid thermodynamic properties based on
    density.

See also
    Foam::basicThermo

SourceFiles
    rhoMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoMulticomponentThermo_H
#define rhoMulticomponentThermo_H

#include "rhoThermo/rhoThermo.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rhoMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoMulticomponentThermo
:
    virtual public rhoThermo,
    virtual public fluidMulticomponentThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("rhoMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoMulticomponentThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );


    // Selectors

        //- Standard selection based onobjectRegistry
        static autoPtr<rhoMulticomponentThermo> New
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~rhoMulticomponentThermo();
};


/*---------------------------------------------------------------------------*\
              Class rhoMulticomponentThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class rhoMulticomponentThermo::implementation
:
    virtual public rhoMulticomponentThermo
{

public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
                 Class rhoMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public rhoThermo::implementation,
    public fluidMulticomponentThermo::implementation,
    public rhoMulticomponentThermo::implementation
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
            fluidThermo::implementation(obr, phaseName),
            rhoThermo::implementation(obr, phaseName),
            fluidMulticomponentThermo::implementation(obr, phaseName),
            rhoMulticomponentThermo::implementation(obr, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
