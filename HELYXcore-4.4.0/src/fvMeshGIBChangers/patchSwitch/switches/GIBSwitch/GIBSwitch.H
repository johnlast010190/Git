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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::GIBSwitch

Description
    Abstract base class for switching GIB on and off.

SourceFiles
    GIBSwitch.C
    GIBSwitchNew.C

\*---------------------------------------------------------------------------*/

#ifndef GIBSwitch_H
#define GIBSwitch_H

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                          Class GIBSwitch Declaration
\*---------------------------------------------------------------------------*/

class GIBSwitch
:
    public IOdictionary
{
    // Private Member Functions

        //- Initialise
        void initialise();

        //- Find the master and slave IDs
        void findPatchIDs();

        //- Create the IOobject for writing the coeffs
        IOobject createIOobject(const fvMesh&,  const dictionary&) const;

        //- Condition for enabling the GIB
        virtual bool enableCondition() = 0;

        //- Condition for disabling the GIB
        virtual bool disableCondition() = 0;

        //- Enable the GIB fvPatch
        void enableBoundary();

        //- Disable the GIB fvPatch
        void disableBoundary();


protected:

        //- Reference to the fvMesh
        const fvMesh& mesh_;

        //- Reference to subDict of dynamicFvMesh
        const dictionary& dict_;

        //- zoneID of gib patch
        label zoneID_;

        //- Master ID
        label masterID_;

        //- Slave ID
        label slaveID_;

        //- initCondition
        Switch startActive_;

        //- When to start checking
        const scalar startTime_;

        //- Timesteps needed to enable or disable
        label condInterval_;

        //- n timesteps (needed for restarts etc)
        label curCondInterval_;


public:

    //- Runtime type information
    TypeName("GIBSwitch");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            GIBSwitch,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct
        explicit GIBSwitch(const fvMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        GIBSwitch(const GIBSwitch&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<GIBSwitch> New(const fvMesh&, const dictionary&);


    //- Destructor
    virtual ~GIBSwitch(){};


    // Member Functions

        //- Check condition and update patch status if needed
        virtual bool update();

        //- Update from dictionary
        virtual void updateDict();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GIBSwitch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
