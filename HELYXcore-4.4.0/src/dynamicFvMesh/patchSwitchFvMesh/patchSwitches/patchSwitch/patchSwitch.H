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
    Foam::patchSwitch

Description
    Abstract base class for switching GIB on and off.

SourceFiles
    patchSwitch.C
    patchSwitchNew.C

\*---------------------------------------------------------------------------*/

#ifndef patchSwitch_H
#define patchSwitch_H

#include "db/dictionary/dictionary.H"
#include "fvMesh/fvMesh.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchSwitch Declaration
\*---------------------------------------------------------------------------*/

class patchSwitch
:
    public IOdictionary
{

    // Private Member

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

        //- reference to fvMesh
        const fvMesh& mesh_;

        //- reference to subDict of dynamicFvMesh
        const dictionary& dict_;

        //- zoneID of gib patch
        label zoneID_;

        //- master ID
        label masterID_;

        //- slave ID
        label slaveID_;

        //- initCondition
        Switch startActive_;

        //- When to start checking
        const scalar startTime_;

        //- timesteps that needed to enable or disable
        label condInterval_;

        //- n timesteps (needed for restarts etc)
        label curCondInterval_;


public:

    //- Runtime type information
    TypeName("patchSwitch");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        patchSwitch,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );


    // Constructors

        //- Construct
        explicit patchSwitch
        (
            const fvMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        patchSwitch(const patchSwitch&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<patchSwitch> New
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~patchSwitch(){};


    // Member Functions

        //- Check condition and update patch status if needed
        virtual bool update();

    // IO

        //- update from dictionary
        virtual void updateDict();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchSwitch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
