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
    (c) 2011 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd.

Class
    Foam::catchRatio

Group
    grpReportFunctionObjects

Description
    Computes local and global catch ratio of disperse Eulerain phases.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    catchRatio.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catchRatio_H
#define functionObjects_catchRatio_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class catchRatio Declaration
\*---------------------------------------------------------------------------*/

class catchRatio
:
    public fvMeshFunctionObject
{
public:

    // Public enumerations

        //- Name list of phases
        wordList phaseNames_;

        //- Scalar list of phases size distribution
        scalarList psd_;

        //- Global rain fall intensity
        scalar Rh_;


    // Private Member Functions

        //- Return phase alpha field
        const volVectorField& Ud(const label idx) const;

        //- Return phase velocity field
        const volScalarField& alphad(const label idx) const;

        //- Create and store catch ratio fields
        void createAndStoreFields();


public:

    //- Runtime type information
    TypeName("catchRatio");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        catchRatio
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        catchRatio(const catchRatio&) = delete;


    //- Destructor
    virtual ~catchRatio();


    // Member Functions

        //- Read the catchRatio data
        virtual bool read(const dictionary&);

        //- Execute, calculate fluid height and write to screen
        virtual bool execute();

        //- Does nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const catchRatio&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
