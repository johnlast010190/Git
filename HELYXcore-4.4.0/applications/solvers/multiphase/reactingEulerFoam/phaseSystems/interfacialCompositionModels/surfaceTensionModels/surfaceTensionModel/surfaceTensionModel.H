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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::surfaceTensionModel

Description

SourceFiles
    surfaceTensionModel.C
    newAspectRatioModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionModel_H
#define surfaceTensionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                           Class surfaceTensionModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("surfaceTensionModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceTensionModel,
        legacyMultiphase,
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        ),
        (dict, pair, registerObject)
    );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimSigma;


    // Constructors

        //- Construct from a dictionary and a phase pair
        surfaceTensionModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~surfaceTensionModel();


    // Selectors

        static autoPtr<surfaceTensionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> sigma() const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
