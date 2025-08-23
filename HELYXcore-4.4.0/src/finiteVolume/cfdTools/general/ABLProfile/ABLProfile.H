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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::ABLProfile

Description
    Abstract base class ABL profiles

SourceFiles
    ABLProfile.C

\*---------------------------------------------------------------------------*/

#ifndef ABLProfile_H
#define ABLProfile_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "cfdTools/general/patchDistanceFunction/patchDistanceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ABLProfile Declaration
\*---------------------------------------------------------------------------*/

class ABLProfile
{
protected:

    // Protected Data

        //- Patch
        const fvPatch& patch_;

        //- Distance field
        autoPtr<patchDistanceFunction> Z_;


public:

    //- Runtime type information
    TypeName("ABLProfile");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ABLProfile,
        dictionary,
        (
            const fvPatch& patch,
            const dictionary& ABLProfileDict
        ),
        (patch, ABLProfileDict)
    );


    // Constructors

        //- Null constructor
        ABLProfile(const fvPatch&);

        //- Construct from components
        ABLProfile(const fvPatch&, const dictionary&);

        //- copy constructor
        ABLProfile(const ABLProfile&);

        //- Construct and return a clone
        virtual ABLProfile* clone() const = 0;


    // Selectors

        //- Return an autoPtr to the selected ABLProfile
        static autoPtr<ABLProfile> New
        (
            const fvPatch& patch,
            const dictionary& ABLProfileDict
        );


    //- Destructor
    virtual ~ABLProfile()
    {}


    // Member Functions

        //- Return distance field
        inline const Field<scalar> Z() const
        {
            return Z_->value();
        }


        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            scalarField& L,
            symmTensorField& R
        ) = 0;

            //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            vectorField& L,
            symmTensorField& R
        ) = 0;

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ABLProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
