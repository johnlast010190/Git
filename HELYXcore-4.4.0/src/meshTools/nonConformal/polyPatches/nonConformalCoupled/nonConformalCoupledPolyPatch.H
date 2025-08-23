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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::nonConformalCoupledPolyPatch

Description
    Non-conformal coupled poly patch. As nonConformalPolyPatch, but this patch
    is coupled to another non-conformal patch. Information about the
    owner/neighbour relationship and transformation are made available by this
    class. Also, this patch allows access to the error patch.

See also
    Foam::nonConformalPolyPatch
    Foam::nonConformalErrorPolyPatch

SourceFiles
    nonConformalCoupledPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCoupledPolyPatch_H
#define nonConformalCoupledPolyPatch_H

#include "nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalErrorPolyPatch;

/*---------------------------------------------------------------------------*\
                 Class nonConformalCoupledPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCoupledPolyPatch
:
    public nonConformalPolyPatch
{
    // Private data

        //- Reference to the coupled polyPatch
        const coupledPolyPatch& patch_;

        //- Name of the error patch associated with this cyclic
        mutable word errorPatchName_;

        //- Index of the error patch associated with this cyclic
        mutable label errorPatchID_;


protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames) override;

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex) override;


public:

    //- Runtime type information
    TypeName("nonConformalCoupled");


    // Constructors

        //- Construct from a patch
        nonConformalCoupledPolyPatch(const polyPatch& patch);

        //- Construct from a patch and an original patch name
        nonConformalCoupledPolyPatch
        (
            const polyPatch& patch,
            const word& origPatchName
        );

        //- Construct from a patch and a dictionary
        nonConformalCoupledPolyPatch
        (
            const polyPatch& patch,
            const dictionary& dict
        );

        //- Construct from a patch and a non-conformal coupled patch
        nonConformalCoupledPolyPatch
        (
            const polyPatch& patch,
            const nonConformalCoupledPolyPatch& nccPatch
        );


    //- Destructor
    virtual ~nonConformalCoupledPolyPatch();


    // Member Functions

        //- Does this side own the patch?
        bool owner() const
        {
            return patch_.owner();
        }

        //- Does the coupled side own the patch?
        bool neighbour() const
        {
            return patch_.neighbour();
        }

        //- Return transformation between the coupled patches
        const transformer& transform() const
        {
            return patch_.transform();
        }

        //- Error patch name
        const word& errorPatchName() const
        {
            return patch_.boundaryMesh()[errorPatchID()].name();
        }

        //- Error patchID
        label errorPatchID() const;

        //- Error patch
        const nonConformalErrorPolyPatch& errorPatch() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
