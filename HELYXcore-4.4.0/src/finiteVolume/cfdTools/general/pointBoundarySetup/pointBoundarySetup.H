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
    (c) 2010-2016 Engys Ltd.

Class
    Foam::pointBoundarySetup

Description
    Interpret a structured input dictionary to create a boundary for
    a pointField. Note this is not the standard boundaryField format dictionary.

SourceFiles
    pointBoundarySetup.C

\*---------------------------------------------------------------------------*/

#ifndef pointBoundarySetup_H
#define pointBoundarySetup_H

#include "meshes/pointMesh/pointMesh.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pointBoundarySetup Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointBoundarySetup
:
    public PtrList<pointPatchField<Type>>
{

    // boundaryField dictionary for CBC usage
    autoPtr<dictionary> dictPtr_;

    //- default fvPatch-type wide boundary dictionary entries
    void defaultBoundaries
        (const pointMesh&);

    //- default PointPatch-type wide boundaries
    void defaultBoundaries
        (const pointMesh&, const DimensionedField<Type, pointMesh>&);

    //- fvPatch-type wide boundary dictionary entries
    void typeBoundaries
    (
        const pointMesh&,
        const dictionary&
    );

    //- PointPatch-type wide boundaries
    void typeBoundaries
    (
        const pointMesh&,
        const DimensionedField<Type, pointMesh>&,
        const dictionary&
    );

    //- boundary dictionary entries for named patches
    void namedBoundaries
    (
        const pointMesh&,
        const dictionary&,
        bool strictPatchNameChecking
    );

    //- boundaries for named patches
    void namedBoundaries
    (
        const pointMesh&,
        const DimensionedField<Type, pointMesh>&,
        const dictionary&,
        bool strictPatchNameChecking
    );

public:

    // Constructors

        //- Construct dictionary
        pointBoundarySetup
        (
            const pointMesh&,
            const word&,
            const dictionary&,
            bool strictPatchNameChecking = true
        );

        //- Construct from components.
        pointBoundarySetup
        (
            const pointMesh&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&,
            bool strictPatchNameChecking = true
        );

    // Destructor

        ~pointBoundarySetup();

    // Member functions

    // Return dictionary
    const dictionary& boundaryDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "cfdTools/general/pointBoundarySetup/pointBoundarySetup.C"
#endif

#endif

// ************************************************************************* //
