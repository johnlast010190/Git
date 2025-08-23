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
    (c) 2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::nonConformalPolyFacesFvsPatchLabelField

SourceFiles
    nonConformalPolyFacesFvsPatchLabelField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalPolyFacesFvsPatchLabelField_H
#define nonConformalPolyFacesFvsPatchLabelField_H

#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nonConformalPolyFacesFvsPatchLabelField Declaration
\*---------------------------------------------------------------------------*/

class nonConformalPolyFacesFvsPatchLabelField
:
    public fvsPatchLabelField
{
public:

    //- Runtime type information
    TypeName("nonConformalPolyFaces");


    // Constructors

        //- Construct from patch and internal field
        nonConformalPolyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalPolyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalPolyFacesFvsPatchLabelField
        (
            const nonConformalPolyFacesFvsPatchLabelField&,
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalPolyFacesFvsPatchLabelField
        (
            const nonConformalPolyFacesFvsPatchLabelField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalPolyFacesFvsPatchLabelField
        (
            const nonConformalPolyFacesFvsPatchLabelField&,
            const DimensionedField<label, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchLabelField> clone
        (
            const DimensionedField<label, surfaceMesh>& iF
        ) const override
        {
            return tmp<fvsPatchLabelField>
            (
                new nonConformalPolyFacesFvsPatchLabelField(*this, iF)
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
