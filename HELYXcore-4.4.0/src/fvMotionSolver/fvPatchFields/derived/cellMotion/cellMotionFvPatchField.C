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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvPatchFields/derived/cellMotion/cellMotionFvPatchField.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "fields/GeometricFields/pointFields/pointFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const cellMotionFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const cellMotionFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const cellMotionFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::cellMotionFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& p = this->patch();
    const polyPatch& pp = p.patch();
    const fvMesh& mesh = this->internalField().mesh();
    const pointField& points = mesh.points();

    word pfName = this->internalField().name();
    pfName.replace("cell", "point");

    const PointField<Type>& pointMotion =
        this->db().objectRegistry::template
            lookupObject<PointField<Type>>
            (pfName);

    forAll(p, i)
    {
        this->operator[](i) = pp[i].average(points, pointMotion);
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::cellMotionFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}

// ************************************************************************* //
