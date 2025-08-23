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
    (c) 2014-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/limitedGradSchemes/cellPartialLimitedGrad/cellPartialLimitedGrad.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "meshQualityMetrics/meshQualityMetrics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::cellPartialLimitedGrad<Type>::blendLimiter
(
    Field<Type>& limiter
) const
{
    const meshQualityMetrics& mqMetrics =
        meshQualityMetrics::New(this->mesh(), meshObjectName_);
    const boolList& markedCells = mqMetrics.badQualityCells();
    forAll(limiter, cI)
    {
        if (!markedCells[cI]) limiter[cI] = pTraits<Type>::one;
    }

    if (this->mesh().time().outputTime() && mqMetrics.writeField())
    {
        this->writeField(markedCells);
    }
}


template<class Type>
void Foam::fv::cellPartialLimitedGrad<Type>::writeField
(
    const boolList& bqMarker
) const
{
    volScalarField outField
    (
        IOobject
        (
            meshObjectName_ + "_badQualityCells",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0),
        "zeroGradient"
    );
    forAll(bqMarker, cI)
    {
        if (bqMarker[cI]) outField[cI] = 1;
    }
    outField.write();
}


// ************************************************************************* //
