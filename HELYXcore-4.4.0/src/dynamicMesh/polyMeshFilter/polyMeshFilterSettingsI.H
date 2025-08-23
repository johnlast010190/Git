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
    (c) 2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dictionary&
Foam::polyMeshFilterSettings::collapseEdgesCoeffDict() const
{
    return collapseEdgesCoeffDict_;
}

inline const Foam::dictionary&
Foam::polyMeshFilterSettings::collapseFacesCoeffDict() const
{
    return collapseFacesCoeffDict_;
}

inline const Foam::dictionary&
Foam::polyMeshFilterSettings::meshQualityCoeffDict() const
{
    return meshQualityCoeffDict_;
}


inline const Foam::Switch&
Foam::polyMeshFilterSettings::controlMeshQuality() const
{
    return controlMeshQuality_;
}

inline const Foam::scalar& Foam::polyMeshFilterSettings::minLen() const
{
    return minLen_;
}

inline const Foam::scalar& Foam::polyMeshFilterSettings::maxCos() const
{
    return maxCos_;
}

inline const Foam::scalar&
Foam::polyMeshFilterSettings::edgeReductionFactor() const
{
    return edgeReductionFactor_;
}

inline const Foam::label& Foam::polyMeshFilterSettings::maxIterations() const
{
    return maxIterations_;
}

inline const Foam::label& Foam::polyMeshFilterSettings::maxSmoothIters() const
{
    return maxSmoothIters_;
}

inline const Foam::scalar&
Foam::polyMeshFilterSettings::initialFaceLengthFactor() const
{
    return initialFaceLengthFactor_;
}

inline const Foam::scalar&
Foam::polyMeshFilterSettings::faceReductionFactor() const
{
    return faceReductionFactor_;
}

inline const Foam::label&
Foam::polyMeshFilterSettings::maxPointErrorCount() const
{
    return maxPointErrorCount_;
}


// ************************************************************************* //
