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
    (c) 2011-2018 OpenFOAM Foundation
    (c) 2016, 2020 OpenCFD Ltd.
    (c) 2017-2024 Engys

\*---------------------------------------------------------------------------*/


inline bool Foam::AMIInterpolation::upToDate() const
{
    return upToDate_;
}


inline bool& Foam::AMIInterpolation::upToDate()
{
    return upToDate_;
}


inline const Foam::labelList&
Foam::AMIInterpolation::overlapTargetFaces() const
{
    return (*overlapTgtAddrPtr_);
}

inline const Foam::labelList&
Foam::AMIInterpolation::nonOverlapTargetFaces() const
{
    return (*nonOverlapTgtAddrPtr_);
}

inline const Foam::labelList&
Foam::AMIInterpolation::overlapSourceFaces() const
{
    return (*overlapSrcAddrPtr_);
}

inline const Foam::labelList&
Foam::AMIInterpolation::nonOverlapSourceFaces() const
{
    return (*nonOverlapSrcAddrPtr_);
}


inline Foam::label Foam::AMIInterpolation::singlePatchProc() const
{
    return singlePatchProc_;
}


inline Foam::scalar Foam::AMIInterpolation::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


inline bool Foam::AMIInterpolation::applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


inline const Foam::scalarField& Foam::AMIInterpolation::srcMagSf() const
{
    return srcMagSf_;
}


inline Foam::scalarField& Foam::AMIInterpolation::srcMagSf()
{
    return srcMagSf_;
}


inline const Foam::vectorField& Foam::AMIInterpolation::srcNf() const
{
    return srcNf_;
}


inline const Foam::labelListList& Foam::AMIInterpolation::srcAddress() const
{
    return srcAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolation::srcAddress()
{
    return srcAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolation::srcWeights() const
{
    return srcWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolation::srcWeights()
{
    return srcWeights_;
}


inline const Foam::scalarList& Foam::AMIInterpolation::srcSymWeights() const
{
    return srcSymWeights_;
}


inline Foam::scalarList& Foam::AMIInterpolation::srcSymWeights()
{
    return srcSymWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolation::srcWeightsSum() const
{
    return srcWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolation::srcWeightsSum()
{
    return srcWeightsSum_;
}


inline const Foam::pointListList& Foam::AMIInterpolation::srcCentroids() const
{
    return srcCentroids_;
}


inline Foam::pointListList& Foam::AMIInterpolation::srcCentroids()
{
    return srcCentroids_;
}


inline const Foam::distributionMap& Foam::AMIInterpolation::srcMap() const
{
    return srcMapPtr_();
}


inline const Foam::scalarField& Foam::AMIInterpolation::tgtMagSf() const
{
    return tgtMagSf_;
}


inline Foam::scalarField& Foam::AMIInterpolation::tgtMagSf()
{
    return tgtMagSf_;
}


inline const Foam::vectorField& Foam::AMIInterpolation::tgtNf() const
{
    return tgtNf_;
}


inline const Foam::labelListList& Foam::AMIInterpolation::tgtAddress() const
{
    return tgtAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolation::tgtAddress()
{
    return tgtAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolation::tgtWeights() const
{
    return tgtWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolation::tgtWeights()
{
    return tgtWeights_;
}


inline const Foam::scalarList& Foam::AMIInterpolation::tgtSymWeights() const
{
    return tgtSymWeights_;
}


inline Foam::scalarList& Foam::AMIInterpolation::tgtSymWeights()
{
    return tgtSymWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolation::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolation::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


inline const Foam::distributionMap& Foam::AMIInterpolation::tgtMap() const
{
    return tgtMapPtr_();
}


// ************************************************************************* //
