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
    (c) 2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::meshPhiPreCorrectInfo

SourceFiles
    meshPhiPreCorrectInfoI.H

\*---------------------------------------------------------------------------*/

#ifndef meshPhiPreCorrectInfo_H
#define meshPhiPreCorrectInfo_H

#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvPatch;
class fvMesh;
class transformer;

// Forward declaration of friend functions and operators
class meshPhiPreCorrectInfo;
Ostream& operator<<(Ostream&, const meshPhiPreCorrectInfo&);
Istream& operator>>(Istream&, meshPhiPreCorrectInfo&);

/*---------------------------------------------------------------------------*\
                    Class meshPhiPreCorrectInfo Declaration
\*---------------------------------------------------------------------------*/

class meshPhiPreCorrectInfo
{
    // Private Data

        //- Layer index
        label layer_;

        //- Weight
        scalar weight_;


public:

    // Constructors

        //- Construct null
        inline meshPhiPreCorrectInfo();

        //- Construct given the layer index and weight
        inline meshPhiPreCorrectInfo(const label layer, const scalar weight);


    // Member Functions

        // Access

            //- Return the layer index
            inline label layer() const;

            //- Return the weight
            inline scalar weight() const;


        // Needed by FaceCellWave

            //- Check whether the meshPhiPreCorrectInfo has been changed at all
            //  or still contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for checking
            //  consistency across cyclics.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const fvMesh& mesh,
                const meshPhiPreCorrectInfo&,
                const scalar,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const fvPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            )
            {}

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh& mesh,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const meshPhiPreCorrectInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const label neighbourCelli,
                const meshPhiPreCorrectInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const meshPhiPreCorrectInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const meshPhiPreCorrectInfo&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const meshPhiPreCorrectInfo&) const;
        inline bool operator!=(const meshPhiPreCorrectInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const meshPhiPreCorrectInfo&);
        friend Istream& operator>>(Istream&, meshPhiPreCorrectInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
constexpr bool contiguous<meshPhiPreCorrectInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moving/meshPhiPreCorrectInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
