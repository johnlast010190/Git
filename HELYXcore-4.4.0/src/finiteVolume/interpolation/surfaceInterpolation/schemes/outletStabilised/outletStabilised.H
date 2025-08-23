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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

Class
    Foam::outletStabilised

Group
    grpFvSurfaceInterpolationSchemes

Description
    Outlet-stabilised interpolation scheme which applies upwind differencing
    to the faces of the cells adjacent to outlets.

    This is particularly useful to stabilise the velocity at entrainment
    boundaries for LES cases using linear or other centred differencing
    schemes.

SourceFiles
    outletStabilised.C

\*---------------------------------------------------------------------------*/

#ifndef outletStabilised_H
#define outletStabilised_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/skewCorrected/skewCorrectionVectors.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchField.H"
#include "fields/fvPatchFields/basic/directionMixed/directionMixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class outletStabilised Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletStabilised
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const surfaceScalarField& faceFlux_;
        tmp<surfaceInterpolationScheme<Type>> tScheme_;


public:

    //- Runtime type information
    TypeName("outletStabilised");


    // Constructors

        //- Construct from mesh and Istream
        outletStabilised
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                db.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            ),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux_, is)
            )
        {}


        //- Construct from mesh, faceFlux and Istream
        outletStabilised
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        outletStabilised(const outletStabilised&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const VolField<Type>& vf
        ) const
        {
            tmp<surfaceScalarField> tw = tScheme_().weights(vf);
            surfaceScalarField& w = tw.ref();

            const fvMesh& mesh_ = this->mesh();
            const cellList& cells = mesh_.cells();

            forAll(vf.boundaryField(), patchi)
            {
                if
                (
                    isA<zeroGradientFvPatchField<Type>>
                        (vf.boundaryField()[patchi])
                 || isA<mixedFvPatchField<Type>>(vf.boundaryField()[patchi])
                 || isA<directionMixedFvPatchField<Type>>
                    (vf.boundaryField()[patchi])
                )
                {
                    const labelList& pFaceCells =
                        mesh_.boundary()[patchi].faceCells();

                    forAll(pFaceCells, pFacei)
                    {
                        const cell& pFaceCell = cells[pFaceCells[pFacei]];

                        forAll(pFaceCell, fi)
                        {
                            label facei = pFaceCell[fi];

                            if (mesh_.isInternalFace(facei))
                            {
                                // Apply upwind differencing
                                w[facei] = pos0(faceFlux_[facei]);
                            }
                        }
                    }
                }
            }

            return tw;
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme_().corrected();
        }

        //- Return the explicit correction to the face-interpolate
        //  set to zero on the near-boundary faces where upwinf is applied
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const
        {
            if (tScheme_().corrected())
            {
                tmp<SurfaceField<Type>> tcorr = tScheme_().correction(vf);

                SurfaceField<Type>& corr = tcorr.ref();

                const fvMesh& mesh_ = this->mesh();
                const cellList& cells = mesh_.cells();

                forAll(vf.boundaryField(), patchi)
                {
                    if
                    (
                        isA<zeroGradientFvPatchField<Type>>
                            (vf.boundaryField()[patchi])
                     || isA<mixedFvPatchField<Type>>
                            (vf.boundaryField()[patchi])
                    )
                    {
                        const labelList& pFaceCells =
                            mesh_.boundary()[patchi].faceCells();

                        forAll(pFaceCells, pFacei)
                        {
                            const cell& pFaceCell = cells[pFaceCells[pFacei]];

                            forAll(pFaceCell, fi)
                            {
                                label facei = pFaceCell[fi];

                                if (mesh_.isInternalFace(facei))
                                {
                                    // Remove correction
                                    corr[facei] = Zero;
                                }
                            }
                        }
                    }
                }

                return tcorr;
            }
            else
            {
                return tmp<SurfaceField<Type>>
                (
                    nullptr
                );
            }
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const outletStabilised&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
