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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::localMin

Group
    grpFvSurfaceInterpolationSchemes

Description
    LocalMin-mean differencing scheme class.

    This scheme interpolates 1/field using a scheme specified at run-time
    and return the reciprocal of the interpolate.

SourceFiles
    localMin.C

\*---------------------------------------------------------------------------*/

#ifndef localMin_H
#define localMin_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localMin Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class localMin
:
    public surfaceInterpolationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("localMin");


    // Constructors

        //- Construct from mesh
        localMin(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        localMin
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        localMin
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights(const VolField<Type>&) const
        {
            NotImplemented;
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<SurfaceField<Type>> interpolate
        (
            const VolField<Type>& vf
        ) const
        {
            const fvMesh& mesh = vf.mesh();

            tmp<SurfaceField<Type>> tvff
            (
                SurfaceField<Type>::New
                (
                    "localMin::interpolate(" + vf.name() + ')',
                    mesh,
                    vf.dimensions()
                )
            );
            SurfaceField<Type>& vff = tvff.ref();

            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();

            forAll(vff, facei)
            {
                vff[facei] = min(vf[own[facei]], vf[nei[facei]]);
            }

            typename SurfaceField<Type>::
                Boundary& bff = vff.boundaryFieldRef();

            forAll(bff, patchi)
            {
                const fvPatchField<Type>& pf = vf.boundaryField()[patchi];
                Field<Type>& pff = bff[patchi];

                if (pf.coupled())
                {
                    tmp<Field<Type>> tpif(pf.patchInternalField());
                    const Field<Type>& pif = tpif();

                    tmp<Field<Type>> tpnf(pf.patchNeighbourField());
                    const Field<Type>& pnf = tpnf();

                    forAll(pff, i)
                    {
                        pff[i] = min(pif[i], pnf[i]);
                    }
                }
                else
                {
                    pff = pf;
                }
            }

            return tvff;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const localMin&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
