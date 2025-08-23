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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::volFieldWeighted

Description
    Interpolation scheme class that generates weightes from a specified scalar
    field looked up from the registry.

    Weights are calculated inverse:

    w = A[neightbour]/(A[owner] + A[neighbour])

SourceFiles
    volFieldWeighted.C

\*---------------------------------------------------------------------------*/

#ifndef volFieldWeighted_H
#define volFieldWeighted_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class volFieldWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volFieldWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const volScalarField& vf_;


public:

    //- Runtime type information
    TypeName("volFieldWeighted");


    // Constructors

        //- Construct from Istream
        volFieldWeighted
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            vf_
            (
                db.objectRegistry::template
                lookupObject<const volScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        volFieldWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            vf_
            (
                this->mesh().objectRegistry::template
                lookupObject<const volScalarField>
                (
                    word(is)
                )
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights(const VolField<Type>&) const
        {
            tmp<surfaceScalarField> volFieldWeights
            (
                surfaceScalarField::New
                (
                    "volFieldWeights",
                    this->mesh(),
                    dimless
                )
            );

            const labelUList& owner = this->mesh().owner();
            const labelUList& neighbour = this->mesh().neighbour();


            // ... and reference to the internal field of the weighting factors
            scalarField& w = volFieldWeights.ref().primitiveFieldRef();

            forAll(owner, facei)
            {
                w[facei]
                    = vf_[neighbour[facei]]
                    /(vf_[owner[facei]] + vf_[neighbour[facei]]);
            }

            forAll(this->mesh().boundary(), patchI)
            {
                if (volFieldWeights->boundaryField()[patchI].coupled())
                {

                    if (Pstream::parRun())
                    {
                        const scalarField vfNei
                        (
                            vf_.boundaryField()[patchI].patchNeighbourField()
                        );

                        const scalarField vfOwn
                        (
                            vf_.boundaryField()[patchI].patchInternalField()
                        );

                        volFieldWeights.ref().boundaryFieldRef()[patchI]
                            = vfNei / (vfNei + vfOwn);

                    }
                    else
                    {
                        volFieldWeights.ref().boundaryFieldRef()[patchI] = 1.0;
                    }

                }
                else
                {
                    volFieldWeights.ref().boundaryFieldRef()[patchI] = 1.0;
                }
            }

            return volFieldWeights;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volFieldWeighted&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
