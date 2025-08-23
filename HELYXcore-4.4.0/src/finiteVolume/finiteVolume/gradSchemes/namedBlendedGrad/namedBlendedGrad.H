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
    Foam::fv::namedBlendedGrad

Description
    blended gradient scheme
    Formal: namedBlended name <nameField> <scheme1> <scheme2>

SourceFiles
    namedBlendedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef namedBlendedGrad_H
#define namedBlendedGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class namedBlendedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class namedBlendedGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions
        word blendingFieldName_;

        //- Scheme 1
        tmp<gradScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<gradScheme<Type>> tScheme2_;


public:

    //- Runtime type information
    TypeName("namedBlended");


    // Constructors

        //- Construct from Istream
        namedBlendedGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
            blendingFieldName_(is),
            tScheme1_
            (
                gradScheme<Type>::New(mesh, db, is)
            ),
            tScheme2_
            (
                gradScheme<Type>::New(mesh, db, is)
            )
        {}

        //- Disallow default bitwise copy construction
        namedBlendedGrad(const namedBlendedGrad&) = delete;


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad(const VolField<Type>& vsf, const word& name) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const namedBlendedGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteVolume/gradSchemes/namedBlendedGrad/namedBlendedGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
