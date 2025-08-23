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
    Foam::fv::maskBlendedGrad

Description
    blended gradient scheme
    Formal: namedBlended name <nameField> <scheme1> <scheme2>

SourceFiles
    maskBlendedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef maskBlendedGrad_H
#define maskBlendedGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class maskBlendedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class maskBlendedGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- Scheme 1
        tmp<gradScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<gradScheme<Type>> tScheme2_;


public:

    //- Runtime type information
    TypeName("maskBlended");


    // Constructors

        //- Construct from Istream
        maskBlendedGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
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
        maskBlendedGrad(const maskBlendedGrad&) = delete;


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad(const VolField<Type>& vsf, const word& name) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const maskBlendedGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteVolume/gradSchemes/maskBlendedGrad/maskBlendedGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
