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
    (c) 2024 Engys Ltd.

Class
    Foam::fv::leastSquaresGrad

Group
    grpFvGradSchemes

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    leastSquaresGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresGrad_H
#define leastSquaresGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"
#include "blockLduSystem/BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class leastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresGrad
:
    public fv::gradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("leastSquares");


    // Constructors

        //- Construct from mesh
        leastSquaresGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream&
        )
        :
            gradScheme<Type>(mesh)
        {}

        //- Disallow default bitwise copy construction
        leastSquaresGrad(const leastSquaresGrad&) = delete;


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad
        (
            const VolField<Type>& vsf,
            const word& name
        ) const;

        tmp
        <
            BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
        > fvmGrad
        (
            const VolField<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/gradSchemes/leastSquaresGrad/leastSquaresGrad.C"
#endif

// Specialisation for scalars
#include "finiteVolume/gradSchemes/leastSquaresGrad/scalarLeastSquaresGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
