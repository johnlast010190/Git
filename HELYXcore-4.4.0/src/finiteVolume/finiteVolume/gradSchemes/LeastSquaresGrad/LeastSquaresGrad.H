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
    (c) 2013-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::LeastSquaresGrad

Group
    grpFvGradSchemes

Description
    Gradient calculated using weighted least-squares on an arbitrary stencil.
    The stencil type is provided via a template argument and any cell-based
    stencil is supported:

    \table
        Stencil                     | Connections     | Scheme name
        centredCFCCellToCellStencil | cell-face-cell  | Not Instantiated
        centredCPCCellToCellStencil | cell-point-cell | pointCellsLeastSquares
        centredCECCellToCellStencil | cell-edge-cell  | edgeCellsLeastSquares
    \endtable

    The first of these is not instantiated by default as the standard
    leastSquaresGrad is equivalent and more efficient.

Usage
    Example of the gradient specification:
    \verbatim
    gradSchemes
    {
        default         pointCellsLeastSquares;
    }
    \endverbatim

See also
    Foam::fv::LeastSquaresVectors
    Foam::fv::leastSquaresGrad

SourceFiles
    LeastSquaresGrad.C
    LeastSquaresVectors.H
    LeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef LeastSquaresGrad_H
#define LeastSquaresGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class LeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Stencil>
class LeastSquaresGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        //- Minimum determinant criterion to choose extra cells
        scalar minDet_;


public:

    //- Runtime type information
    TypeName("LeastSquares");


    // Constructors

        //- Construct from Istream
        LeastSquaresGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            gradScheme<Type>(mesh)
        {}

        //- Disallow default bitwise copy construction
        LeastSquaresGrad(const LeastSquaresGrad&) = delete;


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad
        (
            const VolField<Type>& vsf,
            const word& name
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LeastSquaresGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeLeastSquaresGradTypeScheme(SS, STENCIL, TYPE)                      \
    typedef Foam::fv::LeastSquaresGrad<Foam::TYPE, Foam::STENCIL>              \
        LeastSquaresGrad##TYPE##STENCIL##_;                                    \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
        (LeastSquaresGrad##TYPE##STENCIL##_, #SS, 0);                          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            typedef LeastSquaresGrad<Foam::TYPE, Foam::STENCIL>                \
                LeastSquaresGrad##TYPE##STENCIL##_;                            \
                                                                               \
            gradScheme<Foam::TYPE>::addIstreamConstructorToTable               \
                <LeastSquaresGrad<Foam::TYPE, Foam::STENCIL>>                  \
                add##SS##STENCIL##TYPE##IstreamConstructorToTable_;            \
        }                                                                      \
    }

#define makeLeastSquaresGradScheme(SS, STENCIL)                                \
    typedef Foam::fv::LeastSquaresVectors<Foam::STENCIL>                       \
        LeastSquaresVectors##STENCIL##_;                                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
        (LeastSquaresVectors##STENCIL##_, #SS, 0);                             \
                                                                               \
    makeLeastSquaresGradTypeScheme(SS,STENCIL,scalar)                          \
    makeLeastSquaresGradTypeScheme(SS,STENCIL,vector)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/gradSchemes/LeastSquaresGrad/LeastSquaresGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
