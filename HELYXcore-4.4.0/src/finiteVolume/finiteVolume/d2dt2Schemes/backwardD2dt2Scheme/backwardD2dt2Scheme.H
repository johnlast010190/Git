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

Class
    Foam::fv::backwardD2dt2Scheme

Description
    Second-order backward implicit d2dt2 using the current and four previous
    time-step values.

SourceFiles
    backwardD2dt2Scheme.C

Author
    (c) Zeljko Tukovic
        Some additions by Philip Cardiff (solids4foam)
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef backwardD2dt2Scheme_H
#define backwardD2dt2Scheme_H

#include "finiteVolume/d2dt2Schemes/d2dt2Scheme/d2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class backwardD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class backwardD2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        scalar deltaT0_(const VolField<Type>&) const;


public:

    //- Runtime type information
    TypeName("backward");


    // Constructors

        //- Construct from mesh
        backwardD2dt2Scheme(const fvMesh& mesh)
        :
            d2dt2Scheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        backwardD2dt2Scheme(const fvMesh& mesh, Istream& is)
        :
            d2dt2Scheme<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        backwardD2dt2Scheme(const backwardD2dt2Scheme&) = delete;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

        tmp<VolField<Type>> fvcD2dt2(const VolField<Type>&);

        tmp<VolField<Type>> fvcD2dt2
        (
            const volScalarField&,
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmD2dt2(const VolField<Type>&);

        tmp<fvMatrix<Type>> fvmD2dt2
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmD2dt2
        (
            const volScalarField&,
            const VolField<Type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const backwardD2dt2Scheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
