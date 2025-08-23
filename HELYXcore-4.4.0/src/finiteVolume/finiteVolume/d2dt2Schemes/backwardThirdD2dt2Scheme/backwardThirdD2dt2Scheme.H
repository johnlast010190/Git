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
    Foam::fv::backwardThirdD2dt2Scheme

Description
    Blended Second-Third order implicit backward-differencing d2dt2
    using the current and six previous time-step values.

SourceFiles
    backwardThirdD2dt2Scheme.C

Author
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef backwardThirdD2dt2Scheme_H
#define backwardThirdD2dt2Scheme_H

#include "finiteVolume/d2dt2Schemes/d2dt2Scheme/d2dt2Scheme.H"
#include "finiteVolume/d2dt2Schemes/backwardD2dt2Scheme/backwardD2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class backwardThirdD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class backwardThirdD2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
    // Private Member Functions

        //- Off centering coefficient to blend second and third-order
        //  backward
        autoPtr<Function1<scalar>> ocCoeff_;

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        scalar deltaT0_(const VolField<Type>&) const;


public:

    //- Runtime type information
    TypeName("backwardThird");


    // Constructors

        //- Construct from mesh
        backwardThirdD2dt2Scheme(const fvMesh& mesh);

        //- Construct from mesh and Istream
        backwardThirdD2dt2Scheme(const fvMesh& mesh, Istream& is);

        //- Disallow default bitwise copy construction
        backwardThirdD2dt2Scheme(const backwardThirdD2dt2Scheme&) = delete;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

        //- Return the current off-centreing coefficient
        scalar ocCoeff() const
        {
            return ocCoeff_->value(mesh().time().value());
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
        void operator=(const backwardThirdD2dt2Scheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardThirdD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
