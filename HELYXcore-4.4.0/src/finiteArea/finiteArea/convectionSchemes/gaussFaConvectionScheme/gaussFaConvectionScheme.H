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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016-2024 Engys Ltd.

Class
    gaussConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    gaussFaConvectionScheme.C
    gaussFaConvectionSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gaussFaConvectionScheme_H
#define gaussFaConvectionScheme_H

#include "finiteArea/convectionSchemes/faConvectionScheme/faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

//- Temporary debug switch to provide warning about backward-compatibility
//  issue with setting div schemes for steady-state
extern int warnUnboundedGauss;


/*---------------------------------------------------------------------------*\
                       Class gaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussConvectionScheme
:
    public fa::convectionScheme<Type>
{
    // Private data

        tmp<edgeInterpolationScheme<Type>> tinterpScheme_;


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from flux and interpolation scheme
        gaussConvectionScheme
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            const tmp<edgeInterpolationScheme<Type>>& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        gaussConvectionScheme
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                edgeInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        gaussConvectionScheme(const gaussConvectionScheme&) = delete;


    //- Destructor
    virtual ~gaussConvectionScheme()
    {}


    // Member Functions

        const edgeInterpolationScheme<Type>& interpScheme() const;

        tmp<GeometricField<Type, faePatchField, faEdgeMesh>> flux
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        tmp<faMatrix<Type>> famDiv
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDiv
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussConvectionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/convectionSchemes/gaussFaConvectionScheme/gaussFaConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
