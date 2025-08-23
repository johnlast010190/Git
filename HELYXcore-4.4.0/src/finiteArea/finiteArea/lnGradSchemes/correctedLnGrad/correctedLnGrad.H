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
    correctedLnGrad

Description
    Simple central-difference lnGrad scheme with non-orthogonal correction.

SourceFiles
    correctedLnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef correctedLnGrad_H
#define correctedLnGrad_H

#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                 Class correctedLnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class correctedLnGrad
:
    public lnGradScheme<Type>
{
    // Private data

        const word forceGradSchemeName_;


public:

    //- Runtime type information
    TypeName("corrected");


    // Constructors

        //- Construct from mesh
        correctedLnGrad
        (
            const faMesh& mesh, const word& forceGradSchemeName=word::null
        )
        :
            lnGradScheme<Type>(mesh),
            forceGradSchemeName_(forceGradSchemeName)
        {}


        //- Construct from mesh and data stream
        correctedLnGrad
        (
            const faMesh& mesh,
            Istream&,
            const word& forceGradSchemeName
        )
        :
            lnGradScheme<Type>(mesh),
            forceGradSchemeName_(forceGradSchemeName)
        {}


    // Destructor

        virtual ~correctedLnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> deltaCoeffs
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the correctedSnGrad
        //  for the given field using the gradient of the field
        tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        fullGradCorrection
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        //- Return the explicit correction to the correctedLnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const correctedLnGrad&) = delete;

};
// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<edgeScalarField> correctedLnGrad<scalar>::correction
(
    const areaScalarField& vsf
) const;

template<>
tmp<edgeVectorField> correctedLnGrad<vector>::correction
(
    const areaVectorField& vvf
) const;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
