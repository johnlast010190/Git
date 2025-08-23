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
    fourthLnGrad

Description
    Fourth-order snGrad scheme with non-orthogonal correction.

SourceFiles
    fourthLnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef fourthLnGrad_H
#define fourthLnGrad_H

#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                 Class fourthLnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fourthLnGrad
:
    public lnGradScheme<Type>
{
    // Private data

        const word forceGradSchemeName_;


public:

    //- Runtime type information
    TypeName("fourth");


    // Constructors

        //- Construct from mesh
        fourthLnGrad
        (
            const faMesh& mesh, const word& forceGradSchemeName=word::null
        )
        :
            lnGradScheme<Type>(mesh),
            forceGradSchemeName_(forceGradSchemeName)
        {}


        //- Construct from mesh and data stream
        fourthLnGrad
        (
            const faMesh& mesh, Istream&, const word& forceGradSchemeName
        )
        :
            lnGradScheme<Type>(mesh),
            forceGradSchemeName_(forceGradSchemeName)
        {}


    // Destructor

        virtual ~fourthLnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> deltaCoeffs
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the fourthLnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fourthLnGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/lnGradSchemes/fourthLnGrad/fourthLnGrad.C"
#else
#   ifdef xlC
#   pragma implementation("fourthLnGrad.C")
#   endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
