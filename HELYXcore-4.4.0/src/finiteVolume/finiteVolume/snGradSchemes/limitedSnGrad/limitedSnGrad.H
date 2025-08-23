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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd

Class
    Foam::fv::limitedSnGrad

Group
    grpFvSnGradSchemes

Description
    Run-time selected snGrad scheme with limited non-orthogonal correction.

    The limiter is controlled by a coefficient with a value between 0 and 1
    which when 0 switches the correction off and the scheme behaves as
    uncorrectedSnGrad, when set to 1 the full correction of the selected scheme
    is used and when set to 0.5 the limiter is calculated such that the
    non-orthogonal contribution does not exceed the orthogonal part.

    Format:
        limited \<corrected scheme\> \<coefficient\>;

        or

        limited \<coefficient\>;  // Backward compatibility

SourceFiles
    limitedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef limitedSnGrad_H
#define limitedSnGrad_H

#include "finiteVolume/snGradSchemes/correctedSnGrad/correctedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class limitedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedSnGrad
:
    public snGradScheme<Type>
{
    // Private data

        tmp<snGradScheme<Type>> correctedScheme_;

        scalar limitCoeff_;


    // Private Member Functions

        //- Lookup function for the corrected to support backward compatibility
        //  of dictionary specification
        tmp<snGradScheme<Type>>
        lookupCorrectedScheme
        (
            Istream& schemeData, const word& forceGradSchemeName
        )
        {
            token nextToken(schemeData);

            if (nextToken.isNumber())
            {
                limitCoeff_ = nextToken.number();
                return tmp<snGradScheme<Type>>
                (
                    new correctedSnGrad<Type>
                    (
                        this->mesh(), forceGradSchemeName
                    )
                );
            }
            else
            {
                schemeData.putBack(nextToken);
                tmp<snGradScheme<Type>> tcorrectedScheme
                (
                    fv::snGradScheme<Type>::New
                    (
                        this->mesh(), schemeData, forceGradSchemeName
                    )
                );

                schemeData >> limitCoeff_;

                return tcorrectedScheme;
            }
        }


public:

    //- Runtime type information
    TypeName("limited");


    // Constructors

        //- Construct from mesh
        limitedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh),
            correctedScheme_(new correctedSnGrad<Type>(this->mesh())),
            limitCoeff_(1)
        {}


        //- Construct from mesh and data stream
        limitedSnGrad
        (
            const fvMesh& mesh,
            Istream& schemeData,
            const word& forceGradSchemeName=word::null
        )
        :
            snGradScheme<Type>(mesh),
            correctedScheme_
            (
                lookupCorrectedScheme(schemeData, forceGradSchemeName)
            )
        {
            if (limitCoeff_ < 0 || limitCoeff_ > 1)
            {
                FatalIOErrorInFunction
                (
                    schemeData
                )   << "limitCoeff is specified as " << limitCoeff_
                    << " but should be >= 0 && <= 1"
                    << exit(FatalIOError);
            }
        }


    //- Destructor
    virtual ~limitedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        virtual scalar limitCoeff() const
        {
            return limitCoeff_;
        }

        //- Return the explicit correction to the limitedSnGrad
        //  for the given field
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitedSnGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/snGradSchemes/limitedSnGrad/limitedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
