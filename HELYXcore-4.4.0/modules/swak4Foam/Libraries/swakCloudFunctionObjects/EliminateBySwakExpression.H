/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::EliminateBySwakExpression

Description
    Eliminate parcels for which a swak-Expression is true

SourceFiles
    EliminateBySwakExpression.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@ice-sf.at>

    Based on the CloudTracks-function object

\*---------------------------------------------------------------------------*/

#ifndef EliminateBySwakExpression_H
#define EliminateBySwakExpression_H

#include "submodels/CloudFunctionObjects/CloudFunctionObject/CloudFunctionObject.H"

#include "include/swakCloudTypes.H"

#include "DataStructures/TimelineCollection/TimelineCollection.H"

#include "CloudValueExpressionDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EliminateBySwakExpression Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EliminateBySwakExpression
:
    public CloudFunctionObject<CloudType>
{
    // Private data

    // Typedefs

    //- Convenience typedef for parcel type
    typedef typename CloudType::parcelType parcelType;

    //- Pointer to the cloud storage
    autoPtr<Cloud<parcelType>> eliminatedPrePtr_;

    //- Pointer to the cloud storage
    autoPtr<Cloud<parcelType>> eliminatedPostPtr_;

    bool eliminatePre_;

    bool eliminatePost_;

    //- here we get our files from
    TimelineCollection out_;

    //- the actual expression driver
    CloudValueExpressionDriver driver_;

    //- the logical expression
    exprString expression_;

    //- check whether particles are outside and mark particles to be outside
    label checkExpression(Cloud<parcelType> &eliminate);

    label nrPre_;

protected:

    // Protected member functions

    //- Write post-processing info
    void write();


public:

    //- Runtime type information
    TypeName("eliminateBySwakExpression");


    // Constructors

    //- Construct from dictionary
    EliminateBySwakExpression
    (
        const dictionary& dict,
        CloudType& owner,
        const word& modelName
    );

    //- Construct copy
    EliminateBySwakExpression(const EliminateBySwakExpression<CloudType>& ppm);

    //- Construct and return a clone
    virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
                (
                    new EliminateBySwakExpression<CloudType>(*this)
                );
        }


    //- Destructor
    virtual ~EliminateBySwakExpression();


    // Member Functions

    // Evaluation

    //- Pre-evolve hook
    virtual void preEvolve();

    //- Post-evolve hook
    virtual void postEvolve();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EliminateBySwakExpression.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
