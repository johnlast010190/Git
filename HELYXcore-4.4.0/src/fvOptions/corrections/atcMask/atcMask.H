/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX® : Professional Open-source CFD
|   o   O   o    |
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------

License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM® <http://www.openfoam.org/>.

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
    (c) 2018-2024 Engys Ltd.

Author
    2018. Nikolaos Magoulas (Engys Ltd.). All rights reserved.

InClass
    Foam::fv::atcMask

Description
    Masks the ATC field according to a sensor field.

Usage
    Example usage:
    \verbatim
    mask
    {
        type            atcMask;
        active          yes;

        atcMaskCoeffs
        {
            selectionMode   all;
            fieldName       atcMask_drag
            max             1e10;
        }

        sensor
        {
            type        maxSingularValue;
            fieldName   U;
        }
    }
    \endverbatim

SourceFiles
    atcMask.C

\*---------------------------------------------------------------------------*/

#ifndef atcMask_H
#define atcMask_H

#include "cellSetOption/cellSetOption.H"
#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class atcMask Declaration
\*---------------------------------------------------------------------------*/

class atcMask
:
    public cellSetOption
{
protected:

    // Protected data

        //- Name of field to be masked
        word fieldName_;

        //- Sensor function to be used as input for the masking
        autoPtr<sensor<vector>> sensor_;

        //- Minimum value of sensor field for which masking doesn't an effect
        scalar min_;

        //- Maximum value of sensor field after which masking is fully applied
        scalar max_;

        //- Write mask
        Switch writeMask_;


public:

    //- Runtime type information
    TypeName("atcMask");


    // Constructors

        //- Construct from components
        atcMask
        (
            const word& name,
            const word& modelType,
            const dictionary& optionDict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        atcMask(const atcMask&) = delete;


    //- Destructor
    virtual ~atcMask()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the energy field
        virtual void correct(volScalarField& atcMask);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const atcMask&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
