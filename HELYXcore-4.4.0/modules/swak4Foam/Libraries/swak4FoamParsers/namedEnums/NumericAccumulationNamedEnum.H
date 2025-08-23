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
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::NumericAccumulationNamedEnum

Description

SourceFiles
    NumericAccumulationNamedEnum.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef NumericAccumulationNamedEnum_H
#define NumericAccumulationNamedEnum_H

#include "primitives/enums/NamedEnum.H"
#include "db/dictionary/dictionary.H"

#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class NumericAccumulationNamedEnum Declaration
\*---------------------------------------------------------------------------*/

class NumericAccumulationNamedEnum
{
    // Private data

public:

    enum value {
        numMin,
        numMax,
        numAverage,
        numWeightedAverage,
        numSum,
        numSumMag,
        numWeightedSum,
        numIntegrate,
        numMedian,
        numWeightedMedian,
        numQuantile,
        numWeightedQuantile,
        numRange,
        numWeightedRange,
        numSmaller,
        numWeightedSmaller,
        numBigger,
        numWeightedBigger,
        numWeightSum,
        numSize
    };

    typedef Tuple2<value,scalar> accuSpecification;

    //- String representation of order enums
    static const NamedEnum<value, 20> names;

    static List<accuSpecification> readAccumulations(
        const dictionary &dict,
        const word &name
    );

    static List<accuSpecification> readAccumulations(
        const wordList &names,
        const fileName &name
    );

    static word toString(const accuSpecification &);
};

Ostream &operator<<(
    Ostream &,
    const NumericAccumulationNamedEnum::accuSpecification &
);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
