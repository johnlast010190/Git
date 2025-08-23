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
    (c) 2024 Engys Ltd.

Class
    Foam::Histogram

Description
    Calculates the counts per bin of a list.

SourceFiles
    Histogram.C

\*---------------------------------------------------------------------------*/

#ifndef Histogram_H
#define Histogram_H

#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                          Class Histogram Declaration
\*---------------------------------------------------------------------------*/

template<class List>
class Histogram
{
    // Private data

        //- Counts per bin
        labelList counts_;

        //- Number of <= lowest bin
        label nLow_;

        //- Number of > highest bin
        label nHigh_;


    // Private Member Functions

        void count(const List& bins, const List& l);


public:

    // Constructors

        //- Construct given bin values and input list
        Histogram(const List& bins, const List& l);

        //- Construct given min, max, number of bins and input list
        Histogram
        (
            const typename List::const_reference min,
            const typename List::const_reference max,
            const label nBins,
            const List& l
        );

        //- Disallow default bitwise copy construction
        Histogram(const Histogram&) = delete;


        // Access

            //- Return the counts per bin
            inline const labelList& counts() const
            {
                return counts_;
            }

            //- Return the number of elements <= bins[0]
            //  (so inclusive lowest bin value)
            inline label nLow() const
            {
                return nLow_;
            }

            //- Return the number of elements > bins[bins.size()-1]
            //  (so exclusive highest bin value)
            inline label nHigh() const
            {
                return nHigh_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Histogram&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//    #include "HistogramI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "containers/Lists/Histogram/Histogram.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
