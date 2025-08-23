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
    (c) 2009-2025 Engys Ltd.

Class
    Foam::dataItem

Description
    A database item object for machine-learning.

SourceFiles
    dataItem.C

\*---------------------------------------------------------------------------*/

#ifndef dataItem_H
#define dataItem_H

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dataItem Declaration
\*---------------------------------------------------------------------------*/

class dataItem
{
protected:

    // Protected Data

        word caseName_;
        scalar rhoRef_;
        vector Uref_;
        word caseType_;
        word caseDir_;

        //- Label of the case as training class label
        scalar identity_;

        scalar xmin_;
        scalar ymin_;
        scalar zmin_;
        scalar xmax_;
        scalar ymax_;
        scalar zmax_;

        label m_;
        label n_;
        label k_;


public:

    // Public Data

        std::vector<scalar> geoData_;


    // Constructors

        //- Construct null
        dataItem();

        //- Construct from string
        dataItem(const word& dat);


    // Member Functions

        void setDataInfo(const word& dat);

        label m() const
        {
            return m_;
        }

        label n() const
        {
            return n_;
        }

        label k() const
        {
            return k_;
        }

        void setData
        (
            label m,
            label n,
            label k,
            const std::vector<scalar>& data
        );

        scalar getDist(const std::vector<scalar>& vdata) const;

        scalar gdataSum() const;

        word path() const
        {
            return caseDir_;
        }

        void setData
        (
            const word& data
        );

        label numSolids() const;

        scalar identity() const
        {
            return identity_;
        }

        const word& caseName() const
        {
            return caseName_;
        }

        void identity(scalar d)
        {
            identity_ = d;
        }

        void show() const;

        void clearGeoData()
        {
            geoData_.resize(0);
        }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline void split(wordList& vect, const word& str)
{
	vect.resize(0);

    std::istringstream iss(str);
    for (word w; iss >> w;)
    {
        vect.append(w);
    }
}


inline void split(wordList& elems, const word& s, char delim)
{
    std::stringstream ss(s);

    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.append(item);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
