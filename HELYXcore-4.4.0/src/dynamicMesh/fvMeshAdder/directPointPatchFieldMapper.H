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
    (c) 2013 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::directPointPatchFieldMapper

Description
    direct pointPatchFieldMapper

\*---------------------------------------------------------------------------*/

#ifndef directPointPatchFieldMapper_H
#define directPointPatchFieldMapper_H

#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class directPointPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class directPointPatchFieldMapper
:
    public generalPointPatchFieldMapper
{
    // Private Data

        //- Addressing from new back to old
        const labelUList& directAddressing_;

        //- Does map contain any unmapped values
        bool hasUnmapped_;


public:

    // Constructors

        //- Construct given addressing
        directPointPatchFieldMapper(const labelUList& directAddressing)
        :
            directAddressing_(directAddressing),
            hasUnmapped_(false)
        {
            if (directAddressing_.size() && min(directAddressing_) < 0)
            {
                hasUnmapped_ = true;
            }
        }


    //- Destructor
    virtual ~directPointPatchFieldMapper()
    {}


    // Member Functions

        label size() const
        {
            return directAddressing_.size();
        }

        bool direct() const
        {
            return true;
        }

        bool indirect() const
        {
            return false;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        //- Return the direct map addressing
        const labelUList& directAddressing() const
        {
            return directAddressing_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
