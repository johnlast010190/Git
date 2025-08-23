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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::polyRemoveCell

Description
    Class containing data for cell removal.

\*---------------------------------------------------------------------------*/

#ifndef polyRemoveCell_H
#define polyRemoveCell_H

#include "primitives/ints/label/label.H"
#include "polyTopoChange/polyTopoChange/topoAction/topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class polyRemoveCell Declaration
\*---------------------------------------------------------------------------*/

class polyRemoveCell
:
    public topoAction
{
    // Private data

        //- Cell ID
        label cellID_;

        //- Merge cell ID or -1
        label mergeCellID_;

public:

    // Static data members

        //- Runtime type information
        TypeName("removeCell");


    // Constructors

        //- Construct null.  Used for constructing lists
        polyRemoveCell()
        :
            cellID_(-1),
            mergeCellID_(-1)
        {}

        //- Construct from components
        polyRemoveCell(const label cellID, const label mergeCellID = -1)
        :
            cellID_(cellID),
            mergeCellID_(mergeCellID)
        {}

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyRemoveCell(*this));
        }


    // Default Destructor


    // Member Functions

        //- Return cell ID
        label cellID() const
        {
            return cellID_;
        }

        //- Return cell ID
        label mergeCellID() const
        {
            return mergeCellID_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
