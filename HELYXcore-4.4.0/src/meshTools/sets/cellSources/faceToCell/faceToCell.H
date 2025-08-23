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
    Foam::faceToCell

Description
    A topoSetSource to select cells based on usage in faces.

SourceFiles
    faceToCell.C

\*---------------------------------------------------------------------------*/

#ifndef faceToCell_H
#define faceToCell_H

#include "sets/topoSetSource/topoSetSource.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceToCell Declaration
\*---------------------------------------------------------------------------*/

class faceToCell
:
    public topoSetSource
{
public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            NEIGHBOUR,
            OWNER,
            ANY,
            ALL
        };

private:

        static const NamedEnum<faceAction, 4> faceActionNames_;


        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;

        //- Option
        faceAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("faceToCell");


    // Constructors

        //- Construct from components
        faceToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const faceAction option
        );

        //- Construct from dictionary
        faceToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        faceToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~faceToCell();


    // Member Functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
