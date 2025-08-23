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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::setsToFaceZone

Description
    A topoSetSource to select faces based on usage in a faceSet and cellSet

SourceFiles
    setsToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setsToFaceZone_H
#define setsToFaceZone_H

#include "sets/topoSetSource/topoSetSource.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setsToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setsToFaceZone
:
    public topoSetSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        const word faceSetName_;

        //- Name of set to use
        const word cellSetName_;

        //- Whether cellSet is slave cells or master cells
        const Switch flip_;

public:

    //- Runtime type information
    TypeName("setsToFaceZone");

    // Constructors

        //- Construct from components
        setsToFaceZone
        (
            const polyMesh& mesh,
            const word& faceSetName,
            const word& cellSetName,
            const Switch& flip
        );

        //- Construct from dictionary
        setsToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        setsToFaceZone
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~setsToFaceZone();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACEZONESOURCE;
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
