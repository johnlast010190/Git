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

Class
    Foam::expressionToCell

Description
    A topoSetSource to select the cell from a logical expression

SourceFiles
    expressionToCell.C

Contributors/Copyright:
    2010-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef expressionToCell_H
#define expressionToCell_H

#include "sets/topoSetSource/topoSetSource.H"

#include "include/swak.H"

#include "helpers/exprString.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class expressionToCell Declaration
\*---------------------------------------------------------------------------*/

class expressionToCell
:
    public topoSetSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- The logical expression to use
        exprString expression_;

        //- dictionary with additional data
        autoPtr<dictionary> dict_;

        void combine(topoSet& set, const bool add) const;

public:

    //- Runtime type information
    TypeName("expressionToCell");

    // Constructors

        //- Construct from components
        expressionToCell
        (
            const polyMesh& mesh,
            const exprString& expression
        );

        //- Construct from dictionary
        expressionToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        expressionToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~expressionToCell();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

#ifdef FOAM_TOPOSETSOURCE_HAS_SETTYPE
        virtual sourceType setType() const;
#endif

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
