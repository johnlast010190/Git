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
    (c) 2024 Engys Ltd.

Class
    manipulateFvSolutionFvSchemesFunctionObject

Description
    Basis for functionObjects that manipulate the data from
    - fvScheme
    - fvSolution
    in memory. Provides a framework and utility methods but the actual
    logic is in the subclasses

SourceFiles
    manipulateFvSolutionFvSchemesFunctionObject.C

Contributors/Copyright:
    2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef manipulateFvSolutionFvSchemesFunctionObject_H
#define manipulateFvSolutionFvSchemesFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"
#include "finiteVolume/fvSchemes/fvSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvSolution;

/*---------------------------------------------------------------------------*\
        Class manipulateFvSolutionFvSchemesFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class manipulateFvSolutionFvSchemesFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        fvSolution& fvSolution_;

        fvSchemes& fvSchemes_;

        dictionary fvSolutionBackup_;

        dictionary fvSchemesBackup_;


    // Private Member Functions

        void rereadFvSolution();

        void rereadFvSchemes();


protected:

    // Protected Member Functions

        void resetFvSolution();

        void resetFvSchemes();

        void resetDict(dictionary &dict);

        dictionary& fvSolutionDict();

        dictionary& fvSchemesDict();


public:

    //- Runtime type information
    TypeName("manipulateFvSolutionFvSchemes");


    // Constructors

        //- Construct from components
        manipulateFvSolutionFvSchemesFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        manipulateFvSolutionFvSchemesFunctionObject
        (
            const manipulateFvSolutionFvSchemesFunctionObject&
        ) = delete;


    // Member Functions

        bool start();

        void writeSimple();

        virtual bool manipulateFvSolution(const Time&) = 0;

        virtual bool manipulateFvSchemes(const Time&) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const manipulateFvSolutionFvSchemesFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
