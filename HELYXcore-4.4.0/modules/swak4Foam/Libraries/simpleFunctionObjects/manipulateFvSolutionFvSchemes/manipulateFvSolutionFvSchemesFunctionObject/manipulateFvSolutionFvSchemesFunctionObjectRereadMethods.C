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

Descriptions:
    Encapsulates dirty tricks to access the private read-methods of
    Foam::solution and Foam::fvSchemes

Contributors/Copyright:
    2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/


// not working

#include "include/swak.H"

#ifdef __clang__

#ifndef FOAM_SOLUTION_HAS_NO_READ_WITH_DICT

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "finiteVolume/fvSchemes/fvSchemes.H"
#include "matrices/solution/solution.H"

#undef debug

#endif

#endif

#include "manipulateFvSolutionFvSchemesFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fvMesh/fvMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

void manipulateFvSolutionFvSchemesFunctionObject::rereadFvSolution()
{
#ifndef FOAM_SOLUTION_HAS_NO_READ_WITH_DICT
    if (debug) {
        Info<< "Rereading " << fvSolution_.name() << endl;
    }

    fvSolution_.read(fvSolutionDict());

#endif
}

void manipulateFvSolutionFvSchemesFunctionObject::rereadFvSchemes()
{
#ifndef FOAM_SCHEMES_HAS_NO_READ_WITH_DICT
    if (debug) {
        Info<< "Rereading " << fvSchemes_.name() << endl;
    }

    fvSchemes_.read(fvSchemesDict());
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
