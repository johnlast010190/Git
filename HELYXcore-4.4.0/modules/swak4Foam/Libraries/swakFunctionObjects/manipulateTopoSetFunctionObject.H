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
    manipulateTopoSetFunctionObject

Description
    Function object that manipulates a topoSet

SourceFiles
    manipulateTopoSetFunctionObject.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef manipulateTopoSetFunctionObject_H
#define manipulateTopoSetFunctionObject_H

#include "manipulateTopoSet.H"
#include "OutputFilterFunctionObject.H"

#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef OutputFilterFunctionObject<manipulateTopoSet<cellSet>> manipulateCellSetFunctionObject;
    typedef OutputFilterFunctionObject<manipulateTopoSet<faceSet>> manipulateFaceSetFunctionObject;
    typedef OutputFilterFunctionObject<manipulateTopoSet<pointSet>> manipulatePointSetFunctionObject;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
