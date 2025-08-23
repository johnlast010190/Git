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
    (c) 2016-2017 Wikki Ltd

Description
    Namespace of functions to calculate implicit derivatives returning a matrix.
    Time derivatives are calculated using Euler-implicit, backward differencing
    or Crank-Nicholson. Spatial derivatives are calculated using Gauss' Theorem.

\*---------------------------------------------------------------------------*/

#ifndef fam_H
#define fam_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "finiteArea/fam/famDdt.H"
#include "finiteArea/fam/famD2dt2.H"
#include "finiteArea/fam/famDiv.H"
#include "finiteArea/fam/famLaplacian.H"
#include "finiteArea/fam/famSup.H"
#include "finiteArea/fam/famNDiv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
