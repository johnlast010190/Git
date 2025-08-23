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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2010-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh/surfaceMesh.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"

#include <cassert>

namespace Foam {
    label CellSetValueExpressionDriver::size() const
    {
        return cellSet_->size();
    }

    label CellSetValueExpressionDriver::pointSize() const
    {
        notImplemented("CellSetValueExpressionDriver::pointSize()");
    }

    const fvMesh &CellSetValueExpressionDriver::mesh() const
    {
        return dynamic_cast<const fvMesh&>(cellSet_->db());
        //        return dynamicCast<const fvMesh&>(cellSet_->db()); // doesn't wor with gcc 4.2
    }

    const objectRegistry& CellSetValueExpressionDriver::db() const
    {
        return dynamic_cast<const objectRegistry&>(cellSet_->db());
        //        return dynamicCast<const fvMesh&>(cellSet_->db()); // doesn't wor with gcc 4.2
    }
}


// ************************************************************************* //
