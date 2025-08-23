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

Contributors/Copyright:
    2008-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "patchAverageFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
Field<T> patchAverageFunctionObject::average(const word& fieldName,T unsetVal) const
{
    const VolField<T>& fld =
        obr_.lookupObject<VolField<T>>
        (
            fieldName
        );

    Field<T> vals(patchNames_.size(), unsetVal);

    const fvMesh &mesh=refCast<const fvMesh>(obr_);

    forAll(patchNames_, patchI)
    {
        if (patchIndizes_[patchI] >= 0)
        {
            label index=patchIndizes_[patchI];
            scalar surf=sum(mesh.magSf().boundaryField()[index]);
            reduce(surf,sumOp<scalar>());
            vals[patchI] = sum
                        (
                            mesh.magSf().boundaryField()[index]
                            *fld.boundaryField()[index]
                        );
            reduce(vals[patchI],sumOp<T>());
            vals[patchI]/=surf;
        }
    }

    if (verbose()) {
        Info<< regionString()
            << " Averages of " << fieldName << " :";

        forAll(patchNames_, patchI)
        {
            Info<< "  " << patchNames_[patchI] << " = "
                << vals[patchI];
        }

        Info<< endl;
    }

    //    Pstream::listCombineGather(vals, isNotEqOp<T>());
    //    Pstream::listCombineScatter(vals);

    return vals;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
