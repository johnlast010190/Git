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
    (c) held by original author

Contributors/Copyright:
    2010 Oliver Borm (oli.borm@web.de)
    2011-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "patchMassFlowAverageFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
Field<T> patchMassFlowAverageFunctionObject::average(const word& fieldName,T unsetVal) const
{
    const VolField<T>& fld =
        obr_.lookupObject<VolField<T>>
        (
            fieldName
        );

    Field<T> vals(patchNames_.size(), unsetVal);

    const surfaceScalarField &phi=obr_.lookupObject<surfaceScalarField>(solver_.phi());

    forAll(patchNames_, patchI)
    {
        if (patchIndizes_[patchI] >= 0)
        {
            label index=patchIndizes_[patchI];
            scalar flux=sum(phi.boundaryField()[index]);
            reduce(flux,sumOp<scalar>());
            vals[patchI] = sum
                        (
                            phi.boundaryField()[index]
                            *fld.boundaryField()[index]
                        );
            reduce(vals[patchI],sumOp<T>());
            vals[patchI]/=flux;
        }
    }

    if (verbose()) {
        Info<< regionString()
            << " Mass-Flow-Weighted Averages of " << fieldName << " :";

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
