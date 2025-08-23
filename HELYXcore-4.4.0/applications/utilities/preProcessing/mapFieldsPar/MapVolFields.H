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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef MapConsistentVolFields_H
#define MapConsistentVolFields_H

#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "meshToMesh/meshToMesh.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void evaluateConstraintTypes(VolField<Type>& fld)
{
    typename VolField<Type>::
        Boundary& fldBf = fld.boundaryFieldRef();

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(Pstream::defaultCommsType);
            }
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(Pstream::defaultCommsType);
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (patchSchedule[patchEvali].init)
                {
                    tgtField.initEvaluate(Pstream::commsTypes::scheduled);
                }
                else
                {
                    tgtField.evaluate(Pstream::commsTypes::scheduled);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
void MapVolFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    const meshToMesh& interp,
    const CombineOp& cop
)
{
    typedef VolField<Type> fieldType;

    const fvMesh& meshSource = static_cast<const fvMesh&>(interp.srcRegion());
    const fvMesh& meshTarget = static_cast<const fvMesh&>(interp.tgtRegion());

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            const fieldType fieldSource(*fieldIter(), meshSource);

            IOobject targetIO
            (
                fieldName,
                meshTarget.time().timeName(),
                meshTarget,
                IOobject::MUST_READ
            );

            if (targetIO.typeHeaderOk<fieldType>(true))
            {
                Info<< "    interpolating onto existing field "
                    << fieldName << endl;
                fieldType fieldTarget(targetIO, meshTarget);

                interp.mapSrcToTgt(fieldSource, cop, fieldTarget);

                evaluateConstraintTypes(fieldTarget);

                fieldTarget.write();
            }
            else
            {
                Info<< "    creating new field "
                    << fieldName << endl;

                targetIO.readOpt() = IOobject::NO_READ;

                tmp<fieldType>
                    tfieldTarget(interp.mapSrcToTgt(fieldSource, cop));

                fieldType fieldTarget(targetIO, tfieldTarget);

                evaluateConstraintTypes(fieldTarget);

                fieldTarget.write();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
