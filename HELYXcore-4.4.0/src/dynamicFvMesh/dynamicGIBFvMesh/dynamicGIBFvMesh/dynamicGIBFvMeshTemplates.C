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
    (c) 2015-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/syncTools/syncTools.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, template<class> class PatchField, class GeoMesh>
void StoreOldTimeFields(const fvMesh& mesh)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());
        field.storeOldTimes();

    }
}


template<class Type>
void CorrectBCs(const fvMesh& mesh)
{
    HashTable<const VolField<Type>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<VolField<Type>>()
    );

    for
    (
        typename HashTable<const VolField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        VolField<Type>& field =
            const_cast<VolField<Type>&>
            (*fieldIter());

        if
        (
            (field.name() == "k")
         || (field.name() == "epsilon")
        )
        {
            forAll(field, celli)
            {
                if (field[celli] < pTraits<Type>::zero)
                {
                    field[celli] = pTraits<Type>::zero;
                }
            }
            field.correctBoundaryConditions();
        }
        else
        {
            field.correctBoundaryConditions();
        }
    }
}


template<class Type>
void CorrectProcessorBCs(const fvMesh& mesh)
{
    HashTable<const VolField<Type>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<VolField<Type>>()
    );

    for
    (
        typename HashTable<const VolField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        VolField<Type>& field =
            const_cast<VolField<Type>&>
            (*fieldIter());

        label nReq = Pstream::nRequests();

        forAll(field.boundaryField(), patchi)
        {
            fvPatchField<Type>& pf = field.boundaryFieldRef()[patchi];
            if (isA<processorFvPatchField<Type>>(pf))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
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

        forAll(field.boundaryField(), patchi)
        {
            fvPatchField<Type>& pf = field.boundaryFieldRef()[patchi];
            if (isA<processorFvPatchField<Type>>(pf))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void dynamicGIBFvMesh::shrinkPopFields
(
    const labelList& fPop,
    const surfaceScalarField& phiPop,
    const boolList& popIndi
)
{
    const label dbgCell = -1;
    const label dbgProc = -1;
    const string dbgField = "";

    const boolList& popUpC = this->popUpCells();
    boolList popUpCNbr;
    syncTools::swapBoundaryCellList(*this, popUpC, popUpCNbr);

    const labelList& cReg = cRegion();
    labelList cRegNbr;
    syncTools::swapBoundaryCellList(*this, cReg, cRegNbr);

    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        this->thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        // Have not synced boundary values yet
        List<Type> fieldNbr;
        syncTools::swapBoundaryCellList(*this, field.primitiveField(), fieldNbr);

        const scalarField& cV0 = this->V0();
        scalarField V0 = cV0;

        if
        (
            Pstream::myProcNo() == dbgProc
         && (!dbgField.size() || field.name() == dbgField)
        )
        {
            Pout<< "Before PopFields: " << field.name()
                << "[" << dbgCell << "] = " << field[dbgCell] << endl;
        }

        for (int iter=0; iter < popSubSteps_; iter++)
        {
            scalarField V0Start = V0;
            GeometricField<Type, PatchField, GeoMesh> fieldStart = field;
            // Midpoint rule
            for (label m = 0; m < 2; m++)
            {
                scalar dd = 1.0/popSubSteps_;
                if (m == 0)
                {
                    dd /= 2;
                }
                const label nCells = this->cells().size();
                Field<Type> df(nCells, pTraits<Type>::zero);
                Type ff = pTraits<Type>::zero;
                Field<Type> minff(nCells, pTraits<Type>::one*GREAT);
                Field<Type> maxff(nCells, -pTraits<Type>::one*GREAT);

                scalarField sV(nCells, 0.0);

                forAll(fPop, fI)
                {
                    const label fPopI = fPop[fI];
                    if (fPopI < this->nInternalFaces())
                    {
                        const label own = this->owner()[fPopI];
                        const label nei = this->neighbour()[fPopI];

                        ff = pTraits<Type>::zero;
                        scalar sVf = dd*phiPop[fPopI]*this->time().deltaTValue();

                        if (cRegion()[own] == cRegion()[nei])
                        {
                            if (popUpC[own] == 1 && popUpC[nei] != 1)
                            {
                                ff = field[nei];
                            }
                            else if (popUpC[nei] == 1 && popUpC[own] != 1)
                            {
                                ff = field[own];
                            }
                            else if (popUpC[nei] == 1 && popUpC[own] == 1)
                            {
                                ff = 0.5*(field[own] + field[nei]);
                            }
                            else if (popUpC[nei] != 1 && popUpC[own] != 1)
                            {
                                ff = 0.5*(field[own] + field[nei]);
                            }

                            df[own] += ff*sVf;
                            df[nei] -= ff*sVf;

                            if (boundPopValues_)
                            {
                                minff[own] = min(minff[own], ff);
                                maxff[own] = max(maxff[own], ff);
                                minff[nei] = min(minff[nei], ff);
                                maxff[nei] = max(maxff[nei], ff);
                            }
                        }
                        else
                        {
                            // Don't try to conserve across the regions
                            // We should perhaps do something to ensure global
                            // conservation in each region, but allowing it
                            // to flow between them seems to make no sense.
                            // Setting it to zero here would do this but
                            // locally would probably produce spikes.
                            df[own] += field[own]*sVf;
                            df[nei] -= field[nei]*sVf;

                            if (boundPopValues_)
                            {
                                minff[own] = min(minff[own], field[own]);
                                maxff[own] = max(maxff[own], field[own]);
                                minff[nei] = min(minff[nei], field[nei]);
                                maxff[nei] = max(maxff[nei], field[nei]);
                            }
                        }

                        if
                        (
                            field.name() == dbgField
                         && Pstream::myProcNo() == dbgProc
                         && (own == dbgCell || nei == dbgCell)
                        )
                        {
                            Pout<< own << " " << nei
                                << " face: " << popUpC[own]
                                << " " << popUpC[nei]
                                << " " << cRegion()[own]
                                << " " << cRegion()[nei]
                                << " " << field[own]
                                << " " << field[nei]
                                << " " << sVf << endl;
                        }

                        sV[own] += sVf;
                        sV[nei] -= sVf;
                    }
                    else
                    {
                        const label own = this->faceOwner()[fPopI];
                        const label patchi =
                            this->boundaryMesh().whichPatch(fPopI);

                        if (!isA<indirectPolyPatch>(this->boundary()[patchi].patch()))
                        {
                            const label lpfI =
                                fPopI - this->boundaryMesh()[patchi].start();
                            const label bFacei = fPopI - this->nInternalFaces();
                            const scalar sVf = dd*
                                (
                                    phiPop.boundaryField()[patchi][lpfI]
                                )*
                                this->time().deltaTValue();

                            if (this->boundary()[patchi].coupled())
                            {
                                ff = pTraits<Type>::zero;

                                if (cRegion()[own] == cRegNbr[bFacei])
                                {
                                    if (popUpC[own] == 1 && popUpCNbr[bFacei] != 1)
                                    {
                                        ff = fieldNbr[bFacei];
                                    }
                                    else if (popUpCNbr[bFacei] == 1 && popUpC[own] != 1)
                                    {
                                        ff = field[own];
                                    }
                                    else if (popUpCNbr[bFacei] != 1 && popUpC[own] != 1)
                                    {
                                        ff = 0.5*(field[own] + fieldNbr[bFacei]);
                                    }
                                    else if (popUpCNbr[bFacei] == 1 && popUpC[own] == 1)
                                    {
                                        ff = 0.5*(field[own] + fieldNbr[bFacei]);
                                    }

                                    if
                                    (
                                        field.name() == dbgField
                                     && Pstream::myProcNo() == dbgProc
                                     && own == dbgCell
                                    )
                                    {
                                        Pout<< own << " "
                                            << " pface: " << fPopI
                                            << " " << popUpC[own]
                                            << " " << popUpCNbr[bFacei]
                                            << " " << cRegion()[own]
                                            << " " << cRegNbr[bFacei]
                                            << " " << field[own]
                                            << " " << fieldNbr[bFacei]
                                            << " " << ff
                                            << " " << ff*sVf << endl;
                                    }

                                    df[own] += ff*sVf;

                                    if (boundPopValues_)
                                    {
                                        minff[own] = min(minff[own], ff);
                                        maxff[own] = max(maxff[own], ff);
                                    }
                                }
                            }
                            else
                            {
                                if (popUpC[own] == 1)
                                {
                                    if
                                    (
                                        field.name() == dbgField
                                     && Pstream::myProcNo() == dbgProc
                                     && own == dbgCell
                                    )
                                    {
                                        Pout<< own << " bface: "
                                            << popUpC[own] << endl;
                                    }

                                    ff = field[own];
                                    df[own] += ff*sVf;

                                    if (boundPopValues_)
                                    {
                                        minff[own] = min(minff[own], ff);
                                        maxff[own] = max(maxff[own], ff);
                                    }

                                    sV[own] += sVf;
                                }
                            }
                        }
                    }
                }
                if (m == 1)
                {
                    V0 = V0Start;
                    field.ref() = fieldStart.ref();
                }
                forAll(popIndi, celli)
                {
                    if (popIndi[celli])
                    {
                        scalar newVol = V0[celli] + sV[celli];

                        if (mag(newVol) > SMALL)
                        {
                            field[celli] = fieldStart[celli]*V0[celli] + df[celli];
                            if
                            (
                                field.name() == dbgField
                             && Pstream::myProcNo() == dbgProc
                             && celli == dbgCell
                            )
                            {
                                Pout<< celli << " " << field[celli]/newVol
                                    << " " << fieldStart[celli]
                                    << " " << V0[celli]
                                    << " " << newVol << endl;
                            }
                            V0[celli] = newVol;
                            field[celli] /= V0[celli];
                        }
                        else
                        {
                          field[celli] = pTraits<Type>::zero;
                          V0[celli] = newVol;
                        }

                        if (boundPopValues_)
                        {
                            // Clamp the value so that it can't exceed the min/max
                            // of face values. Can occur if overall change in volume small
                            // compared to sum of abs values of face volume sweeps - then we
                            // end up with a large imbalance of incoming/outgoing fluxes due to slight
                            // variations in face values.
                            // TODO: redistribute conservatively rather
                            // Failsafe in case no faces had any influence
                            if (minff[celli] <= maxff[celli])
                            {
                                field[celli] = max(minff[celli], field[celli]);
                                field[celli] = min(maxff[celli], field[celli]);
                                if
                                (
                                    field.name() == dbgField
                                 && Pstream::myProcNo() == dbgProc
                                 && celli == dbgCell
                                )
                                {
                                    Pout<< celli << " " << minff[celli]
                                        << " " << maxff[celli] << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
        if
        (
            debugMode_
         && (
                field.name() == "T"
             || field.name() == "rho"
            )
        )
        {
            forAll(field, i)
            {
                if (refCast<scalarField>(field)[i] < 0)
                {
                    Pout<< "Neg " << field.name()
                        << " after pop: " << i
                        << " " << field[i] << endl;
                }
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void dynamicGIBFvMesh::growPopFields
(
    const labelList& fPop,
    const surfaceScalarField& phiPop,
    const boolList& popIndi,
    const scalarField& oldVolgrow
)
{
    const label dbgCell = -1;
    const label dbgProc = -1;
    const string dbgField = "";

    const boolList& popUpC = this->popUpCells();
    boolList popUpCNbr;
    syncTools::swapBoundaryCellList(*this, popUpC, popUpCNbr);

    const labelList& cReg = cRegion();
    labelList cRegNbr;
    syncTools::swapBoundaryCellList(*this, cReg, cRegNbr);

    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        this->thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        if
        (
            Pstream::myProcNo() == dbgProc
         && (!dbgField.size() || field.name() == dbgField)
        )
        {
            Pout<< "Before PopFields2: " << field.name()
                << "[" << dbgCell << "] = " << field[dbgCell] << endl;
        }

        // Initialise the faces of the pop cells by sweeping values
        // from non-pop cells in the same region (should maybe use FaceCellWave?)
        boolList popUpToDo(popUpC);
        while (true)
        {
            bool somethingChanged = false;
            bool allDone = true;
            boolList popUpToDoNext(popUpToDo);
            boolList popUpToDoNbr;
            syncTools::swapBoundaryCellList(*this, popUpToDo, popUpToDoNbr);
            // Have not synced boundary values yet
            List<Type> fieldNbr;
            syncTools::swapBoundaryCellList(*this, field.primitiveField(), fieldNbr);
            forAll(popUpToDo, celli)
            {
                if (popUpToDo[celli])
                {
                    label count = 0;
                    Type val = pTraits<Type>::zero;
                    forAll(this->cells()[celli], cFacei)
                    {
                        const label facei = this->cells()[celli][cFacei];
                        if (facei < this->nInternalFaces())
                        {
                            const label own = this->owner()[facei];
                            const label nei = this->neighbour()[facei];
                            if
                            (
                                field.name() == dbgField
                             && celli == dbgCell
                             && Pstream::myProcNo() == dbgProc
                            )
                            {
                                Pout<< facei << " " << cReg[own]
                                    << " " << cReg[nei]
                                    << " " << popUpToDo[own]
                                    << " " << popUpToDo[nei] << endl;
                            }
                            if (own == celli && !popUpToDo[nei] && cReg[nei] == cReg[celli])
                            {
                                val += field[nei];
                                count++;
                            }
                            else if (nei == celli && !popUpToDo[own] && cReg[own] == cReg[celli])
                            {
                                val += field[own];
                                count++;
                            }
                        }
                        else
                        {
                            const label patchi = this->boundaryMesh().whichPatch(facei);
                            if (this->boundary()[patchi].coupled())
                            {
                                const label bFacei = facei - this->nInternalFaces();
                                if (!popUpToDoNbr[bFacei] && cRegNbr[bFacei] == cReg[celli])
                                {
                                    val += fieldNbr[bFacei];
                                    count++;
                                }
                            }
                        }
                    }
                    if
                    (
                        field.name() == dbgField
                     && celli == dbgCell
                     && Pstream::myProcNo() == dbgProc
                    )
                    {
                        Pout<< count << endl;
                    }
                    if (count)
                    {
                        val /= count;
                        field[celli] = val;
                        popUpToDoNext[celli] = false;
                        somethingChanged = true;
                    }
                    else
                    {
                        if (debugMode_)
                        {
                            Pout<< celli << tab
                                << this->C()[celli] << tab
                                << cReg[celli] << tab
                                << this->cells()[celli] << tab
                                << endl;
                         }

                        allDone = false;

                        // Patch fix for sudden pop cells (boundary or internal)
                        // without any neighbours to take information from
                        // This happens mostly in the solid region
                        // Now we do not care about the solution inside the
                        // solid region. This happens rarely
                        // Most probably wrong value will be assign.
                        // However mean value will be a "safe" value to avoid
                        // assigning a zero density in compressible cases
                        // Otherwise another information has to be imported from
                        // the user or derived with another method
                        field[celli] = val;
                        field[celli] = average(field.primitiveField());
                    }
                }
            }
            popUpToDo = popUpToDoNext;
            reduce(allDone, andOp<bool>());
            if (!allDone)
            {
                reduce(somethingChanged, orOp<bool>());
                if (!somethingChanged)
                {
                    WarningInFunction
                        << "Unable to propagate initial values "
                        << "to all pop cells" << nl << endl;
                    break;
                }
                else if (debugMode_)
                {
                    Pout<< "Need another sweep for pop cell initialisation" << endl;
                }
            }
            else
            {
                break;
            }
        }

        // Have not synced boundary values yet
        List<Type> fieldNbr;
        syncTools::swapBoundaryCellList(*this, field.primitiveField(), fieldNbr);

        if (field.name() == dbgField && Pstream::myProcNo() == dbgProc)
        {
            Pout<< field[dbgCell] << endl;
        }

        scalarField oldVolg = oldVolgrow;

        for (int iter=0; iter < popSubSteps_; iter++)
        {
            scalarField volStart = oldVolg;
            GeometricField<Type, PatchField, GeoMesh> fieldStart = field;
            // Midpoint rule
            for (label m = 0; m < 2; m++)
            {
                scalar dd = 1.0/popSubSteps_;
                if (m == 0)
                {
                    dd /= 2;
                }
                const label nCells = this->cells().size();
                Field<Type> df(nCells, pTraits<Type>::zero);
                Type ff = pTraits<Type>::zero;
                Field<Type> minff(nCells, pTraits<Type>::one*GREAT);
                Field<Type> maxff(nCells, -pTraits<Type>::one*GREAT);

                scalarField sV(nCells, 0.0);

                forAll(fPop, fI)
                {
                    const label fPopI = fPop[fI];

                    if (fPopI < this->nInternalFaces())
                    {
                        const label own = this->owner()[fPopI];
                        const label nei = this->neighbour()[fPopI];

                        ff = pTraits<Type>::zero;
                        scalar sVf = dd*phiPop[fPopI]*this->time().deltaTValue();

                        if (cRegion()[own] == cRegion()[nei])
                        {
                            if (popUpC[own] == 1 && popUpC[nei] != 1)
                            {
                                ff = field[nei];
                            }
                            else if (popUpC[nei] == 1 && popUpC[own] != 1)
                            {
                                ff = field[own];
                            }
                            else if (popUpC[nei] != 1 && popUpC[own] != 1)
                            {
                                ff = 0.5*(field[own] + field[nei]);
                            }
                            else if (popUpC[nei] == 1 && popUpC[own] == 1)
                            {
                                ff = 0.5*(field[own] + field[nei]);
                            }

                            if
                            (
                                field.name() == dbgField
                             && Pstream::myProcNo() == dbgProc
                             && (own == dbgCell || nei == dbgCell)
                            )
                            {
                                Pout<< own << " " << nei
                                    << " face: " << fPopI
                                    << " " << popUpC[own]
                                    << " " << popUpC[nei]
                                    << " " << cRegion()[own]
                                    << " " << cRegion()[nei]
                                    << " " << field[own]
                                    << " " << field[nei]
                                    << " " << ff
                                    << " " << ff*sVf << endl;
                            }

                            df[own] += ff*sVf;
                            df[nei] -= ff*sVf;

                            if (boundPopValues_)
                            {
                                minff[own] = min(minff[own], ff);
                                maxff[own] = max(maxff[own], ff);
                                minff[nei] = min(minff[nei], ff);
                                maxff[nei] = max(maxff[nei], ff);
                            }
                        }
                        else
                        {
                            // Don't try to conserve across the regions
                            // We should perhaps do something to ensure global
                            // conservation in each region, but allowing it
                            // to flow between them seems to make no sense.
                            // Setting it to zero here would do this but
                            // locally would probably produce spikes.
                            if
                            (
                                field.name() == dbgField
                             && Pstream::myProcNo() == dbgProc
                             && (own == dbgCell || nei == dbgCell)
                            )
                            {
                                Pout<< own << " " << nei
                                    << " face: " << fPopI
                                    << " " << popUpC[own]
                                    << " " << popUpC[nei]
                                    << " " << cRegion()[own]
                                    << " " << cRegion()[nei]
                                    << " " << field[own]
                                    << " " << field[nei]
                                    << " " << field[own]*sVf
                                    << " " << field[nei]*sVf << endl;
                            }

                            df[own] += field[own]*sVf;
                            df[nei] -= field[nei]*sVf;

                            if (boundPopValues_)
                            {
                                minff[own] = min(minff[own], field[own]);
                                maxff[own] = max(maxff[own], field[own]);
                                minff[nei] = min(minff[nei], field[nei]);
                                maxff[nei] = max(maxff[nei], field[nei]);
                            }
                        }

                        sV[own] += sVf;
                        sV[nei] -= sVf;
                    }
                    else
                    {
                        const label own = this->faceOwner()[fPopI];
                        label patchi = this->boundaryMesh().whichPatch(fPopI);

                        if (!isA<indirectPolyPatch>(this->boundary()[patchi].patch()))
                        {
                            const label lpfI = fPopI - this->boundaryMesh()[patchi].start();
                            const label bFacei = fPopI - this->nInternalFaces();
                            const scalar sVf = dd*
                                (
                                    phiPop.boundaryField()[patchi][lpfI]
                                )*
                                this->time().deltaTValue();

                            if (this->boundary()[patchi].coupled())
                            {
                                ff = pTraits<Type>::zero;

                                if (cRegion()[own] == cRegNbr[bFacei])
                                {
                                    if (popUpC[own] == 1 && popUpCNbr[bFacei] != 1)
                                    {
                                        ff = fieldNbr[bFacei];
                                    }
                                    else if (popUpCNbr[bFacei] == 1 && popUpC[own] != 1)
                                    {
                                        ff = field[own];
                                    }
                                    else if (popUpCNbr[bFacei] != 1 && popUpC[own] != 1)
                                    {
                                        ff = 0.5*(field[own] + fieldNbr[bFacei]);
                                    }
                                    else if (popUpCNbr[bFacei] == 1 && popUpC[own] == 1)
                                    {
                                        ff = 0.5*(field[own] + fieldNbr[bFacei]);
                                    }

                                    if
                                    (
                                        field.name() == dbgField
                                     && Pstream::myProcNo() == dbgProc
                                     && own == dbgCell
                                    )
                                    {
                                        Pout<< own << " "
                                            << " pface: " << fPopI
                                            << " " << popUpC[own]
                                            << " " << popUpCNbr[bFacei]
                                            << " " << cRegion()[own]
                                            << " " << cRegNbr[bFacei]
                                            << " " << field[own]
                                            << " " << fieldNbr[bFacei]
                                            << " " << ff
                                            << " " << ff*sVf << endl;
                                    }

                                    df[own] += ff*sVf;
                                    sV[own] += sVf;

                                    if (boundPopValues_)
                                    {
                                        minff[own] = min(minff[own], ff);
                                        maxff[own] = max(maxff[own], ff);
                                    }
                                }
                            }
                            else
                            {
                                if (popUpC[own] == 1)
                                {
                                    ff = field[own];
                                    df[own] += ff*sVf;
                                    if (boundPopValues_)
                                    {
                                        minff[own] = min(minff[own], ff);
                                        maxff[own] = max(maxff[own], ff);
                                    }
                                    sV[own] += sVf;
                                }
                            }
                        }
                    }
                }
                if (m == 1)
                {
                    oldVolg = volStart;
                    field.ref() = fieldStart.ref();
                }
                forAll(popIndi, celli)
                {
                    if (popIndi[celli])
                    {
                        scalar newVol = oldVolg[celli] + sV[celli];
                        scalar oldVol = oldVolg[celli];

                        if (mag(newVol) > SMALL)
                        {
                            field[celli] = (field[celli]*oldVol + df[celli])/newVol;
                        }
                        else
                        {
                            // Bad fix for now at proc
                            const labelList& cellFaces(this->cells()[celli]);
                            forAll(cellFaces, fI)
                            {
                                label facei = cellFaces[fI];
                                if (facei < this->nInternalFaces())
                                {
                                    const label own = this->owner()[facei];
                                    const label nei = this->neighbour()[facei];
                                    if (popUpC[own] != 1)
                                    {
                                        field[celli] = field[own];
                                    }
                                    if (popUpC[nei] != 1)
                                    {
                                        field[celli] = field[nei];
                                    }
                                }
                            }
                        }

                        if (boundPopValues_)
                        {
                            // Clamp the value so that it can't exceed the min/max
                            // of face values. Can occur if overall change in volume small
                            // compared to sum of abs values of face volume sweeps - then we
                            // end up with a large imbalance of incoming/outgoing fluxes due to slight
                            // variations in face values.
                            // TODO: redistribute conservatively rather
                            // Failsafe in case no faces had any influence
                            if (minff[celli] <= maxff[celli])
                            {
                                field[celli] = max(minff[celli], field[celli]);
                                field[celli] = min(maxff[celli], field[celli]);
                            }
                        }

                        if
                        (
                            field.name() == dbgField &&
                            Pstream::myProcNo() == dbgProc &&
                            celli == dbgCell
                        )
                        {
                            Pout<< oldVol << tab
                                << "newVol: "
                                << newVol << tab
                                << field[celli] << tab
                                << df[celli] << endl;
                        }

                        oldVolg[celli] = newVol;
                    }
                }
            }
        }

        // Deal with ordinary boundaries attached to pop
        // cells: Their current value will not be meaningful
        // since they have changed regions; therefore
        // we initialise with extrapolated values
        forAll(fPop, fI)
        {
            const label fPopI = fPop[fI];
            if (fPopI >= this->nInternalFaces())
            {
                const label own = this->faceOwner()[fPopI];
                if (popIndi[own])
                {
                    const label patchi = this->boundaryMesh().whichPatch(fPopI);
                    if (!isA<indirectPolyPatch>(this->boundary()[patchi].patch()))
                    {
                        const label lpfI = fPopI - this->boundaryMesh()[patchi].start();
                        // Copy and reassign the whole patch in order to
                        // respect any overriding of operator=
                        Field<Type> bf(field.boundaryField()[patchi]);
                        bf[lpfI] = field[own];
                        field.boundaryFieldRef()[patchi] = bf;
                    }
                }
            }
        }

        const labelList& fsc = fullSnapCells();
        // Have not synced final boundary values yet
        syncTools::swapBoundaryCellList
        (
            *this,
            field.primitiveField(),
            fieldNbr
        );
        forAll(fsc, celli)
        {
            label gcI = fsc[celli];
            if (popUpC[gcI] == 1)
            {
                field[gcI] = pTraits<Type>::zero;
                scalar counter = 0;
                forAll(this->cells()[gcI], cFacei)
                {
                    const label facei = this->cells()[gcI][cFacei];
                    if (facei < this->nInternalFaces())
                    {
                        const label own = this->owner()[facei];
                        const label nei = this->neighbour()[facei];
                        if (cReg[nei] == cReg[own])
                        {
                            if (own == gcI)
                            {
                                if
                                (
                                    field.name() == dbgField
                                 && Pstream::myProcNo() == dbgProc
                                 && gcI == dbgCell
                                )
                                {
                                    Pout<< nei << " " << field[nei] << endl;
                                }

                                field[gcI] += field[nei];
                                counter += 1;
                            }
                            else
                            {

                                if
                                (
                                    field.name() == dbgField
                                 && Pstream::myProcNo() == dbgProc
                                 && gcI == dbgCell
                                )
                                {
                                    Pout<< own << " " << field[own] << endl;
                                }

                                field[gcI] += field[own];
                                counter += 1;
                            }
                        }
                    }
                    else
                    {
                        const label patchi =
                            this->boundaryMesh().whichPatch(facei);
                        if (this->boundaryMesh()[patchi].coupled())
                        {
                            const label bfaceI = facei - this->nInternalFaces();
                            if (cReg[gcI] == cRegNbr[bfaceI])
                            {
                                field[gcI] += fieldNbr[bfaceI];
                                counter++;
                            }
                        }
                    }
                }
                if (!counter)
                {
                    if (cRegion()[gcI] == 0)
                    {
                        FatalErrorInFunction
                            << "Found no usable neighbours for full snap pop cell: "
                            << endl
                            << "label: " << gcI << endl
                            << "C(): " << this->C()[gcI] << endl
                            << "CRefion: " << cRegion()[gcI] << endl
                            << "CRefion0: " << cRegion0()[gcI]
                            << exit(FatalError);
                    }
                    else
                    {
                        if (debugMode_)
                        {
                            Pout<< "Found no usable neighbours for full snap pop cell: "
                                << endl
                                << "However the cell is inside inactive zone" << endl
                                << "Initialization with average" << endl
                                << "label: " << gcI << endl
                                << "C(): " << this->C()[gcI] << endl
                                << "CRefion: " << cRegion()[gcI] << endl
                                << "CRefion0: " << cRegion0()[gcI]
                                << endl;
                        }
                        field[gcI] = average(field.primitiveField());
                        counter = 1;
                    }
                }

                field[gcI] /= counter;

                if
                (
                    (field.name() == dbgField)
                 && (Pstream::myProcNo() == dbgProc)
                 && (gcI == dbgCell)
                )
                {
                    Pout<< "fsc: " << field[gcI] << endl;
                }
            }
        }
        if
        (
            debugMode_
         && (
                (field.name() == "T")
             || (field.name() == "rho")
            )
        )
        {
            forAll(field, i)
            {
                if (refCast<scalarField>(field)[i] < 0)
                {
                    Pout<< "Neg" << tab
                        << field.name() << tab
                        << "after pop2: "
                        << i << tab
                        << field[i]
                        << endl;
                }
            }
        }

        if
        (
            Pstream::myProcNo() == dbgProc
         && (!dbgField.size() || field.name() == dbgField)
        )
        {
            Pout<< "After PopFields2: " << field.name()
                << "[" << dbgCell << "] = " << field[dbgCell] << endl;
        }
    }
}


template<class Type>
void Foam::dynamicGIBFvMesh::writeScalarField
(
    word name,
    const List<Type>& iF,
    bool writeNow
) const
{
    volScalarField outField
    (
        IOobject
        (
            name,
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimless, 0),
        "calculated"
    );

    // forAll loop because some = operations are
    // not general for instance boolList->scalarField
    // outField.internalField() = iF;
    forAll(outField, celli)
    {
        outField[celli] = iF[celli];
    }

    if (time().outputTime()||writeNow)
    {
        outField.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// end namespace FOAM
}
// ************************************************************************* //
