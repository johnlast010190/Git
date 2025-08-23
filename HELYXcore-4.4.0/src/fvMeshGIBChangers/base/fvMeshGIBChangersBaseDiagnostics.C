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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "base/fvMeshGIBChangersBase.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "sets/topoSets/faceSet.H"
#include "motionSmoother/motionSmoother.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangersBase::writeProblematicCells(label cI) const
{
    if (debugMode_)
    {
        Pout<< "In cell" << tab << cI << tab << "at location"
             << tab << mesh().C()[cI] << endl;

        const labelList& cellsI = mesh().cells()[cI];

        Pout<< "cell faces:" << endl;

        forAll(cellsI, fI)
        {
            const label& faceI = cellsI[fI];
            if (faceI < mesh().nInternalFaces())
            {
                Pout<< fI << tab << faceI << tab << mesh().Cf()[faceI] << endl;
            }
            else
            {
                label patchI = mesh().boundaryMesh().whichPatch(faceI);
                label lfI = faceI - mesh().boundaryMesh()[patchI].start();
                Pout<< fI << tab << faceI << tab
                     << patchI << tab << mesh().boundary()[patchI].name() << tab
                     << mesh().Cf().boundaryField()[patchI][lfI] << endl;
            }
        }

        OFstream str
        (
            mesh().time().path()/
            "cells_zeroOldVol_" + mesh().time().timeName() + ".obj"
        );

        meshTools::writeOBJ
        (
            str,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            labelList(1, cI)
        );
    }
}


void Foam::fvMeshGIBChangersBase::writeProblematicCells
(
    label cI,
    const surfaceScalarField& ssf
) const
{
    if (debugMode_)
    {
        Pout<< "In cell" << tab << cI << tab << "at location"
             << tab << mesh().C()[cI] << endl;

        const labelList& cellsI = mesh().cells()[cI];

        Pout<< "cell faces:" << endl;

        forAll(cellsI, fI)
        {
            const label& faceI = cellsI[fI];
            if (faceI < mesh().nInternalFaces())
            {
                Pout<< fI << tab << faceI << tab
                     << mesh().Cf()[faceI] << tab
                     << ssf[faceI] << endl;
            }
            else
            {
                label patchI = mesh().boundaryMesh().whichPatch(faceI);
                Pout<< fI << tab << faceI << tab << patchI << tab
                    << mesh().boundary()[patchI].name() << endl;
            }
        }

        OFstream str
        (
            mesh().time().path()/
            "cells_zeroOldVol_" + mesh().time().timeName() + ".obj"
        );

        meshTools::writeOBJ
        (
            str,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            labelList(1, cI)
        );
    }
}


void Foam::fvMeshGIBChangersBase::writeProblematicCells() const
{
    if (debugMode_)
    {
        DynamicList<label> dlV0(mesh().cells().size());
        DynamicList<label> dlV(mesh().cells().size());
        const scalarField& cV0 = mesh().V0();
        const scalarField& cV = mesh().V();
        forAll(mesh().cells(), cI)
        {
            if (cV0[cI]<=0)
            {
                dlV0.append(cI);
                Pout<< "cell " << cI << tab << "negative V0 found" << endl;
            }
            if (cV[cI]<=0)
            {
                dlV.append(cI);
                Pout<< "cell " << cI << tab << "negative V found" << endl;
            }
        }
        dlV0.shrink();
        dlV.shrink();

        const labelList lV0 = labelList(dlV0);
        const labelList lV = labelList(dlV);

        label sizeV0 = lV0.size();
        label sizeV = lV.size();

        Foam::reduce(
            std::tie(sizeV0, sizeV),
            ParallelOp<sumOp<label>, sumOp<label>>{},
            mesh().comm()
        );

        if (sizeV0 != 0)
        {
            OFstream str
            (
                mesh().time().path()/
                "probCells_V0_" + mesh().time().timeName() + ".obj"
            );

            meshTools::writeOBJ
            (
                str,
                mesh().cells(),
                mesh().faces(),
                mesh().points(),
                lV0
            );
        }
        if (sizeV != 0)
        {
            OFstream str
            (
                mesh().time().path()/
                "probCells_V_" + mesh().time().timeName() + ".obj"
            );

            meshTools::writeOBJ
            (
                str,
                mesh().cells(),
                mesh().faces(),
                mesh().points(),
                lV
            );
        }
    }
}


void Foam::fvMeshGIBChangersBase::writeProblematicCells
(
    const boolList& interPoints
) const
{
    if (debugMode_)
    {
        const pointField& basePoints = *basePoints_;
        //const vectorField& baseCc = *baseCC_;
        DynamicList<label> dliP(mesh().cells().size());
        forAll(mesh().cells(), cI)
        {
            const labelList& cP = mesh().cellPoints()[cI];
            bool foundUnsnapped = false;
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                if (!interPoints[cPI])
                {
                    foundUnsnapped = true;
                }
            }
            if (!foundUnsnapped)
            {
                dliP.append(cI);
            }
        }
        dliP.shrink();

        const labelList liP = labelList(dliP);

        label sizeiP = liP.size();
        reduce(sizeiP, sumOp<label>());

        if (sizeiP != 0)
        {
            OFstream str
            (
                mesh().time().path()/
                "probCells_snapPoints_" + mesh().time().timeName() + ".obj"
            );

            meshTools::writeOBJ
            (
                str,
                mesh().cells(),
                mesh().faces(),
                basePoints,
                liP
            );
        }
    }
}


void Foam::fvMeshGIBChangersBase::checkRegion(const labelList& cR) const
{
    if (debugMode_)
    {
        bool regFound = false;
        forAll(cR, cI)
        {
            if
            (
                (cR[cI]!=0) &&
                (regFound==false)
            )
            {
                regFound = true;
            }
        }
        reduceToMaster(regFound, orOp<bool>());
        if (!regFound)
        {
            WarningInFunction
                << "in fvMeshGIBChangersBase::modifyRegionLabels"
                << "(labelList& cellIndi): "
                << " different regions not found" << endl;
        }
    }
}


void Foam::fvMeshGIBChangersBase::writeRegionInterface
(
    word name,
    const boolList& regFace
) const
{
    DynamicList<label> dregFaceL(regFace.size());
    const pointField& basePoints = *basePoints_;
    forAll(regFace, fI)
    {
        if (regFace[fI])
        {
            dregFaceL.append(fI);
        }
    }
    dregFaceL.shrink();

    labelList regFaceL(dregFaceL);
    simpleVTKWriter
    (
        mesh().faces(),
        regFaceL,
        basePoints
    ).write(name + "Base_" + mesh().time().timeName() + ".vtk");

    simpleVTKWriter
    (
        mesh().faces(),
        regFaceL,
        mesh().points()
    ).write(name + mesh().time().timeName() + ".vtk");
}


void Foam::fvMeshGIBChangersBase::writeRegionInterface
(
    word name,
    const boolList& regFace,
    label cI
) const
{
    DynamicList<label> dregFaceL(regFace.size());
    const pointField& basePoints = *basePoints_;

    if (Pstream::myProcNo()==6)
    {
        const labelList& cFaces = mesh().cells()[cI];
        forAll(cFaces, cFI)
        {
            const label& cFacesI = cFaces[cFI];
            if (regFace[cFacesI])
            {
                dregFaceL.append(cFacesI);
            }
        }
    }

    dregFaceL.shrink();
    labelList regFaceL(dregFaceL);
    simpleVTKWriter
    (
        mesh().faces(),
        regFaceL,
        basePoints
    ).write(name + "Base_" + mesh().time().timeName() + ".vtk");
}


void Foam::fvMeshGIBChangersBase::qualityControl() const
{
    const dictionary& dmCoeffs = dynamicMeshCoeffs_.
        subDict("meshQualityControls");

    faceSet errorFaces
    (
        mesh(),
        "errorFaces",
        mesh().nFaces()
    );

    bool hasErrors = motionSmoother::checkMesh
    (
        false,  // report
        mesh(),
        dmCoeffs,
        errorFaces
    );
    if (hasErrors)
    {
        Info<< "Mesh has errors" << endl;
    }

    errorFaces.write();
}


// ************************************************************************* //
