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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermalBaffle/derivedFvPatchFields/thermalBaffle/thermalBaffleFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "mappedPatches/mappedPolyPatch/mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF),
    owner_(false),
    baffle_(),
    dict_(dictionary::null),
    extrudeMeshPtr_()
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField
    (
        ptf,
        p,
        iF,
        mapper
    ),
    owner_(ptf.owner_),
    baffle_(),
    dict_(ptf.dict_),
    extrudeMeshPtr_()
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF, dict),
    owner_(false),
    baffle_(),
    dict_(dict),
    extrudeMeshPtr_()
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (mesh.name() == polyMesh::defaultRegion)
    {
        const word regionName = dict_.lookupOrDefault<word>("region", "none");

        const word baffleName("3DBaffle" + regionName);

        if
        (
            !mesh.time().foundObject<fvMesh>(regionName)
         && regionName != "none"
        )
        {
            if (extrudeMeshPtr_.empty())
            {
                initBaffleMesh();
            }

            baffle_.reset
            (
                new regionModels::thermalBaffle
                (
                    "thermalBaffle",
                    mesh,
                    dict
                )
            );
            owner_ = true;
            baffle_->rename(baffleName);
        }
    }
}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(ptf, iF),
    owner_(ptf.owner_),
    baffle_(),
    dict_(ptf.dict_),
    extrudeMeshPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void thermalBaffleFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void thermalBaffleFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void thermalBaffleFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);
}


void thermalBaffleFvPatchScalarField::initBaffleMesh()
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const wordList patchNames({"bottom", "top", "side"});

    const wordList patchTypes
    ({
        mappedWallPolyPatch::typeName,
        mappedWallPolyPatch::typeName,
        dict_.lookup<bool>("columnCells")
      ? emptyPolyPatch::typeName
      : directPolyPatch::typeName
    });

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const wordList coupleGroup({mpp.coupleGroup()});
    const wordList coupleGroupSlave
    ({
        coupleGroup[0](0, coupleGroup[0].find('_')) + "_slave"
    });

    // Initialise the bottom and top patches
    List<dictionary> patchDicts(3);
    forAll(patchDicts, patchi)
    {
        patchDicts[patchi].set("nFaces", 0);
        patchDicts[patchi].set("startFace", 0);
    }
    patchDicts[0].add("coupleGroup", coupleGroup[0]);
    patchDicts[0].add("inGroups", coupleGroup);
    patchDicts[0].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);
    patchDicts[1].add("coupleGroup", coupleGroupSlave[0]);
    patchDicts[1].add("inGroups", coupleGroupSlave);
    patchDicts[1].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);

    List<polyPatch*> patchPtrs(3);
    forAll(patchPtrs, patchi)
    {
        patchPtrs[patchi] =
            polyPatch::New
            (
                patchTypes[patchi],
                patchNames[patchi],
                patchDicts[patchi],
                patchi,
                mesh.boundaryMesh()
            ).ptr();
    }

    extrudeMeshPtr_.reset
    (
        new extrudePatchMesh
        (
            mesh,
            patch(),
            dict_,
            dict_.lookup("region"),
            patchPtrs
        )
    );

    if (extrudeMeshPtr_.empty())
    {
        WarningInFunction
            << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " patchMeshPtr not set."
            << endl;
    }
}


void thermalBaffleFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (owner_ && mesh.name() == polyMesh::defaultRegion)
    {
        baffle_->evolve();
    }

    turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs();
}


void thermalBaffleFvPatchScalarField::write(Ostream& os) const
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::write(os);

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (owner_ && mesh.name() == polyMesh::defaultRegion)
    {
        forAllConstIter(dictionary, dict_, iter)
        {
            os << *iter;
        }
   }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thermalBaffleFvPatchScalarField
);

// Backward compatibility
addSpecialNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    thermalBaffleFvPatchScalarField,
    compressible::thermalBaffleWallFunction,
    compressible__thermalBaffleWallFunction
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
