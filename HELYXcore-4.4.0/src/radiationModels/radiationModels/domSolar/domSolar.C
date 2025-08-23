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
    (c) 2010-2024 Engys Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "radiationModels/domSolar/domSolar.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/Vector/lists/vectorList.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "boundaryRadiationProperties/boundaryRadiationProperties.H"
#include "primitives/strings/wordRes/wordRes.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "primitives/strings/lists/stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(domSolar, 0);
        addToRadiationRunTimeSelectionTables(domSolar);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiationModels::domSolar::initialise(const dictionary& dict)
{
    //find and mark solar load source patches
    patchSourceCoeffs_.setSize(mesh_.boundaryMesh().size());
    forAll(patchSourceCoeffs_, patchI)
    {
        patchSourceCoeffs_.set
        (
            patchI,
            new scalarField(mesh_.boundaryMesh()[patchI].size(), 0.0)
        );
    }

    List<wordRe> sourcePatches(dict.lookup("sourcePatches"));
    wordRes spm(sourcePatches);

    const fvBoundaryMesh& bm = mesh_.boundary();
    const boundaryRadiationProperties& boundaryRadiation
    (
        boundaryRadiationProperties::New(mesh_, TRef_.value())
    );

    forAll(bm, patchI)
    {
        if (spm.match(bm[patchI].name()))
        {
            if
            (
                !boundaryRadiation.radBoundaryProperties()[patchI].empty()
            )
            {
                patchSourceCoeffs_[patchI]
                    = boundaryRadiation.transmissivity(patchI, 0, true);
            }
            else
            {
                patchSourceCoeffs_[patchI] = 1.0;
            }
        }
    }

    // populate directional source fields
    PtrList<entry> sources
    (
        (dict.lookupEntryPtr("sources",false,false))->stream()
    );

    IenvPtr_.setSize(sources.size());

    // count sources
    label sourcesSize = 0;
    forAllIter(PtrList<entry>, sources, iter)
    {
        if (iter().isDict())
        {
            sourcesSize++;
        }
    }
    patchLocalSourceCoeffs_.setSize(sourcesSize);

    label index = 0;
    forAllIter(PtrList<entry>, sources, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const dictionary& envRadDict = iter().dict();

        const word& name = iter().keyword();

        // use local source patches if specified
        if (envRadDict.found("localSourcePatches"))
        {
            // init
            patchLocalSourceCoeffs_.set
            (
                index,
                new FieldField<Field, scalar>(mesh_.boundaryMesh().size())
            );
            forAll(patchLocalSourceCoeffs_[index], patchI)
            {
                patchLocalSourceCoeffs_[index].set
                (
                    patchI,
                    new scalarField(mesh_.boundaryMesh()[patchI].size(), 0.0)
                );
            }

            // local source patches for i-source (index)
            List<wordRe> localSourcePatches(envRadDict.lookup("localSourcePatches"));
            wordRes lspm(localSourcePatches);

            forAll(bm, patchI)
            {
                if (lspm.match(bm[patchI].name()))
                {
                    if
                    (
                        !boundaryRadiation.radBoundaryProperties()[patchI].empty()
                    )
                    {
                        patchLocalSourceCoeffs_[index][patchI]
                            = boundaryRadiation.transmissivity(patchI, 0, true);
                    }
                    else
                    {
                        patchLocalSourceCoeffs_[index][patchI] = 1.0;
                    }
                }
            }

            IenvPtr_.set
            (
                index,
                new radiantField
                (
                    word("I"+name),
                    mesh_,
                    envRadDict,
                    patchLocalSourceCoeffs_[index]
                )
            );
            index ++;
        }
        else
        {
            IenvPtr_.set
            (
                index++,
                new radiantField
                (
                    word("I"+name),
                    mesh_,
                    envRadDict,
                    patchSourceCoeffs_
                )
            );
        }
    } // end source iter

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::domSolar::domSolar
(
    const volScalarField& T,
    dimensionedScalar TRef
)
:
    radiationModel(typeName, T, TRef),
    IenvPtr_(),
    Qenv_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    patchSourceCoeffs_(0),
    patchLocalSourceCoeffs_(),
    primary_(true)
{
    initialise(coeffs_);
}


Foam::radiationModels::domSolar::domSolar
(
    const dictionary& dict,
    const volScalarField& T,
    dimensionedScalar TRef
)
:
    radiationModel(typeName, dict, T, TRef),
    IenvPtr_(),
    Qenv_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    patchSourceCoeffs_(0),
    patchLocalSourceCoeffs_(),
    primary_(true)
{
    initialise(coeffs_);
}


Foam::radiationModels::domSolar::domSolar
(
    const dictionary& dict,
    const volScalarField& T,
    const word radWallFieldName,
    dimensionedScalar TRef
)
:
    radiationModel("none", T, TRef),
    IenvPtr_(),
    Qenv_
    (
        IOobject
        (
            radWallFieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    patchSourceCoeffs_(0),
    patchLocalSourceCoeffs_(),
    primary_(false)
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::domSolar::~domSolar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::domSolar::calculate()
{
    // reset environment heat flux to zero
    forAll(Qenv_.boundaryField(), bI)
    {
        const fvPatch& p = mesh_.boundary()[bI];

        if (!p.coupled())
        {
            Qenv_.boundaryFieldRef()[bI] = 0;
        }
    }

    // update radiant fields
    forAll(IenvPtr_, iI)
    {
        IenvPtr_[iI].calculate(Qenv_);
    }

    if (primary_) //apply absorptivity to scale incident radiation flux
    {
        const boundaryRadiationProperties& boundaryRadiation
        (
            boundaryRadiationProperties::New(mesh_, TRef_.value())
        );

        forAll(Qenv_.boundaryField(), patchID)
        {
            if
            (
                !boundaryRadiation.radBoundaryProperties()[patchID].empty()
            )
            {
                //only single band support
                Qenv_.boundaryFieldRef()[patchID]
                    *= boundaryRadiation.absorptivity(patchID);
            }
        }
    }
}


bool Foam::radiationModels::domSolar::read()
{
    initialise(coeffs_);

    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}



Foam::tmp<Foam::volScalarField> Foam::radiationModels::domSolar::Rp() const
{
    //this model has no volumetric contribution
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar
        (
            dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
            0
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::domSolar::Ru() const
{
    //this model has no volumetric contribution
    return tmp<DimensionedField<Foam::scalar, Foam::volMesh>>
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
        )
    );
}

// ************************************************************************* //
