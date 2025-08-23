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
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "SLGThermo/SLGThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SLGThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SLGThermo::SLGThermo(const fvMesh& mesh, fluidThermo& thermo)
:
    regIOobject
    (
        IOobject
        (
            SLGThermo::typeName,
            mesh.polyMesh::instance(),
            mesh
        )
    ),
    thermo_(thermo),
    carrier_(nullptr),
    liquids_(nullptr),
    solids_(nullptr)
{
    Info<< "Creating component thermo properties:" << endl;

    if (isA<basicSpecieMixture>(thermo))
    {
        basicSpecieMixture& mcThermo =
            dynamic_cast<basicSpecieMixture&>(thermo);
        carrier_ = &mcThermo;
        if (basicThermo::dictName == basicThermo::matDictName)
        {
            const label nSpecies(mcThermo.Y().size());
            haModels_.resize(nSpecies);
            CpModels_.resize(nSpecies);
            kappaModels_.resize(nSpecies);
            hsModels_.resize(nSpecies);
            hfModels_.resize(nSpecies);
            muModels_.resize(nSpecies);
            const word phaseName(mcThermo.Y().first().group());
            forAll(mcThermo.Y(), speciei)
            {
                const word specieName(mcThermo.Y()[speciei].member());
                matScalarTable& speciesModels =
                    thermo.materials().sTable(phaseName, specieName);
                haModels_.set(speciei, speciesModels[haModel::typeName]);
                CpModels_.set(speciei, speciesModels[CpModel::typeName]);
                kappaModels_.set(speciei, speciesModels[kappaModel::typeName]);
                hsModels_.set(speciei, speciesModels[hsModel::typeName]);
                hfModels_.set(speciei, speciesModels[hfModel::typeName]);
                muModels_.set(speciei, speciesModels[muModel::typeName]);
            }
        }
        Info<< "    multi-component carrier - " << mcThermo.species().size()
            << " species" << endl;
    }
    else
    {
        Info<< "    single component carrier" << endl;
    }

    if (thermo.properties().found("liquids"))
    {
        liquids_ =
            liquidMixtureProperties::New
            (
                thermo.properties().subDict("liquids")
            );
        Info<< "    liquids - " << liquids_->components().size()
            << " components" << endl;
    }
    else
    {
        Info<< "    no liquid components" << endl;
    }

    if (thermo.properties().found("solids"))
    {
        solids_  =
            solidMixtureProperties::New
            (
                thermo.properties().subDict("solids")
            );
        Info<< "    solids - " << solids_->components().size()
            << " components" << endl;
    }
    else
    {
        Info<< "    no solid components" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SLGThermo::~SLGThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fluidThermo& Foam::SLGThermo::thermo() const
{
    return thermo_;
}


Foam::scalar Foam::SLGThermo::Cp(label speciei, scalar p, scalar T)
{
    if (basicThermo::dictName != basicThermo::matDictName)
    {
        return carrier().Cp(speciei, p, T);
    }
    return CpModels_[speciei].value(p, T);
}


Foam::scalar Foam::SLGThermo::ha(label speciei, scalar p, scalar T)
{
    if (basicThermo::dictName != basicThermo::matDictName)
    {
        return carrier().ha(speciei, p, T);
    }
    return haModels_[speciei].value(p, T);
}


Foam::scalar Foam::SLGThermo::kappa(label speciei, scalar p, scalar T)
{
    if (basicThermo::dictName != basicThermo::matDictName)
    {
        return carrier().kappa(speciei, p, T);
    }
    return kappaModels_[speciei].value(p, T);
}


Foam::scalar Foam::SLGThermo::hs(label speciei, scalar p, scalar T)
{
    if (basicThermo::dictName != basicThermo::matDictName)
    {
        return carrier().hs(speciei, p, T);
    }
    return hsModels_[speciei].value(p, T);
}


Foam::scalar Foam::SLGThermo::mu(label speciei, scalar p, scalar T)
{
    if (basicThermo::dictName != basicThermo::matDictName)
    {
        return carrier().mu(speciei, p, T);
    }
    return muModels_[speciei].value(p, T);
}


Foam::scalar Foam::SLGThermo::hf(label speciei)
{
    if (basicThermo::dictName != basicThermo::matDictName)
    {
        return carrier().hf(speciei);
    }
    return hfModels_[speciei][0];
}


const Foam::basicSpecieMixture& Foam::SLGThermo::carrier() const
{
    if (carrier_ == nullptr)
    {
        FatalErrorInFunction
            << "carrier requested, but object is not allocated"
            << abort(FatalError);
    }

    return *carrier_;
}


const Foam::liquidMixtureProperties& Foam::SLGThermo::liquids() const
{
    if (!liquids_.valid())
    {
        FatalErrorInFunction
            << "liquids requested, but object is not allocated"
            << abort(FatalError);
    }

    return liquids_();
}


const Foam::solidMixtureProperties& Foam::SLGThermo::solids() const
{
    if (!solids_.valid())
    {
        FatalErrorInFunction
            << "solids requested, but object is not allocated"
            << abort(FatalError);
    }

    return solids_();
}


Foam::label Foam::SLGThermo::carrierId
(
    const word& cmptName,
    bool allowNotfound
) const
{
    forAll(carrier().species(), i)
    {
        if (cmptName == carrier_->species()[i])
        {
            return i;
        }
    }

    if (!allowNotfound)
    {
        FatalErrorInFunction
            << "Unknown carrier component " << cmptName
            << ". Valid carrier components are:" << nl
            << carrier_->species() << exit(FatalError);
    }

    return -1;
}


Foam::label Foam::SLGThermo::liquidId
(
    const word& cmptName,
    bool allowNotfound
) const
{
    forAll(liquids().components(), i)
    {
        if (cmptName == liquids_->components()[i])
        {
            return i;
        }
    }

    if (!allowNotfound)
    {
        FatalErrorInFunction
            << "Unknown liquid component " << cmptName << ". Valid liquids are:"
            << nl << liquids_->components() << exit(FatalError);
    }

    return -1;
}


Foam::label Foam::SLGThermo::solidId
(
    const word& cmptName,
    bool allowNotfound
) const
{
    forAll(solids().components(), i)
    {
        if (cmptName == solids_->components()[i])
        {
            return i;
        }
    }

    if (!allowNotfound)
    {
        FatalErrorInFunction
            << "Unknown solid component " << cmptName << ". Valid solids are:"
            << nl << solids_->components() << exit(FatalError);
    }

    return -1;
}


bool Foam::SLGThermo::hasMulticomponentCarrier() const
{
    return (carrier_ != nullptr);
}


bool Foam::SLGThermo::hasLiquids() const
{
    return liquids_.valid();
}


bool Foam::SLGThermo::hasSolids() const
{
    return solids_.valid();
}


// ************************************************************************* //
