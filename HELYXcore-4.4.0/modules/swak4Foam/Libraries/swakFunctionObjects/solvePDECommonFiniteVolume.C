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
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "solvePDECommonFiniteVolume.H"

#include "meshes/polyMesh/polyMesh.H"

#include "fields/volFields/volFields.H"

#include "FieldValueExpressionDriver.H"

#include "fvMatrices/fvScalarMatrix/fvScalarMatrix.H"

#include "finiteVolume/fvm/fvm.H"

namespace Foam {
    defineTypeNameAndDebug(solvePDECommonFiniteVolume,0);
}

Foam::solvePDECommonFiniteVolume::solvePDECommonFiniteVolume
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
):
    solvePDECommon(
        name,
        obr,
        dict,
        loadFromFiles
    )
{
    if (!isA<polyMesh>(obr))
    {
        active_=false;
        WarningIn("solvePDECommonFiniteVolume::solvePDECommonFiniteVolume")
            << "Not a polyMesh. Nothing I can do"
                << endl;
    }
}

Foam::solvePDECommonFiniteVolume::~solvePDECommonFiniteVolume()
{}

Foam::volScalarField &Foam::solvePDECommonFiniteVolume::theField()
{
    // either the field was created by someone else ... then it should be
    // in the registry. Or we created it.
    if (theField_.valid()) {
        return theField_();
    } else {
        return const_cast<volScalarField&>(
            obr_.lookupObject<volScalarField>(
                fieldName_
            )
        );
    }
}

void Foam::solvePDECommonFiniteVolume::read(const dictionary& dict)
{
    solvePDECommon::read(dict);

    if (active_) {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (
            theField_.valid()
            &&
            fieldName_!=theField_->name()
        ) {
            WarningIn("Foam::solvePDECommonFiniteVolume::read(const dictionary& dict)")
                << "Throwing out field " << theField_->name()
                    << " and loading " << fieldName_ << ". "
                    << "This might lead to unpredicatable behaviour" << endl;
            theField_.clear();
        }
        if (!theField_.valid()) {
            if (obr_.foundObject<volScalarField>(fieldName_)) {
                if (!dict.found("useFieldFromMemory")) {
                    FatalErrorIn("Foam::solvePDECommonFiniteVolume::read(const dictionary& dict)")
                        << "Field " << fieldName_ << " alread in memory. "
                            << "Set 'useFieldFromMemory true;' to use it or "
                            << "use different name"
                            << endl
                            << exit(FatalError);

                }
                bool useFieldFromMemory=readBool(
                    dict.lookup("useFieldFromMemory")
                );
                if (!useFieldFromMemory) {
                    FatalErrorIn("Foam::solvePDECommonFiniteVolume::read(const dictionary& dict)")
                        << "Field " << fieldName_ << " alread in memory. "
                            << "Use different name"
                            << endl
                            << exit(FatalError);
                }
            } else {
                theField_.set(
                    new volScalarField(
                        IOobject (
                            fieldName_,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh
                    )
                );
            }
        }

        driver_.set(
            new FieldValueExpressionDriver(
                mesh.time().timeName(),
                mesh.time(),
                mesh,
                false, // no caching. No need
                true,  // search fields in memory
                false  // don't look up files on disc
            )
        );

        driver_->readVariablesAndTables(dict);

        driver_->createWriterAndRead(name_+"_"+fieldName_+"_"+type());
    }
}

void Foam::solvePDECommonFiniteVolume::writeData()
{
    theField().write();

    driver_->tryWrite();
}

void Foam::solvePDECommonFiniteVolume::writeNewField()
{
    theField().write();
}

void Foam::solvePDECommonFiniteVolume::writeOldField()
{
    volScalarField temp(fieldName_+".presolve",theField());
    temp.write();
}

// ************************************************************************* //
