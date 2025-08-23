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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "explicitPorositySource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(explicitPorositySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        explicitPorositySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::explicitPorositySource::explicitPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    porosityPtr_(nullptr),
    adjointVelocityPrefix_("Ua_"),
    UName_("U")
{
    read(dict);

    if (selectionMode_ != smCellZone)
    {
        FatalErrorInFunction
            << "selection mode is " << selectionModeTypeNames_[selectionMode_]
            << exit(FatalError);
    }

    porosityPtr_.reset
    (
        porosityModel::New
        (
            name_,
            obr_,
            mesh_,
            coeffs_,
            cellSetName_
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::explicitPorositySource::sourceFields
(
    wordList& fieldNames
)
{
    if (coeffs_.found("UNames"))
    {
        fieldNames = coeffs_.lookup<wordList>("UNames");
    }
    else if (coeffs_.found("U"))
    {
        UName_ = coeffs_.lookup<word>("U");
        fieldNames = wordList(1, UName_);
    }
    else
    {
        fieldNames = wordList(1, UName_);
    }
}


void Foam::fv::explicitPorositySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (fieldNames()[fieldi].find(adjointVelocityPrefix_, 0) == string::npos)
    {
        porosityPtr_->addResistance(porosityEqn);
    }
    else
    {
        const volVectorField& Uprimal =
            obr_.lookupObject<volVectorField>(UName_);

        porosityPtr_->addAdjointResistance(porosityEqn, Uprimal);
    }

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    //- Here Assuming the dimensions are the same inside the matrix which is
    //  true for the porosity classes. They are applied on the momentum block
    //  only. It is done to avoid coding more constructors in block
    //  Code can be easily extended if needed in the future
    fvBlockMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensionSets()[0]);

    porosityPtr_->addResistance(porosityEqn);

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (fieldNames()[fieldi].find(adjointVelocityPrefix_, 0) == string::npos)
    {
        porosityPtr_->addResistance(porosityEqn);
    }
    else
    {
        const volVectorField& Uprimal =
            obr_.lookupObject<volVectorField>(UName_);

        porosityPtr_->addAdjointResistance(porosityEqn, Uprimal);
    }

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    //- Here Assuming the dimensions are the same inside the matrix which is
    //  true for the porosity classes. They are applied on the momentum block
    //  only. It is done to avoid coding more constructors in block
    //  Code can be easily extended if needed in the future
    fvBlockMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensionSets()[0]);

    porosityPtr_->addResistance(porosityEqn);

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (fieldNames()[fieldi].find(adjointVelocityPrefix_, 0) == string::npos)
    {
        porosityPtr_->addResistance(porosityEqn);
    }
    else
    {
        const volVectorField& Uprimal =
            obr_.lookupObject<volVectorField>(UName_);

        porosityPtr_->addAdjointResistance(porosityEqn, Uprimal);
    }

    eqn -= alpha*porosityEqn;
}


Foam::label Foam::fv::explicitPorositySource::applyToField
(
    const word& fieldName, const word& regionName
) const
{
    //inject adjoint velocity fields into lists if they are called

    if
    (
        fieldName.find(adjointVelocityPrefix_, 0) != string::npos &&
        findIndex(fieldNames(), fieldName) == -1
    )
    {
        wordList& cfieldNames = const_cast<wordList&>(fieldNames());
        wordList& cregionNames = const_cast<wordList&>(regionNames());
        boolList& capplied = const_cast<boolList&>(applied());
        label fsize = fieldNames().size();

        cfieldNames.setSize(fsize+1);
        cfieldNames[fsize] = fieldName;

        if (cregionNames.size())
        {
            cregionNames.setSize(fsize+1);
            cregionNames[fsize] =
                (obr_.name() == word::null ? mesh_.name() : obr_.name());
        }

        capplied.setSize(fsize+1);
        capplied[fsize] = true;
    }

    return option::applyToField(fieldName, regionName);
}

// ************************************************************************* //
