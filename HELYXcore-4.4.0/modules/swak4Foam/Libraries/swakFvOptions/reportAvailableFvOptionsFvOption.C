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
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "reportAvailableFvOptionsFvOption.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(reportAvailableFvOptionsFvOption, 0);

    addToRunTimeSelectionTable
    (
        option,
        reportAvailableFvOptionsFvOption,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fv::reportAvailableFvOptionsFvOption::alwaysApply() const
{
    return true;
}

Foam::label Foam::fv::reportAvailableFvOptionsFvOption::applyToField(const word &fieldName) const
{
    return 0;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::reportAvailableFvOptionsFvOption::reportAvailableFvOptionsFvOption
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr)
{}

Foam::fv::reportAvailableFvOptionsFvOption::~reportAvailableFvOptionsFvOption()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::reportAvailableFvOptionsFvOption::sourceFields(wordList& fieldNames)
{
    fieldNames = List<word>(1,"dummy");
}


void Foam::fv::reportAvailableFvOptionsFvOption::report(
    const string &message
) {
    Info<< name() << " called: " << message.c_str() << endl;
}

void Foam::fv::reportAvailableFvOptionsFvOption::correct(volVectorField& U)
{
    report("correct(volVectorField& "+U.name()+")");
}


void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    report(
        "addSup(fvMatrix<vector>& "+eqn.psi().name()+")"
    );
}

// #ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+rho.name()+
        ",fvMatrix<vector>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+alpha.name()+
        ",const volScalarField& "+rho.name()+
        ",fvMatrix<vector>& "+eqn.psi().name()+")"
    );
}
// #endif

void Foam::fv::reportAvailableFvOptionsFvOption::setValue
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    report(
        "setValue(fvMatrix<vector>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    report(
        "constrain(fvMatrix<vector>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::correct(volScalarField& U)
{
    report("correct(volScalarField& "+U.name()+")");
}


void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    report(
        "addSup(fvMatrix<scalar>& "+eqn.psi().name()+")"
    );
}

// #ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+rho.name()+
        ",fvMatrix<scalar>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+alpha.name()+
        ",const volScalarField& "+rho.name()+
        ",fvMatrix<scalar>& "+eqn.psi().name()+")"
    );
}
// #endif

void Foam::fv::reportAvailableFvOptionsFvOption::setValue
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    report(
        "setValue(fvMatrix<scalar>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    report(
        "constrain(fvMatrix<scalar>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::correct(volTensorField& U)
{
    report("correct(volTensorField& "+U.name()+")");
}


void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(fvMatrix<tensor>& "+eqn.psi().name()+")"
    );
}

// #ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+rho.name()+
        ",fvMatrix<tensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+alpha.name()+
        ",const volScalarField& "+rho.name()+
        ",fvMatrix<tensor>& "+eqn.psi().name()+")"
    );
}
// #endif

void Foam::fv::reportAvailableFvOptionsFvOption::setValue
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    report(
        "setValue(fvMatrix<tensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::constrain
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    report(
        "constrain(fvMatrix<tensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::correct(volSymmTensorField& U)
{
    report("correct(volSymmTensorField& "+U.name()+")");
}


void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(fvMatrix<symmTensor>& "+eqn.psi().name()+")"
    );
}

// #ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+rho.name()+
        ",fvMatrix<symmTensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+alpha.name()+
        ",const volScalarField& "+rho.name()+
        ",fvMatrix<symmTensor>& "+eqn.psi().name()+")"
    );
}
// #endif

void Foam::fv::reportAvailableFvOptionsFvOption::setValue
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    report(
        "setValue(fvMatrix<symmTensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::constrain
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    report(
        "constrain(fvMatrix<symmTensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::correct(volSphericalTensorField& U)
{
    report("correct(volSphericalTensorField& "+U.name()+")");
}


void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(fvMatrix<sphericalTensor>& "+eqn.psi().name()+")"
    );
}

// #ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+rho.name()+
        ",fvMatrix<sphericalTensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    report(
        "addSup(const volScalarField& "+alpha.name()+
        ",const volScalarField& "+rho.name()+
        ",fvMatrix<sphericalTensor>& "+eqn.psi().name()+")"
    );
}
// #endif


void Foam::fv::reportAvailableFvOptionsFvOption::setValue
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    report(
        "setValue(fvMatrix<sphericalTensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::constrain
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    report(
        "constrain(fvMatrix<sphericalTensor>& "+eqn.psi().name()+")"
    );
}

void Foam::fv::reportAvailableFvOptionsFvOption::makeRelative(
    surfaceScalarField& phi
) const
{
    const_cast<reportAvailableFvOptionsFvOption&>(*this).
            report(
                "makeRelative(surfaceScalarField& "+phi.name()+")"
            );
}

void Foam::fv::reportAvailableFvOptionsFvOption::makeRelative(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    const_cast<reportAvailableFvOptionsFvOption&>(*this).
        report(
            "makeRelative(const surfaceScalarField& "+rho.name()+",surfaceScalarField& "+phi.name()+")"
        );
}

void Foam::fv::reportAvailableFvOptionsFvOption::makeRelative(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    const_cast<reportAvailableFvOptionsFvOption&>(*this).
            report(
                "makeRelative(FieldField<fvsPatchField, scalar>&)"
            );
}


void Foam::fv::reportAvailableFvOptionsFvOption::makeAbsolute(
    surfaceScalarField& phi
) const
{
    const_cast<reportAvailableFvOptionsFvOption&>(*this).
            report("makeAbsolute(surfaceScalarField& "+phi.name()+")");
}

void Foam::fv::reportAvailableFvOptionsFvOption::makeAbsolute(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    const_cast<reportAvailableFvOptionsFvOption&>(*this).
        report(
            "makeAbsolute(const surfaceScalarField& "+rho.name()+",surfaceScalarField& "+phi.name()+")"
        );
}


// ************************************************************************* //
