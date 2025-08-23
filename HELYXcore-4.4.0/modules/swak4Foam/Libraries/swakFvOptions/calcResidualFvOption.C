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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2014-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "calcResidualFvOption.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(calcResidualFvOption, 0);

    addToRunTimeSelectionTable
    (
        option,
        calcResidualFvOption,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::calcResidualFvOption::calcResidualFvOption
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr),
    verbose_(coeffs_.lookup<bool>("verbose")),
    doAddSup_(coeffs_.lookup<bool>("doAddSup")),
    doSetValue_(coeffs_.lookup<bool>("doSetValue"))
{}

Foam::fv::calcResidualFvOption::~calcResidualFvOption()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::calcResidualFvOption::sourceFields(wordList& fieldNames)
{
    fieldNames = coeffs_.lookup<wordList>("fieldNames");
}


template<class Type>
void Foam::fv::calcResidualFvOption::calcResidual(
    const word &prefix,
    const fvMatrix<Type> &mat
) {
    const word &fName=mat.psi().name();

    if (!mat.hasDiag()) {
        WarningIn("Foam::fv::calcResidualFvOption::calcResidual")
            << "Matrix for " << fName << " has no diagonal"
                << endl;
        return;
    }
    if (!mat.hasLower() && !mat.hasUpper()) {
        WarningIn("Foam::fv::calcResidualFvOption::calcResidual")
            << "Matrix for " << fName << " has no upper and no lower part"
                << endl;
        return;
    }
    tmp<Field<Type>> pRes(mat.residual());
    Field<Type> &res=const_cast<Field<Type>&>(pRes());

    scalarField tmpField(res.size());

    Field<Type> mult((mat & mat.psi())().internalField());

    Type normFactor=pTraits<Type>::one;
    // for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++) {
    //     setComponent(
    //         normFactor,
    //         cmpt
    //     ) =
    //         gSum(
    //             mag(mult.component(cmpt))
    //             +
    //             mag(mat.source().component(cmpt))
    //         )
    //         +
    //         SolverPerformance<Type>::small_;
    // }

    Type residual=cmptDivide(gSumCmptMag(res), normFactor);
    Type residual2=cmptDivide(gSumCmptMag(mult), normFactor);

    if (verbose_) {
        Info<< name() << ": "
            << prefix << "-Residual for " << fName << ": " << residual
            << " " << residual2
            << "Norm factor: " << normFactor
            << endl;
    }
}


void Foam::fv::calcResidualFvOption::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (doAddSup_) {
        calcResidual("addSup",eqn);
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP

void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::calcResidualFvOption::setValue
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (doSetValue_) {
        calcResidual("setValue",eqn);
    }
}

void Foam::fv::calcResidualFvOption::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (doAddSup_) {
        calcResidual("addSup",eqn);
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::calcResidualFvOption::setValue
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (doSetValue_) {
        calcResidual("setValue",eqn);
    }
}


void Foam::fv::calcResidualFvOption::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    if (doAddSup_) {
        calcResidual("addSup",eqn);
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::calcResidualFvOption::setValue
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    if (doSetValue_) {
        calcResidual("setValue",eqn);
    }
}

void Foam::fv::calcResidualFvOption::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    if (doAddSup_) {
        calcResidual("addSup",eqn);
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::calcResidualFvOption::setValue
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    if (doSetValue_) {
        calcResidual("setValue",eqn);
    }
}


void Foam::fv::calcResidualFvOption::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    if (doAddSup_) {
        calcResidual("addSup",eqn);
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::calcResidualFvOption::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::calcResidualFvOption::setValue
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    if (doSetValue_) {
        calcResidual("setValue",eqn);
    }
}

// ************************************************************************* //
