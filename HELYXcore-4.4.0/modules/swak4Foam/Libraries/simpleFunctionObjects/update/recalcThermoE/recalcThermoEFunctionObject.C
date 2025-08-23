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
    2008-2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "foamVersion4swak.H"

#if FOAM_VERSION4SWAK_MAJOR<2 || (FOAM_VERSION4SWAK_MAJOR==2 && FOAM_VERSION4SWAK_MINOR<2)

#include "recalcThermoEFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "basicThermo/basicThermo.H"

#include "fixedInternalEnergyFvPatchScalarField.H"
#include "gradientInternalEnergyFvPatchScalarField.H"
#include "mixedInternalEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(recalcThermoEFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        recalcThermoEFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

recalcThermoEFunctionObject::recalcThermoEFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    updateSimpleFunctionObject(name,t,dict)
{
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}

void recalcThermoEFunctionObject::recalc()
{
    basicThermo &thermo=const_cast<basicThermo&>(
        obr_.lookupObject<basicThermo>("thermophysicalProperties")
    );
    Info<< "Recalculating internal energy e" << endl;

    const volScalarField &T=thermo.T();
    volScalarField &e=thermo.e();

    labelList allCells(T.size());
    forAll(allCells,cellI) {
        allCells[cellI]=cellI;
    }
    e.internalField()=thermo.e(
        T.internalField(),
        allCells
    );
    forAll(e.boundaryField(), patchi)
    {
        e.boundaryField()[patchi] ==
            thermo.e(
                T.boundaryField()[patchi],
                patchi
            );
    }

    // eBoundaryCorrection
    volScalarField::GeometricBoundaryField& ebf = e.boundaryField();

    forAll(ebf, patchi)
    {
        if (isA<gradientInternalEnergyFvPatchScalarField>(ebf[patchi]))
        {
            refCast<gradientInternalEnergyFvPatchScalarField>(ebf[patchi]).gradient()
                = ebf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedInternalEnergyFvPatchScalarField>(ebf[patchi]))
        {
            refCast<mixedInternalEnergyFvPatchScalarField>(ebf[patchi]).refGrad()
                = ebf[patchi].fvPatchField::snGrad();
        }
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

#endif

// ************************************************************************* //
