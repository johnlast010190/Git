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

#include "include/swak.H"

#ifdef FOAM_PATCHFIELDTYPE_IN_GEOFIELD_IS_NOW_PATCH
#define PatchFieldType Patch
#endif

#ifdef FOAM_HAS_ENERGY_HE

#include "recalcThermoHeFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "basicThermo/basicThermo.H"

#include "derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(recalcThermoHeFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        recalcThermoHeFunctionObject,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

recalcThermoHeFunctionObject::recalcThermoHeFunctionObject
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

void recalcThermoHeFunctionObject::recalc()
{
    //TODO: NO support for material library?
    basicThermo& thermo =
        const_cast<basicThermo&>
        (
            obr_.lookupObject<basicThermo>("thermophysicalProperties")
        );
    Info<< "Recalculating enthalpy h" << endl;

    const volScalarField& T = thermo.T();
    const volScalarField& p = thermo.p();
    volScalarField& h = thermo.he();

    labelList allCells(T.size());
    forAll(allCells, cellI)
    {
        allCells[cellI] = cellI;
    }
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
    const_cast<scalarField&>(h.internalField().field())
#else
    h.internalField()
#endif
    = thermo.he(p.internalField(), T.internalField(), allCells);
    forAll(h.boundaryField(), patchi)
    {
        const_cast<volScalarField::PatchFieldType&>
        (
            h.boundaryField()[patchi]
        ).forceAssign
        (
            thermo.he(T.boundaryField()[patchi], patchi)
        );
    }

    // hBoundaryCorrection
    volScalarField::Boundary& hbf =
        const_cast<volScalarField::Boundary&>(h.boundaryField());

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient() =
                hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad() =
                hbf[patchi].fvPatchField::snGrad();
        }
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

#endif

// ************************************************************************* //
