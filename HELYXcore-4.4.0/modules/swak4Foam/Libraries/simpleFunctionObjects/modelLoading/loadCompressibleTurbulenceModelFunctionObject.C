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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "loadCompressibleTurbulenceModelFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

#include "include/swakThermoTypes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(loadCompressibleTurbulenceModelFunctionObject, 0);

    addNamedToRunTimeSelectionTable
    (
        functionObject,
        loadCompressibleTurbulenceModelFunctionObject,
        dictionary,
        loadCompressibleTurbulenceModel
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

loadCompressibleTurbulenceModelFunctionObject::loadCompressibleTurbulenceModelFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    modelLoadingFunctionObject<compressible::turbulenceModel>(name,t,dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    autoPtr<compressible::turbulenceModel> loadCompressibleTurbulenceModelFunctionObject::initModel()
{
    return compressible::turbulenceModel::New(
        obr().lookupObject<volScalarField>(
            dict_.lookup("rhoName")
        ),
        obr().lookupObject<volVectorField>(
            dict_.lookup("UName")
        ),
        obr().lookupObject<surfaceScalarField>(
            dict_.lookup("phiName")
        ),
        obr().lookupObject<swakFluidThermoType>(
            dict_.lookup("thermoName")
        )
    );
}


} // namespace Foam

// ************************************************************************* //
