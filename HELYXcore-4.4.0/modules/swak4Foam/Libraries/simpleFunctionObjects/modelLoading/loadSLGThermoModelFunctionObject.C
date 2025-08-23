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

#include "include/swak.H"

#ifndef FOAM_NO_SLG_THERMOPHYSICS

#include "loadSLGThermoModelFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

#include "include/swakThermoTypes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    template <>
    void modelLoadingFunctionObject<SLGThermo>::writeSimple()
    {
        if (correctModel_) {
            if (model_.valid()) {
                FatalErrorIn("modelLoadingFunctionObject::start()")
                    << "SLGThermo has no correct method"
                        << endl
                        << exit(FatalError);
            } else {
                FatalErrorIn("modelLoadingFunctionObject::start()")
                    << "Model has never been intialized"
                        << endl
                        << exit(FatalError);
            }
        }
    }

    defineTypeNameAndDebug(loadSLGThermoModelFunctionObject, 0);

    addNamedToRunTimeSelectionTable
    (
        functionObject,
        loadSLGThermoModelFunctionObject,
        dictionary,
        loadSLGThermoModel
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

loadSLGThermoModelFunctionObject::loadSLGThermoModelFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    modelLoadingFunctionObject<SLGThermo>(name,t,dict)
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<SLGThermo> loadSLGThermoModelFunctionObject::initModel()
{
    autoPtr<SLGThermo> result(
        new SLGThermo(
            dynamicCast<const fvMesh &>(
                obr()
            ),
            const_cast<swakFluidThermoType &>(
                obr().lookupObject<swakFluidThermoType>(
                    dict_.lookup("thermoName")
                )
            )
        )
    );

    return result;
}


} // namespace Foam

#endif

// ************************************************************************* //
