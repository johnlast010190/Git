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
    2008-2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "recalcPhiFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "surfaceMesh/surfaceMesh.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "cfdTools/general/adjustPhi/adjustPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(recalcPhiFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        recalcPhiFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

recalcPhiFunctionObject::recalcPhiFunctionObject
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


bool recalcPhiFunctionObject::start()
{
    UName_ = word(dict_.lookup("UName"));
    phiName_ = word(dict_.lookup("phiName"));
    pName_ = word(dict_.lookup("pName"));
    rhoName_ = dict_.lookupOrDefault<word>("rhoName","none");
    writeOldFields_ = dict_.lookup<bool>("writeOldFields");
    writeFields_ = dict_.lookup<bool>("writeFields");

    return updateSimpleFunctionObject::start();
}


void recalcPhiFunctionObject::recalc()
{
    Info<< "Calculating flux field " << phiName_
        << " for velocity " << UName_ << endl;

    const fvMesh &mesh=dynamicCast<const fvMesh&>(obr_);

    const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);
    volScalarField& p =
        mesh.lookupObjectRef<volScalarField>(pName_);

    surfaceScalarField& phi =
        mesh.lookupObjectRef<surfaceScalarField>(phiName_);

    if (writeOldFields_)
    {
        Info<< "Writing copy of old " << phiName_ << endl;
        surfaceScalarField oldPhi(phiName_+".old",phi);
        oldPhi.write();
    }

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0, 0, 0))
    {
        phi = fvc::interpolate(U) & mesh.Sf();
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0, 0, 0))
    {
        if (rhoName_=="none")
        {
            // force read
            rhoName_=word(dict_.lookup("rhoName"));
        }
        const volScalarField &rho=mesh.lookupObject<volScalarField>(rhoName_);
        phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
     }
     else
     {
        FatalErrorIn("recalcPhiFunctionObject::calcPhi()")
            << "Can't deal with a flux field " << phiName_
                << " with dimensions " << phi.dimensions()
                << endl
                << exit(FatalError);

    }
    adjustPhi(phi, U, p);

    if (writeFields_)
    {
        Info<< "Writing new value of " << phiName_ << endl;
        phi.write();
    }
}

} // namespace Foam

// ************************************************************************* //
