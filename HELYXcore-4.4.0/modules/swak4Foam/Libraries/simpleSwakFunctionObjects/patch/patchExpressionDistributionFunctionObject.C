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
    2008-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "patchExpressionDistributionFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"

#include "primitives/strings/lists/wordReList.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/strings/lists/stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchExpressionDistributionFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchExpressionDistributionFunctionObject,
        dictionary
    );



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchExpressionDistributionFunctionObject::patchExpressionDistributionFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    distributionFunctionObject(name,t,dict),
    expression_(
        dict.lookup("expression"),
        dict
    ),
    data_(dict)
{
    const fvMesh &mesh=refCast<const fvMesh>(obr_);

    wordReList newPatches(dict.lookup("patchNames"));
    HashSet<word> patchNamesNew;
    wordList allPatches(mesh.boundaryMesh().names());
    forAll(newPatches,i) {
        labelList IDs=findStrings(newPatches[i],allPatches);
        forAll(IDs,j) {
            patchNamesNew.insert(allPatches[IDs[j]]);
        }
    }

    patchNames_=patchNamesNew.toc();
    patchIDs_.setSize(patchNames_.size());

    forAll(patchNames_,i) {
        const word &name=patchNames_[i];
        patchIDs_[i]=mesh.boundaryMesh().findPatchID(name);
        if (patchIDs_[i]<0) {
            FatalErrorIn("patchExpressionDistributionFunctionObject::patchExpressionDistributionFunctionObject")
                << "Patch name " << name << " is not a valid patch"
                    << endl
                    << exit(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word patchExpressionDistributionFunctionObject::dirName()
{
    word result(typeName);
    forAll(patchNames_,i) {
    result+="_"+patchNames_[i];
    }

    return result;
}

word patchExpressionDistributionFunctionObject::baseName()
{
    word result("expression");
    forAll(patchNames_,i) {
    result+="_"+patchNames_[i];
    }

    return result;
}

void patchExpressionDistributionFunctionObject::getDistribution()
{
    if (drivers_.size()<=0) {
        WarningIn("patchExpressionDistributionFunctionObject::getDistribution")
            << "No drivers/patches specified"
                << endl;
        return;
    }

    getDistributionInternal(distScalar_);
    getDistributionInternal(distVector_);
    getDistributionInternal(distTensor_);
    getDistributionInternal(distSymmTensor_);
    getDistributionInternal(distSphericalTensor_);
}

void patchExpressionDistributionFunctionObject::writeSimple()
{
    forAll(patchIDs_,i) {
        PatchValueExpressionDriver &driver=drivers_[i];

        if (verbose()) {
            Info<< "DistributionExpression " << name() << " on " << patchNames_[i] << ": ";
        }

        driver.clearVariables();
        driver.parse(expression_);
    }

    distributionFunctionObject::writeSimple();

    if (verbose()) {
        Info<< endl;
    }

    // make sure that the stored Variables are consistently written
    forAll(drivers_,i) {
        drivers_[i].tryWrite();
    }
}

bool patchExpressionDistributionFunctionObject::start()
{
    bool result=distributionFunctionObject::start();
    const fvMesh &mesh=refCast<const fvMesh>(obr_);


    drivers_.clear();
    drivers_.resize(patchIDs_.size());

    forAll(drivers_,i) {
        drivers_.set(
            i,
            new PatchValueExpressionDriver(
                data_,
                mesh.boundary()[patchIDs_[i]],
                obr_
            )
        );
        drivers_[i].createWriterAndRead(
            name()+"_"+mesh.boundary()[patchIDs_[i]].name()+"_"+type()
        );
    }
    return result;
}

} // namespace Foam

// ************************************************************************* //
