/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "foamMeshes.H"

#include "postDict/postDictKeys.H"

#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"

#undef Log
#include "vtkLogger.h"

namespace Foam::functionObjects::runTimeVis
{

void FoamMeshes::updateDomainRangeForColourFields()
{
    domainRanges_.clear();

    if (colourFields_.hasAllFieldsMarker())
    {
        FoamFields allFields;
        for (const MeshAndFields& meshAndFields : meshes_)
        {
            allFields.merge(meshAndFields.listAllFields());
        }
        updateDomainRangeForList(allFields);
    }
    else
    {
        updateDomainRangeForList(colourFields_);
    }
}

void FoamMeshes::updateDomainRangeForList(const FoamFields &fields)
{
    for (const foamField& field : fields.getColorFields())
    {
        vtkLogF(INFO, "Updating range for %s", field.c_str());
        scalarMinMax domainRange;

        for (const MeshAndFields& meshAndFields : meshes_)
        {
            domainRange += meshAndFields.calculateDomainRangeForField(field);
        }

        reduce(domainRange, minMaxOp<scalar>());
        domainRanges_.set(field, domainRange);
    }
}

void FoamMeshes::updateDomainBounds()
{
    boundBox domainBoundBox;
    for (const MeshAndFields& meshAndFields : meshes_)
    {
        boundBox currentBounds = meshAndFields.getMesh().bounds();
        domainBoundBox.add(currentBounds);
    }
    domainBoundBox.reduce();
    domainBounds_ = domainBoundBox;
}

void FoamMeshes::addColourField(const foamField& colourField)
{
    colourFields_.addField(colourField);
}

const objectRegistry &FoamMeshes::getFieldObjectRegistry(
    const Time &time,
    const word &meshRegionName,
    const word &instanceName
)
{
    const objectRegistry & fieldsRegistry = time.lookupObject<fvMesh>(meshRegionName);
    if (!instanceName.empty())
    {
        if (!time.foundObject<solutionInstanceRegistry>(instanceName))
        {
            FatalError << "Instance " << instanceName << " was specified, but was not found in the instance registry."
                       << abort(FatalError);
        }
        const auto &solReg =
            time.lookupObject<solutionInstanceRegistry>(instanceName);

        const auto& meshNames = solReg.meshNames();
        forAll(meshNames, n)
        {
            if (meshNames[n] == meshRegionName)
            {
                const auto& regionNames = solReg.regionNames();
                if (fieldsRegistry.foundObject<objectRegistry>(regionNames[n]))
                {
                    return  fieldsRegistry.lookupObject<objectRegistry>(regionNames[n]);
                }
            }
        }
    }

    return fieldsRegistry;
}

// * * * * * * * * * * * * *       Constructors        * * * * * * * * * * * //

FoamMeshes::FoamMeshes(const Time& time, const Dictionaries &dictionaries)
{
    const word instanceName = dictionaries.getRtppDict().lookupOrDefault<word>(postDictKeys::INSTANCE_KEY, "");

    for (const word& region : time.toc())
    {
        if(time.foundObject<objectRegistry>(region))
        {
            addFoamMesh(time, region, instanceName);
        }
    }
}

void FoamMeshes::addFoamMesh(const Time &time, const word &meshRegionName, const word& instanceName)
{
    // This is similar to fvMeshFunctionObject::selectMesh
    word meshName;
    const objectRegistry* obr = &(time.lookupObject<objectRegistry>(meshRegionName));

    // The dictionary operatingPointDict only exists when obr behaves as an
    // operating point of the multi-point of HELYX-Adjoint.
    if(obr->foundObject<IOdictionary>(postDictKeys::OPERATING_POINT_SUBDICT_KEY))
    {
        meshName =
            obr->lookupObject<IOdictionary>(postDictKeys::OPERATING_POINT_SUBDICT_KEY).
                lookup<word>(operatingPointKeys::PARENT_MESH_NAME_KEY);
    }
    else
    {
        meshName = meshRegionName;
    }

    if (time.foundObject<fvMesh>(meshName))
    {
        const auto &currentMesh = time.lookupObject<fvMesh>(meshName);
        const objectRegistry &currentFieldRegistry = getFieldObjectRegistry(
            time,
            meshRegionName,
            instanceName
        );

        MeshAndFields meshAndFields(currentMesh, currentFieldRegistry);
        meshes_.set(meshName, meshAndFields);
    }
}

scalarMinMax FoamMeshes::getDomainRangeForField(const foamField &field) const
{
    return domainRanges_(field.lessAssociation(), Foam::scalarMinMax());
}

// ************************************************************************* //
} // End namespace Foam
