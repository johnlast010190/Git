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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "db/error/messageStream.H"
#include "rtppActor.H"
#include "Utils/Utils.H"
#include "colourMaps/colours.H"
#include "postDict/postDictKeys.H"
#include "rendering/renderManager.H"

// VTK includes
#include "vtkActor.h"
#include "vtkFeatureEdges.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkAlgorithmOutput.h"
#include "vtkScalarsToColors.h"
#include "engysOutlineFilter.h"
#include "vtkMultiProcessController.h"
#include "vtkBitArray.h"
#include "vtkIncrementalOctreePointLocator.h"
#include "vtkCompositeSurfaceLICMapper.h"
#include "vtkSurfaceLICInterface.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellCenters.h"
#include "vtkMaskPoints.h"
#include "vtkArrowSource.h"
#include "vtkGlyph3D.h"

#define DEFAULT_AMBIENT 0.0
#define DEFAULT_DIFFUSE 0.9
#define DEFAULT_SPECULAR 0.1
#define DEFAULT_SPECULAR_POWER 0.8

static Foam::label coincidentResolutionFactor = 0;

namespace Foam::functionObjects::runTimeVis
{

rtppActor::rtppActor(const Visualisation &visualisation, int index)
    : visualisation_(visualisation), index(index)
{
    initialiseActor();
    initialiseBackfaceStyling();
    initialiseMapper();
    initialiseRepresentation();
    initialiseActiveMapper();
    initialiseNormalsActor();
}

void
rtppActor::processInputData(vtkDataSet *input)
{
    vtkSmartPointer<vtkPolyData> processed;

    vtkNew<vtkFieldData> inputFieldData;
    inputFieldData->DeepCopy(input->GetFieldData());

    if (normalsFilter_)
    {
        normalsFilter_->SetInputData(nullptr);
    }

    if (isOutline())
    {
        FatalError << "Outline data should be set with the setProcessedInputData method" << endl << abort(FatalError);
    }
    else
    {
        vtkSmartPointer<vtkPolyData> surfaceInput = Utils::getSurfacePolyDataFrom(input);
        if (visualisation_.representation.isProfile())
        {
            processed = Utils::mergePolyData(surfaceInput);
            edgeFilter_->SetInputData(processed);
            edgeFilter_->Update();
            processed = edgeFilter_->GetOutput();
        }
        else
        {
            processed = surfaceInput;
            processed->UpdateCellGhostArrayCache();
            processed->RemoveGhostCells();
            if (normalsFilter_)
            {
                normalsFilter_->SetInputData(processed);
            }
        }
    }

    processed->SetFieldData(inputFieldData);
    setProcessedInputData(processed);
}

void rtppActor::setProcessedInputData(vtkPolyData *input)
{
    processedInput_ = input;
    if (visualisation_.showActiveGIBBoundary)
    {
        updateActiveMapperAndRepresentation(input);
    }
    else if (mapper_)
    {
        mapper_->SetInputData(processedInput_);
    }
}

vtkSmartPointer<vtkPolyData> rtppActor::getProcessedInputData()
{
    return processedInput_;
}

void rtppActor::initialiseBackfaceStyling()
{
    switch (visualisation_.backfaceStyling.getValue())
    {
        default:
        case BackfaceStyling::FOLLOW_FRONTFACE:
            actor_->GetProperty()->FrontfaceCullingOff();
            actor_->GetProperty()->BackfaceCullingOff();
            break;
        case BackfaceStyling::CULL_BACKFACE:
            actor_->GetProperty()->FrontfaceCullingOff();
            actor_->GetProperty()->BackfaceCullingOn();
            break;
        case BackfaceStyling::CULL_FRONTFACE:
            actor_->GetProperty()->FrontfaceCullingOn();
            actor_->GetProperty()->BackfaceCullingOff();
            break;
    }
}

void rtppActor::initialiseActor()
{
    actor_ = vtkSmartPointer<vtkActor>::New();
    actor_->GetProperty()->SetInterpolationToPhong();
    actor_->GetProperty()->SetOpacity(visualisation_.opacity);
    actor_->GetProperty()->SetRepresentationToSurface();
    actor_->GetProperty()->EdgeVisibilityOff();
    actor_->GetProperty()->LightingOn();
    actor_->GetProperty()->SetAmbient(DEFAULT_AMBIENT);
    actor_->GetProperty()->SetDiffuse(DEFAULT_DIFFUSE);
    actor_->GetProperty()->SetSpecular(DEFAULT_SPECULAR);
    actor_->GetProperty()->SetSpecularPower(DEFAULT_SPECULAR_POWER);
    actor_->GetProperty()->SetColor(ENGYS_CYAN_DARK[0], ENGYS_CYAN_DARK[1], ENGYS_CYAN_DARK[2]);
    actor_->GetProperty()->SetEdgeColor(BLACK[0], BLACK[1], BLACK[2]);
    actor_->GetProperty()->SetLineWidth(visualisation_.lineThickness);
    actor_->SetVisibility(visualisation_.visible);

    if (visualisation_.opacity == 0)
    {
        actor_->SetVisibility(false);
    }
}

void rtppActor::initialiseMapper()
{
    if (visualisation_.surfaceLic.isVisible())
    {
        mapper_ = createSurfaceLICMapper();
        mapper_->InterpolateScalarsBeforeMappingOff();
    }
    else
    {
        mapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper_->InterpolateScalarsBeforeMappingOn();
    }

    scalar relativeTopologyShift = -static_cast<scalar>(coincidentResolutionFactor) / static_cast<scalar>(64.0);
    coincidentResolutionFactor++;
    if (coincidentResolutionFactor >= 32)
    { coincidentResolutionFactor = 0; }

    scalar polyFactor = static_cast<scalar>(visualisation_.representation.isSurfaceWithEdges()) * static_cast<scalar>(0.5);
    mapper_->SetRelativeCoincidentTopologyPolygonOffsetParameters(polyFactor, relativeTopologyShift + polyFactor + 1);
    mapper_->SetRelativeCoincidentTopologyLineOffsetParameters(0, relativeTopologyShift);
    mapper_->SetRelativeCoincidentTopologyPointOffsetParameter(relativeTopologyShift);

    if (visualisation_.colourField.isSolidColour())
    {
        mapper_->ScalarVisibilityOff();
        colourBySolidColour();
    }
    else
    {
        mapper_->UseLookupTableScalarRangeOn();
        mapper_->ScalarVisibilityOn();
        mapper_->SetColorModeToMapScalars();

        // Use either point or cell data
        // - if both point and cell data exists, preferentially choose
        //   point data.  This is often the case when using glyphs.
        if (visualisation_.colourField.isPointAssociation())
        {
            mapper_->SetScalarModeToUsePointFieldData();
        }
        else if (visualisation_.colourField.isCellAssociation())
        {
            mapper_->SetScalarModeToUseCellFieldData();
        }
        else
        {
            WarningInFunction << "Unknown field association \""
                              << visualisation_.colourField.getAssociation()
                              << "\" when adding colour field \""
                              << visualisation_.colourField << "\" to mapper,"
                              << " assuming point data" << endl;
        }

        mapper_->SelectColorArray(visualisation_.colourField.getFoamName().c_str());
    }
    actor_->SetMapper(mapper_);
}

vtkSmartPointer<vtkPolyDataMapper> rtppActor::createSurfaceLICMapper() const
{
    vtkNew<vtkCompositeSurfaceLICMapper> surfaceLicMapper;
    surfaceLicMapper->SetInputArrayToProcess(
        0,
        0,
        0,
        vtkDataObject::FIELD_ASSOCIATION_POINTS,
        visualisation_.surfaceLic.vectorField.c_str()
    );
    surfaceLicMapper->GetLICInterface()->SetStepSize(visualisation_.surfaceLic.stepSize);
    surfaceLicMapper->GetLICInterface()->SetNumberOfSteps(visualisation_.surfaceLic.numberOfSteps);
    return surfaceLicMapper;
}

void rtppActor::initialiseActiveMapper()
{
    if (visualisation_.showActiveGIBBoundary)
    {
        activeGIBBoundaryMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
        activeGIBBoundaryMapper_->ScalarVisibilityOff();
        double units, factor;
        mapper_->GetRelativeCoincidentTopologyPolygonOffsetParameters(factor, units);
        scalar polyFactor =
            static_cast<scalar>(visualisation_.representation.isSurface() || visualisation_.representation.isSurfaceWithEdges()) * static_cast<scalar>(0.5);
        activeGIBBoundaryMapper_->SetRelativeCoincidentTopologyPolygonOffsetParameters(polyFactor, units);
        mapper_->GetRelativeCoincidentTopologyLineOffsetParameters(factor, units);
        activeGIBBoundaryMapper_->SetRelativeCoincidentTopologyLineOffsetParameters(factor, units);
    }
}

void rtppActor::initialiseRepresentation()
{
    int representation = visualisation_.representation.getValue();
    if (visualisation_.surfaceLic.isVisible())
    {
        representation = Representation::Value::SURFACE;
    }
    switch (representation)
    {
        case Representation::Value::NONE:
            actor_->VisibilityOff();
            break;
        case Representation::Value::WIREFRAME:
            // Note: colour is set using general SetColor, not SetEdgeColor
            actor_->GetProperty()->SetRepresentationToWireframe();
            break;
        default:
        case Representation::Value::SURFACE:
            actor_->GetProperty()->SetRepresentationToSurface();
            break;
        case Representation::Value::SURFACE_WITH_EDGES:
            actor_->GetProperty()->SetRepresentationToSurface();
            actor_->GetProperty()->EdgeVisibilityOn();
            break;
        case Representation::Value::PROFILE:
        {
            //actor_->GetProperty()->SetRepresentationToSurface();
            actor_->GetProperty()->EdgeVisibilityOff();

            edgeFilter_ = vtkSmartPointer<vtkFeatureEdges>::New();
            edgeFilter_->BoundaryEdgesOn();
            edgeFilter_->FeatureEdgesOn();
            edgeFilter_->ManifoldEdgesOff();
            edgeFilter_->NonManifoldEdgesOff();
            edgeFilter_->SetFeatureAngle(30);
            vtkNew<vtkIncrementalOctreePointLocator> locator;
            edgeFilter_->SetLocator(locator);
            edgeFilter_->ColoringOff();
            edgeFilter_->RemoveGhostInterfacesOff();
            break;
        }
        case Representation::Value::OUTLINE:
            actor_->GetProperty()->SetRepresentationToSurface();
            actor_->GetProperty()->EdgeVisibilityOff();
            break;
    }
}

void rtppActor::initialiseNormalsActor()
{
    if (visualisation_.showNormals)
    {
        this->normalsFilter_ = vtkSmartPointer<vtkPolyDataNormals>::New();
        normalsFilter_->ComputePointNormalsOff();
        normalsFilter_->ComputeCellNormalsOn();
        normalsFilter_->SplittingOff();
        normalsFilter_->FlipNormalsOff();
        normalsFilter_->ConsistencyOn();
        normalsFilter_->AutoOrientNormalsOff();

        vtkNew<vtkCellCenters> normalsCenters;
        normalsCenters->SetInputConnection(normalsFilter_->GetOutputPort());

        vtkNew<vtkMaskPoints> normalsMask;
        normalsMask->SetInputConnection(normalsCenters->GetOutputPort());
        normalsMask->GenerateVerticesOn();
        normalsMask->RandomModeOff();

        int onRatio;
        if (visualisation_.normalsRatio <= 0)
        {
            onRatio = -1;
        }
        else if (visualisation_.normalsRatio >= 1)
        {
            onRatio = 1;
        }
        else
        {
            onRatio = static_cast<int>(std::round(std::pow(2, (1 - visualisation_.normalsRatio) * 10)));
        }

        normalsMask->SetOnRatio(onRatio);

        vtkNew<vtkArrowSource> normalsArrow;
        normalsArrow->SetTipResolution(6);
        normalsArrow->SetShaftResolution(6);
        normalsArrow->SetTipRadius(0.1);
        normalsArrow->SetShaftRadius(0.03);
        normalsArrow->SetTipLength(0.35);
        normalsArrow->Update();

        vtkNew<vtkGlyph3D> normalsSource;
        normalsSource->SetInputConnection(normalsMask->GetOutputPort());
        normalsSource->SetSourceData(normalsArrow->GetOutput());
        normalsSource->SetVectorModeToUseNormal();
        normalsSource->SetScaleModeToDataScalingOff();
        normalsSource->SetScaleFactor(visualisation_.normalsLength);
        normalsSource->OrientOn();

        vtkNew<vtkPolyDataMapper> normalsMapper;
        normalsMapper->SetInputConnection(normalsSource->GetOutputPort());
        normalsMapper->ScalarVisibilityOff();

        double factor, units;
        mapper_->GetRelativeCoincidentTopologyPolygonOffsetParameters(factor, units);
        normalsMapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(factor, units + 1);
        mapper_->GetRelativeCoincidentTopologyLineOffsetParameters(factor, units);
        normalsMapper->SetRelativeCoincidentTopologyLineOffsetParameters(factor, units + 1);
        mapper_->GetRelativeCoincidentTopologyPointOffsetParameter(units);
        normalsMapper->SetRelativeCoincidentTopologyPointOffsetParameter(units + 1);

        normalsActor_ = vtkSmartPointer<vtkActor>::New();
        normalsActor_->SetMapper(normalsMapper);

        vtkNew<vtkProperty> normalsProperty;
        point RGB = visualisation_.normalsColor;
        normalsProperty->SetColor(RGB[0], RGB[1], RGB[2]);
        normalsProperty->SetOpacity(visualisation_.normalsOpacity);
        normalsActor_->SetProperty(normalsProperty);
    }
}

void rtppActor::colourBySolidColour()
{
    point RGB;
    if (visualisation_.colourField.isIndexedColour())
    {
        RGB = getIndexedColor(this->index);
    }
    else
    {
        RGB = visualisation_.colour;
    }
    actor_->GetProperty()->SetColor(RGB[0], RGB[1], RGB[2]);
    RGB = visualisation_.lineColour;
    actor_->GetProperty()->SetEdgeColor(RGB[0], RGB[1], RGB[2]);
}

point rtppActor::getIndexedColor(int index)
{
    static const std::vector<point> indexedColors = initIndexedColorSequence();
    point color = indexedColors[index % indexedColors.size()];
    return color;
}

std::vector<point> rtppActor::initIndexedColorSequence()
{
    float degrees[] = {360, 180, 90, 270, 45, 135};
    std::vector<point> colors(6);
    int i = 0;
    for (point& color : colors) {
        color = HSBtoRGB((degrees[i] / DEGREES), 0.66f, 1.0f) / 256.0F;
        i++;
    }
    return colors;
}

point rtppActor::HSBtoRGB(float hue, float saturation, float brightness)
{
    if (saturation == 0.0F) {
        return {brightness * 255.0F + 0.5F, brightness * 255.0F + 0.5F, brightness * 255.0F + 0.5F};
    } else {
        float h = (hue - (float)std::floor((double)hue)) * 6.0F;
        float f = h - (float)std::floor((double)h);
        float p = brightness * (1.0F - saturation);
        float q = brightness * (1.0F - saturation * f);
        float t = brightness * (1.0F - saturation * (1.0F - f));
        switch ((int)h) {
            case 0:
                return {brightness * 255.0F + 0.5F, t * 255.0F + 0.5F, p * 255.0F + 0.5F};
            case 1:
                return {q * 255.0F + 0.5F, brightness * 255.0F + 0.5F, p * 255.0F + 0.5F};
            case 2:
                return {p * 255.0F + 0.5F, brightness * 255.0F + 0.5F, t * 255.0F + 0.5F};
            case 3:
                return {p * 255.0F + 0.5F, q * 255.0F + 0.5F, brightness * 255.0F + 0.5F};
            case 4:
                return {t * 255.0F + 0.5F, p * 255.0F + 0.5F, brightness * 255.0F + 0.5F};
            case 5:
                return {brightness * 255.0F + 0.5F, p * 255.0F + 0.5F, q * 255.0F + 0.5F};
        }
    }
    return {0, 0, 0};
}

void rtppActor::colourByActiveColour()
{
    point RGB = visualisation_.activeColour;
    actor_->GetProperty()->SetColor(RGB[0], RGB[1], RGB[2]);
    actor_->GetProperty()->SetEdgeColor(RGB[0], RGB[1], RGB[2]);
}


// * * * * * * * * * * * * * * * * Public methods  * * * * * * * * * * * * * * //

void rtppActor::updateLookupTable(vtkScalarsToColors *lookupTable)
{
    if (lookupTable && !visualisation_.colourField.isSolidColour())
    {
        mapper_->SetLookupTable(lookupTable);
    }
}

void rtppActor::updateActiveMapperAndRepresentation(vtkPolyData *pData)
{
    vtkBitArray *activeField = vtkBitArray::SafeDownCast(pData->GetFieldData()->GetArray(activePatchKeys::ACTIVE_FIELD_KEY));
    if (activeGIBBoundaryMapper_ && activeField && activeField->GetValue(0))
    {
        activeGIBBoundaryMapper_->SetInputData(pData);
        actor_->SetMapper(activeGIBBoundaryMapper_);
        if (visualisation_.representation.isSurface())
        {
            actor_->GetProperty()->EdgeVisibilityOn();
        }
        colourByActiveColour();
    }
    else if (mapper_)
    {
        mapper_->SetInputData(pData);
        actor_->SetMapper(mapper_);
        if (visualisation_.representation.isSurface())
        {
            actor_->GetProperty()->EdgeVisibilityOff();
        }
        colourBySolidColour();
    }
}

std::vector<vtkActor*> rtppActor::getVTKActors()
{
    std::vector<vtkActor*> vtkActors;
    if (actor_)
    {
        vtkActors.push_back(this->actor_);
        if (normalsActor_)
        {
            vtkActors.push_back(this->normalsActor_);
        }
    }
    return vtkActors;
}

std::vector<std::string> rtppActor::getVTKActorNameSuffixes()
{
    std::vector<std::string> vtkActorSuffixes;
    if (actor_)
    {
        vtkActorSuffixes.emplace_back("");
        if (normalsActor_)
        {
            vtkActorSuffixes.emplace_back("_normals");
        }
    }
    return vtkActorSuffixes;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeVis::rtppActor::~rtppActor() = default;

}
// ************************************************************************* //
