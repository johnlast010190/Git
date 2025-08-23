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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::rtppActor

Description
    A wrapper around the vtkActor class.

SourceFiles
    rtppActor.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/enums/NamedEnum.H"
#include "baseClasses/visualisation/visualisation.H"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkActor;
class vtkPolyDataMapper;
class vtkDataSet;
class vtkScalarsToColors;
class vtkPolyData;
class engysOutlineFilter;
class vtkFeatureEdges;
class vtkPolyDataNormals;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                          Class rtppActor Declaration
\*---------------------------------------------------------------------------*/

class rtppActor
{
protected:

    // Protected data

    //- Debug switch
    int debug = 0;

    const Visualisation &visualisation_;
    int index = 0;

    vtkSmartPointer<vtkActor> actor_;
    vtkSmartPointer<vtkActor> normalsActor_;
    vtkSmartPointer<vtkPolyDataMapper> mapper_;
    vtkSmartPointer<vtkPolyDataMapper> activeGIBBoundaryMapper_;
    vtkSmartPointer<vtkPolyData> processedInput_;
    vtkSmartPointer<vtkPolyDataNormals> normalsFilter_;

    vtkSmartPointer<vtkFeatureEdges> edgeFilter_;

public:

    // Constructors

    explicit rtppActor(const Visualisation &visualisation, int index);

    //- Destructor
    virtual ~rtppActor();


    // Member Functions

    [[nodiscard]] bool isTransparent() const
    { return visualisation_.isTransparent(); };

    [[nodiscard]] bool isOutline() const
    { return visualisation_.representation.isOutline(); };

    void processInputData(vtkDataSet *input);

    vtkSmartPointer<vtkPolyData> getProcessedInputData();

    void setProcessedInputData(vtkPolyData *input);

    void updateLookupTable(vtkScalarsToColors *lookupTable);

    std::vector<vtkActor*> getVTKActors();
    /**
     * In case there is more than one actor, to differentiate their names
     * each one should get a suffix
     * @return The suffix to be added to the name of each vtkActor
     */
    std::vector<std::string> getVTKActorNameSuffixes();

    //- No copy assignment
    rtppActor &operator=(const rtppActor &) = delete;

protected:

    void initialiseActor();
    void initialiseBackfaceStyling();
    void initialiseMapper();
    void initialiseRepresentation();

    void initialiseNormalsActor();

    void colourBySolidColour();
    void colourByActiveColour();
    void initialiseActiveMapper();
    void updateActiveMapperAndRepresentation(vtkPolyData *pData);
    [[nodiscard]] vtkSmartPointer<vtkPolyDataMapper> createSurfaceLICMapper() const;

    static constexpr int DEGREES = 360;
    static std::vector<point> initIndexedColorSequence();
    static point HSBtoRGB(float hue, float saturation, float brightness);
    static point getIndexedColor(int index);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
