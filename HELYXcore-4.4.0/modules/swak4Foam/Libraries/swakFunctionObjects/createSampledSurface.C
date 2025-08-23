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
    2010, 2013-2014, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "createSampledSurface.H"

#include "repositories/SurfacesRepository.H"

namespace Foam {
    defineTypeNameAndDebug(createSampledSurface,0);
}

Foam::createSampledSurface::createSampledSurface
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
):
    active_(true),
    obr_(obr),
    surfaceName_("noSurfaceHereIn_"+name)
{
    if (!isA<fvMesh>(obr))
    {
        active_=false;
        WarningIn("createSampledSurface::createSampledSurface")
            << "Not a fvMesh. Nothing I can do"
                << endl;
    }
    read(dict);
    execute();
}

Foam::createSampledSurface::~createSampledSurface()
{}

void Foam::createSampledSurface::timeSet()
{
    // Do nothing
}

void Foam::createSampledSurface::read(const dictionary& dict)
{
    if (active_) {
        surfaceName_=word(dict.lookup("surfaceName"));

        SurfacesRepository::getRepository(obr_).getSurface(
            dict,
            dynamic_cast<const fvMesh &>(obr_)
        );
    }
}

void Foam::createSampledSurface::execute()
{
}


void Foam::createSampledSurface::end()
{
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::createSampledSurface::write()
{
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}

void Foam::createSampledSurface::clearData()
{
}

//- Update for changes of mesh
void Foam::createSampledSurface::topoChange(const polyTopoChangeMap&)
{
    redoSurface();
}

//- Update for changes of mesh
void Foam::createSampledSurface::mapMesh(const polyMeshMap&)
{
    redoSurface();
}

//- Update for changes of mesh
#ifdef FOAM_MOVEPOINTS_GETS_POLYMESH
void Foam::createSampledSurface::movePoints(const polyMesh&)
#else
void Foam::createSampledSurface::movePoints(const pointField&)
#endif
{
    redoSurface();
}

void Foam::createSampledSurface::redoSurface()
{
    Info<< "Forcing regeneration of surface " << surfaceName_ << endl;

    bool status=SurfacesRepository::getRepository(obr_).updateSurface(
        surfaceName_,
        dynamic_cast<const fvMesh &>(obr_)
    );
    Info<< "Forcing status for " << surfaceName_ << ": " << status << endl;
}

// ************************************************************************* //
