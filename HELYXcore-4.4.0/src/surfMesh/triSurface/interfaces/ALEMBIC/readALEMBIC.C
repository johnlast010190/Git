/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "alembic/AlembicReader.H"
#include "triSurface/triSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readALEMBIC
(
    const fileName& ALEMBICfileName,
    const dictionary& dict
)
{
    // create alembic reader
    fileFormats::AlembicReader reader(ALEMBICfileName);

    // read alembic setup from file (target time)
    dictionary abcsetup;
    if (dict.toc().size())
    {
        abcsetup = dict;
    }
    else
    {
        fileName alembicSetupFile(ALEMBICfileName.path()/"alembicSetup");
        IFstream setupFile(alembicSetupFile);
        if (!setupFile.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << alembicSetupFile
                << exit(FatalError);
        }
        //abcsetup = dictionary(setupFile).subDict(ALEMBICfileName.nameLessExt());
        abcsetup = dictionary(setupFile).subDict(ALEMBICfileName.name());
    }

    if (abcsetup.found("time"))
    {
        scalar t(abcsetup.lookup<scalar>("time"));
        Info<< "Access .abc data for time " << t << endl;
        reader.setTimeData((Abc::chrono_t) t);
    }
    else if (abcsetup.found("timeIndex"))
    {
        label timeIdx(abcsetup.lookup<label>("timeIndex"));
        Info<< "Access .abc data for time index " << timeIdx << endl;
        reader.setTimeData((index_t) timeIdx);
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find entry time or timeIndex in file "
            << abcsetup.name() << exit(FatalError);
    }

    // access mesh data
    DynamicList<point>& points = reader.points();
    DynamicList<labelledTri>& faces = reader.faces();
    HashTable<label>& groupToPatch = reader.groupToPatch();

    // Convert groupToPatch to patchList.
    label maxGroupID = groupToPatch.size();
    geometricSurfacePatchList patches(maxGroupID);

    if (maxGroupID == 0)
    {
        // Add single (default) patch
        patches = { geometricSurfacePatch("patch0", 0) };
    }
    else
    {
        forAllConstIters(groupToPatch, iter)
        {
            const label patchIdx = iter.object();
            patches[patchIdx] = geometricSurfacePatch
            (
                iter.key(),
                patchIdx
            );
        }
    }

    // Transfer DynamicLists to straight ones.
    pointField allPoints(std::move(points));
    *this = triSurface(faces, patches, allPoints, true);

    return true;
}


// ************************************************************************* //
