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
    (c) 2022 OpenFOAM Foundation
    (c) 2014-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "checkMeshQuality.H"
#include "checkTools.H"
#include "sets/topoSets/faceSet.H"
#include "motionSmoother/motionSmoother.H"
#include "sampledSurface/writers/surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::checkMeshQuality
(
    const fvMesh& mesh,
    const dictionary& dict,
    const wordList& metrics,
    const bool& report,
    const bool& writeCombinedMetric,
    const autoPtr<surfaceWriter>& surfWriter
)
{
    label noFailedChecks = 0;

    {
        PtrList<meshStatistics> metricsFields(0);

        if (metrics.size())
        {
            metricsFields.setSize(metrics.size());
            forAll(metrics, metricI)
            {
                word metricName = metrics[metricI];
                dictionary histSettings;
                if (dict.found("histograms"))
                {
                    dictionary allHist = dict.subDict("histograms");
                    if (allHist.found(metricName))
                    {
                        histSettings = allHist.subDict(metricName);
                    }
                }

                metricsFields.set
                (
                    metricI,
                    new meshStatistics
                    (
                        mesh,
                        metricName,
                        histSettings
                    )
                );
            }
            if (writeCombinedMetric)
            {
                metricsFields.setSize(metrics.size()+1);
                word metricName("meshQuality");
                dictionary histSettings;
                if (dict.found("histograms"))
                {
                    dictionary allHist = dict.subDict("histograms");
                    if (allHist.found(metricName))
                    {
                        histSettings = allHist.subDict(metricName);
                    }
                }

                metricsFields.set
                (
                    metrics.size(),
                    new meshStatistics
                    (
                        mesh,
                        metricName,
                        histSettings
                    )
                );
            }
        }

        faceSet faces(mesh, "meshQualityFaces", mesh.nFaces()/100+1);
        motionSmoother::checkMesh
        (
            false,
            mesh,
            dict,
            faces,
            report,
            &metricsFields
        );

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            noFailedChecks++;

            Info<< "  <<Writing " << nFaces
                << " faces in error to set " << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter.valid())
            {
                mergeAndWrite(surfWriter(), faces);
            }
        }

        if (metrics.size())
        {
            Info<<"\nOutputting mesh quality fields... "<<nl<<endl;

            forAll(metrics, metricI)
            {
                Info<<"Writing field: "<<metrics[metricI]<<endl;
                metricsFields[metricI].field().write();
                metricsFields[metricI].writeStats();
            }

            if (writeCombinedMetric)
            {
                volScalarField& mf = metricsFields[metrics.size()].field();
                Info<<"Writing field: "<<mf.name()<<endl;
                mf /= metrics.size();
                mf.write();
                metricsFields[metrics.size()].calcStats(mf);
                metricsFields[metrics.size()].writeStats();
            }
        }
    }

    return noFailedChecks;
}


// ************************************************************************* //
