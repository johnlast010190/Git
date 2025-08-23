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
    (c) held by original author
    (c) 2022 Engys Ltd.

Application
    faceSetToSTL

Description
    A small utility to create simple STL-surfaces from a pointField and a
    faceList

    It can also based on a faceList with one entry perform a simple trans-
    lation of this face and create a closed STL. See the example file
    "stlDefinitions".

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "foamTime.H"
    #else
        #include "db/Time/Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "db/Time/Time.H"
#else
    #include "db/Time/Time.H"
#endif

//#if EXTBRANCH==1 && OFVERSION>310
//    #include "foamTime.H"
//#else
//    #include "db/Time/Time.H"
//#endif

#include "fvMesh/fvMesh.H"
#include "algorithms/polygonTriangulate/polygonTriangulate.H"
#include "triSurface/triSurface.H"

using namespace Foam;

void extrudeFacesAndPoints
(
    const dictionary&,
    faceList&,
    pointField&
);

int main(int argc, char *argv[])
{

#include "include/addTimeOptions.H"
#include "include/setRootCase.H"

#include "include/createTime.H"
#include "include/createMesh.H"

    // Create a triangulation engine
    polygonTriangulate triEngine;

    IOdictionary stlDefs
    (
        IOobject
        (
            "stlDefinitions",
            runTime.constant(),
            "triSurface",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList toc = stlDefs.toc();

    forAll(toc, item)
    {
        if (stlDefs.isDict(toc[item]))
        {
            Info<< "\nCreates the STL surface for " << toc[item] << endl;

            pointField pp(stlDefs.subDict(toc[item]).lookup("points"));
            faceList faces(stlDefs.subDict(toc[item]).lookup("faces"));

            triFaceList tfl(0);
            label count(0);

            if (
                   stlDefs.subDict(toc[item])
                   .lookupOrDefault<Switch>("extrude", false)
               )
            {
                if (faces.size() <= 1)
                {
                    extrudeFacesAndPoints
                        (
                            stlDefs.subDict(toc[item]),
                            faces,
                            pp
                        );
                }
                else
                {
                    Info<< "\nWARNING: Using extrude, but"
                         << " multiple faces are defined\n" << endl;
                }
            }

            forAll(faces, facei)
            {
                const face& f = faces[facei];

                triEngine.triangulate
                (
                    UIndirectList<point>(pp, f)
                );

                tfl.setSize(count + triEngine.triPoints().size());

                forAll(triEngine.triPoints(), triI)
                {
                    tfl[count++] = triEngine.triPoints(triI, f);
                }
            }

            triSurface ts(tfl, pp);

            Info<< "Writes the STL surface for " << toc[item] << endl;

            ts.write( "constant/triSurface/"+toc[item]+".stl" );
        }
    }

    Info<< nl << "End" << endl;

    return 0;
}

void extrudeFacesAndPoints
(
    const dictionary& dict, faceList& fL, pointField& pp
)
{
	// Get the extrude vector
	vector extrude( dict.lookup("extrudeVector") );

	// Check for correct (positive) orientation of the final shape
    vector faceNormal(fL[0].areaNormal(pp));

    scalar projection = (extrude & faceNormal)
    	/(Foam::mag(extrude)*Foam::mag(faceNormal));

    // If the projection is positive (!), the face is ordered correctly.
    // Note that this is the case, since the original face is swapped in
    // orientation at the end of this method.
    if (projection < -SMALL)
    {
        Info<< "Orientation of the STL-surface is swapped to\n"
        	 << "yield outward pointing normals." << endl;

        // Swap the direction of the original face
        face dummy = fL[0];
        face& target = fL[0];

        forAll(target, pointi)
        {
          	target[target.size() - 1 - pointi] = dummy[pointi];
        }
    }


	// Get the number of points and extrude the points
    label N = pp.size();
    pp.setSize(2*N);

    for (int i=0; i < N; i++)
    {
        pp[N + i] = pp[i] + extrude;
    }

    label M = fL[0].size();

    fL.setSize(2 + M);

    fL[1].setSize(M);

    face& fOrg(fL[0]);
    face& fExt(fL[1]);

    forAll(fOrg, pointi)
    {
//        fExt[N - 1 - pointi] = fOrg[pointi] + N;
        fExt[pointi] = fOrg[pointi] + N;
    }

    for (int i = 0; i < M ; i++)
    {
        face& f(fL[i + 2]);
        f.setSize(4);

        f[0] = fOrg[i];
        f[1] = fOrg[(i + 1)%M];
        f[2] = fExt[(i + 1)%M];
        f[3] = fExt[i];
    }

    // Swap the direction of the original face
    face dummy = fOrg;

    forAll(fOrg, pointi)
    {
//        fExt[N - 1 - pointi] = fOrg[pointi] + N;
    	fOrg[N - 1 - pointi] = dummy[pointi];
//        fExt[pointi] = fOrg[pointi] + N;
    }
}
