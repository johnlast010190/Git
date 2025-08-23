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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2019 Engys Ltd.

Description
    Abaqus surface writer.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeINP(Ostream& os) const
{
    //Write out header
    os << "**" <<endl;
    os << "** OpenFOAM generated ABAQUS file" <<endl;
    os << "**" <<endl;

    labelList faceMap;

    surfacePatchList myPatches(calcPatches(faceMap));

    label faceIndex = 0;

    const pointField& ps = points();

    os.setf(ios_base::scientific);
    // Write vertex coords
    os  <<setw(5)<<"*NODE"<<endl;
    forAll(ps, pointi)
    {
        //Abaqus does not like index starting
        //from 0, this is why +1
        os  << pointi+1<<", "
            << ps[pointi].x()<<", "
            << ps[pointi].y()<<", "
            << ps[pointi].z()
            << endl;
    }

    // Write elements
    os  <<"**"<<endl;
    os  <<"*************E L E M E N T S *************"<<endl;
    os  <<"**"<<endl;

    forAll(myPatches, patchI)
    {
        // Print all faces belonging to this region
        const surfacePatch& patch = myPatches[patchI];
        os  << setw(8) << "*ELEMENT, "<<"TYPE=S3"
            <<", ELSET="<<patch.name()<<endl;
        for
        (
            label patchFaceI = 0;
            patchFaceI < patch.size();
            patchFaceI++
        )
        {
            const label faceI = faceMap[faceIndex++];
            const labelledTri& f = (*this)[faceI];

            os  << faceIndex<<", "
                << f[0]+1<<", "
                << f[1]+1<<", "
                << f[2]+1
                << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
