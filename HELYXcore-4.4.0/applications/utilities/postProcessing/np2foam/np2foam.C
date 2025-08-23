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


Description
    Change Numpy array data into OpenFoam format for visualization in Paraview

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "helyxMap.H"
#include "sample/helyxSample/helyxSample.H"
#include "sample/npGrid/npGrid.H"

int main(int argc, char *argv[])
{
	DynamicList<word> argNames;
	DynamicList<word> argValues;

    argList::addOption
    (
        "fileNames",
        "p,U:0,U:1,blockDef",
        "Specify the field and block definition file names to be converted."
    );

    argList::addOption
    (
        "blockDef",
        "blockDef",
        "Specify the block definition file name."
    );

	argList::addOption
    (
        "scalars",
        "p",
        "Specify the scalar field names, separated by comma."
    );

	argList::addOption
    (
        "vectors",
        "U",
        "Specify the vector field names, separated by comma."
    );


	#include "include/setRootCase.H"

	if (argc>1)
    {
		for (label i=1;i<argc-1;i+=2)
		{
			argNames.append(argv[i]);
			argValues.append(argv[i+1]);
		}
	}
	word flist="p,U:0,U:1,blockDef";
	word blockDef="blockDef";

	if (argNames.size()>=1)
	{
		for (label i=0;i<argNames.size();i++)
		{
			if (argNames[i]=="-fileNames")
			{
				flist=argValues[i];
				continue;
			}

			if (argNames[i]=="-blockDef")
			{
				blockDef=argValues[i];
				continue;
			}

		}
	}

	npGrid grid(blockDef,flist);
	grid.writeVisualization();

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
