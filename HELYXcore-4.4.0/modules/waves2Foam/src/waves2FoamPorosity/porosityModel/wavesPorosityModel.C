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

\*---------------------------------------------------------------------------*/

#include "wavesPorosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    defineTypeNameAndDebug(wavesPorosityModel, 0);
    defineRunTimeSelectionTable(wavesPorosityModel, wavesPorosityModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::wavesPorosityModel::wavesPorosityModel
(
	const fvMesh& mesh
)
:
    porosity_
    (
        IOobject
        (
            "porosity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1),
        "zeroGradient"
    )
{
}


Foam::wavesPorosityModel::~wavesPorosityModel()
{}


autoPtr<wavesPorosityModel> wavesPorosityModel::New
(
    const fvMesh& mesh
)
{
    word wavesPorosityModelTypeName;

    // Enclose the creation of the dictionary to ensure it is deleted before
    // the actual porosity model is created
    {
        if (mesh.thisDb().foundObject<IOdictionary>("waveProperties"))
        {
        	wavesPorosityModelTypeName =
                mesh.thisDb().lookupObject<IOdictionary>
                (
                    "waveProperties"
                ).lookup<word>("porosityModel");
        }
        else
        {
        	IOdictionary wp
        	(
        	    IOobject
        	    (
        	        "waveProperties",
        	        mesh.time().constant(),
        	        mesh,
        	        IOobject::MUST_READ,
        	        IOobject::NO_WRITE
        	    )
        	);

        	wavesPorosityModelTypeName = wp.lookup<word>("porosityModel");
        }
    }

    const auto ctor =
        ctorTableLookup
        (
            "porosity model type",
            wavesPorosityModelConstructorTable_(),
            wavesPorosityModelTypeName
        );
    return autoPtr<wavesPorosityModel>(ctor(mesh));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
