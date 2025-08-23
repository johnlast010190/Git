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
    (c) 2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "coalCloudList/coalCloudList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalCloudList::coalCloudList
(
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& slgThermo
)
:
    PtrList<coalCloud>(),
    mesh_(rho.mesh())
{
    IOdictionary props
    (
        IOobject
        (
            "coalCloudList",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ
        )
    );

    const wordHashSet cloudNames(wordList(props.lookup("clouds")));

    setSize(cloudNames.size());

    label i = 0;
    forAllConstIter(wordHashSet, cloudNames, iter)
    {
        const word& name = iter.key();

        Info<< "creating cloud: " << name << endl;

        set
        (
            i++,
            new coalCloud
            (
                name,
                rho,
                U,
                g,
                slgThermo
            )
        );

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalCloudList::evolve()
{
    forAll(*this, i)
    {
        operator[](i).evolve();
    }
}


// ************************************************************************* //
