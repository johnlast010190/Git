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

#include "porosityZones.H"

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

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<porosityZone>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityZones::porosityZones
(
    const fvMesh& mesh
)
:
    IOPtrList<porosityZone>
    (
        IOobject
        (
            "porosityZones",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        porosityZone::iNew(mesh)
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> Foam::porosityZones::porosity() const
{
    // Create the return field
    tmp<volScalarField> tporosity
    (
        volScalarField::New
        (
            "porosity",
            mesh_,
            dimensionedScalar(dimless, 1),
            "zeroGradient"
        )
    );

    // Make reference to the return field (stripping tmp-nature)
    volScalarField & poro( tporosity() );

    // Loop over all zones
    forAll(*this, i)
    {
        operator[](i).porosity( poro );
    }

    return tporosity;
}


void Foam::porosityZones::addResistance(fvVectorMatrix& UEqn) const
{
    forAll(*this, i)
    {
        operator[](i).addResistance(UEqn);
    }
}

void Foam::porosityZones::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    // addResistance for each zone, delaying the correction of the
    // precessor BCs of AU
    forAll(*this, i)
    {
        operator[](i).addResistance(UEqn, AU, false);
    }

    // Correct the boundary conditions of the tensorial diagonal to ensure
    // processor bounaries are correctly handled when AU^-1 is interpolated
    // for the pressure equation.
    AU.correctBoundaryConditions();
}


bool Foam::porosityZones::readData(Istream& is)
{
    clear();

    IOPtrList<porosityZone> newLst
    (
        IOobject
        (
            "porosityZones",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false     // Don't re-register new zones with objectRegistry
        ),
        porosityZone::iNew(mesh_)
    );

    transfer(newLst);

    return is.good();
}


bool Foam::porosityZones::writeData(Ostream& os, bool subDict) const
{
    // Write size of list
    os << nl << size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, subDict);
    }

    // Write end of contents
    os << token::END_LIST << nl;

    // Check state of IOstream
    return os.good();
}


// ************************************************************************* //
