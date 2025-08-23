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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fv/fv.H"
#include "db/objectRegistry/objectRegistry.H"
#include "matrices/solution/solution.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fv::gradScheme<Type>> Foam::fv::gradScheme<Type>::New
(
    const fvMesh& mesh,
    const objectRegistry& db,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing gradScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Grad scheme not specified" << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    const auto ctor =
        ctorTableLookup("grad scheme", IstreamConstructorTable_(), schemeName);
    return ctor(mesh, db, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::gradScheme<Type>::~gradScheme()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::grad
(
    const VolField<Type>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef VolField<GradType> GradFieldType;

    if (!this->mesh().changing() && this->mesh().solution().cache(name))
    {
        if (!vsf.db().objectRegistry::template foundObject<GradFieldType>(name))
        {
            solution::cachePrintMessage("Calculating and caching", name, vsf);
            tmp<GradFieldType> tgGrad = calcGrad(vsf, name);
            tgGrad->checkIn();
            regIOobject::store(tgGrad.ptr());
        }

        solution::cachePrintMessage("Retrieving", name, vsf);
        GradFieldType& gGrad =
            vsf().db().objectRegistry::template lookupObjectRef<GradFieldType>
            (
                name
            );

        if (gGrad.upToDate(vsf))
        {
            return gGrad;
        }
        else
        {
            solution::cachePrintMessage("Deleting", name, vsf);
            // Deletes object owned by registry
            gGrad.checkOut();

            solution::cachePrintMessage("Recalculating", name, vsf);
            tmp<GradFieldType> tgGrad = calcGrad(vsf, name);

            solution::cachePrintMessage("Storing", name, vsf);
            tgGrad->checkIn();
            regIOobject::store(tgGrad.ptr());
            GradFieldType& gGrad =
                vsf().db().objectRegistry::template
                lookupObjectRef<GradFieldType>(name);

            return gGrad;
        }
    }
    else
    {
        if
        (
            vsf().db().objectRegistry::template foundObject<GradFieldType>
            (
                name
            )
        )
        {
            GradFieldType& gGrad =
            (
                vsf.db().objectRegistry::template
                lookupObjectRef<GradFieldType>(name)
            );

            if (gGrad.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vsf);
                gGrad.checkOut();
            }
        }

        solution::cachePrintMessage("Calculating", name, vsf);
        return calcGrad(vsf, name);
    }
}


template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::grad(const VolField<Type>& vsf) const
{
    return grad(vsf, "grad(" + vsf.name() + ')');
}


template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::grad(const tmp<VolField<Type>>& tvsf) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef VolField<GradType> GradFieldType;

    tmp<GradFieldType> tgrad = grad(tvsf());
    tvsf.clear();
    return tgrad;
}


template<class Type>
Foam::tmp
<
    Foam::BlockLduSystem
    <
        Foam::vector, typename Foam::outerProduct<Foam::vector, Type>::type
    >
>
Foam::fv::gradScheme<Type>::fvmGrad(const VolField<Type>& vf) const
{
    FatalErrorInFunction
        << "Only Gauss linear is supported"
        << exit(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType>> tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
}


// ************************************************************************* //
