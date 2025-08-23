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
    2012-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "readAndUpdateFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class FType>
void Foam::readAndUpdateFields::correctBoundaryConditions(
    PtrList<FType> &flst
) {
    forAll(flst,i)
    {
        Info<< " " << flst[i].name() << ":" << FType::typeName << flush;

    flst[i].correctBoundaryConditions();
        if (
            this->obr_.time().outputTime()
            &&
            flst[i].writeOpt()==IOobject::AUTO_WRITE
        ) {
            Info<< "(writing)" << flush;
            flst[i].write();
        }
    }
}

template<class Type>
bool Foam::readAndUpdateFields::loadField
(
    const word& fieldName,
    PtrList<VolField<Type>>& vflds,
    PtrList<SurfaceField<Type>>& sflds,
    PtrList<PointField<Type>>& pflds
)
{
    typedef VolField<Type> vfType;
    typedef SurfaceField<Type> sfType;
    typedef PointField<Type> pfType;

    if (obr_.foundObject<vfType>(fieldName))
    {
        Info<< this->name() << "  : Field " << fieldName << " of type "
            << vfType::typeName << " already in database" << endl;
        return true;
    }
    else if (obr_.foundObject<sfType>(fieldName))
    {
        Info<< this->name() << "  : Field " << fieldName << " of type "
            << sfType::typeName << " already in database" << endl;
        return true;
    }
    else if (obr_.foundObject<pfType>(fieldName))
    {
        Info<< this->name() << "  : Field " << fieldName << " of type "
            << pfType::typeName << " already in database" << endl;
        return true;
    }
    else
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if
        (
#ifdef FOAM_HAS_TYPE_HEADER_OK
            fieldHeader.typeHeaderOk<vfType>()
#else
            fieldHeader.headerOk()
            && fieldHeader.headerClassName() == vfType::typeName
#endif
        )
        {
            // store field locally
            Info<< this->name() << "  : Field " << fieldName << " of type "
                << vfType::typeName << " being read" << endl;

            label sz = vflds.size();
            vflds.setSize(sz+1);
            vflds.set(sz, new vfType(fieldHeader, mesh));
            //            vflds[sz].writeOpt()=IOobject::AUTO_WRITE;

            return true;
        }
        else if
        (
#ifdef FOAM_HAS_TYPE_HEADER_OK
            fieldHeader.typeHeaderOk<pfType>()
#else
            fieldHeader.headerOk()
            && fieldHeader.headerClassName() == pfType::typeName
#endif
        )
        {
            // store field locally
            Info<< this->name() << "  : Field " << fieldName << " of type "
                << pfType::typeName << " being read" << endl;
            label sz = pflds.size();
            pflds.setSize(sz+1);
            pflds.set(sz, new pfType(fieldHeader, this->pMesh(mesh)));
            //            pflds[sz].writeOpt()=IOobject::AUTO_WRITE;

            return true;
        }
        else if
        (
#ifdef FOAM_HAS_TYPE_HEADER_OK
            fieldHeader.typeHeaderOk<sfType>()
#else
            fieldHeader.headerOk()
            && fieldHeader.headerClassName() == sfType::typeName
#endif
        )
        {
            // store field locally
            Info<< this->name() << "  : Field " << fieldName << " of type "
                << sfType::typeName << " being read" << endl;
            label sz = sflds.size();
            sflds.setSize(sz+1);
            sflds.set(sz, new sfType(fieldHeader, mesh));
            //            pflds[sz].writeOpt()=IOobject::AUTO_WRITE;

            return true;
        }
    }
    return false;
}


// ************************************************************************* //
