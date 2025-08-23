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
    (c) 1991-2005 OpenCFD Ltd. ICE Stroemungsfoschungs GmbH

Description
    Load a field into memory if it exists

Contributors/Copyright:
    2010-2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

template<class FieldType>
bool loadFieldFunction(
    fvMesh &mesh,
    const word &fName,
    SLPtrList<FieldType> &fieldList
) {
    const Time &runTime=mesh.time();

    IOobject f
        (
            fName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
#ifdef FOAM_HAS_TYPE_HEADER_OK
    f.typeHeaderOk<FieldType>(false);
#else
    f.headerOk();
#endif

    word className=f.headerClassName();

    if (className==FieldType::typeName) {
        Info<< "Reading field " << fName << " of type "
            << FieldType::typeName << endl;
        fieldList.append(
            new FieldType
            (
                IOobject
                (
                    fName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        return true;
    } else {
        return false;
    }
}

template<class FieldType>
bool loadPointFieldFunction(
    fvMesh &mesh,
    pointMesh &pMesh,
    const word &fName,
    SLPtrList<FieldType> &fieldList
) {
    const Time &runTime=mesh.time();

    IOobject f
        (
            fName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
#ifdef FOAM_HAS_TYPE_HEADER_OK
    f.typeHeaderOk<FieldType>(false);
#else
    f.headerOk();
#endif

    word className=f.headerClassName();

    if (className==FieldType::typeName) {
        Info<< "Reading field " << fName << " of type "
            << FieldType::typeName << endl;
        fieldList.append(
            new FieldType
            (
                IOobject
                (
                    fName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pMesh
            )
        );
        return true;
    } else {
        return false;
    }
}

template<class FieldType>
string reloadFieldsFunction(
    fvMesh &mesh,
    SLPtrList<FieldType> &fieldList
) {
    if (fieldList.size()<=0) {
        return "";
    }
    const Time &runTime=mesh.time();
    string result="Reloading type "+FieldType::typeName+": ";

    forAllIter(typename SLPtrList<FieldType>,fieldList,i){
        FieldType &field=(*i);
        const word &fName=field.name();
        result+=fName+" ";

        IOobject f
            (
                fName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );
        if (
#ifdef FOAM_HAS_TYPE_HEADER_OK
            f.typeHeaderOk<FieldType>(false)
#else
            f.headerOk()
#endif
        ) {
            word className=f.headerClassName();

            if (className==FieldType::typeName) {
                result+="ok ";
                FieldType newField(
                    IOobject
                    (
                        fName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                );
                field==newField;
            } else {
                result+="wrong:"+className+" ";
            }
        } else {
            result+="missing ";
        }
        result+="...";
    }
    result+="\n";

    return result;
}


// ************************************************************************* //
