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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2010-2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "CellSetValueExpressionDriver.H"
#include "plugins/CellSetValuePluginFunction.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(CellSetValueExpressionDriver, 0);

word CellSetValueExpressionDriver::driverName_="cellSet";

addNamedToRunTimeSelectionTable(CommonValueExpressionDriver, CellSetValueExpressionDriver, dictionary, cellSet);
addNamedToRunTimeSelectionTable(CommonValueExpressionDriver, CellSetValueExpressionDriver, idName, cellSet);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


    CellSetValueExpressionDriver::CellSetValueExpressionDriver(const cellSet &set,const CellSetValueExpressionDriver& orig)
:
        SetSubsetValueExpressionDriver(orig),
        cellSet_(
            new cellSet(
                dynamic_cast<const fvMesh&>(set.db()),
                //            dynamicCast<const fvMesh&>(set.db()), // doesn't work with f++ 4.2
                //            set.name()+"_copy",
                set.name(),
                set
            )
        )
{}

CellSetValueExpressionDriver::CellSetValueExpressionDriver(const cellSet &set)
:
    SetSubsetValueExpressionDriver(set.name(),NEW),
    cellSet_(
        new cellSet
        (
            dynamic_cast<const fvMesh&>(set.db()),
            //            dynamicCast<const fvMesh&>(set.db()), // doesn't work with gcc 4.2
            //            set.name()+"_copy",
            set.name(),
            set
        )
    )
{}

CellSetValueExpressionDriver::CellSetValueExpressionDriver(const word& id,const fvMesh&mesh, const objectRegistry& obr)
 :
    SetSubsetValueExpressionDriver(id,INVALID),
    cellSet_(
        getTopoSet<cellSet>(
            mesh,
            id,
            origin_
        )
    )
{
}

CellSetValueExpressionDriver::CellSetValueExpressionDriver(const dictionary& dict,const fvMesh&mesh)
 :
    SetSubsetValueExpressionDriver(dict,dict.lookup("setName"),NEW),
    cellSet_(
        getTopoSet<cellSet>(
            regionMesh(
                dict,
                mesh,
                searchOnDisc()
            ),
            dict.lookup("setName"),
            origin_
        )
    )
{
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CellSetValueExpressionDriver::~CellSetValueExpressionDriver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<scalar>> CellSetValueExpressionDriver::getScalarField(
        const word &name,bool oldTime
)
{
    return getFieldInternal<volScalarField,cellSet,scalar>(
        name,cellSet_,oldTime
    );
}

tmp<Field<vector>> CellSetValueExpressionDriver::getVectorField(
    const word &name,bool oldTime
)
{
    return getFieldInternal<volVectorField,cellSet,vector>(
        name,cellSet_,oldTime
    );
}

tmp<Field<tensor>> CellSetValueExpressionDriver::getTensorField(
    const word &name,bool oldTime
)
{
    return getFieldInternal<volTensorField,cellSet,tensor>(
        name,cellSet_,oldTime
    );
}

tmp<Field<symmTensor>> CellSetValueExpressionDriver::getSymmTensorField(
    const word &name,bool oldTime
)
{
    return getFieldInternal<volSymmTensorField,cellSet,symmTensor>(
        name,cellSet_,oldTime
    );
}

tmp<Field<sphericalTensor>> CellSetValueExpressionDriver::getSphericalTensorField(
    const word &name,bool oldTime
)
{
    return getFieldInternal<volSphericalTensorField,cellSet,sphericalTensor>(
        name,cellSet_,oldTime
);
}

tmp<vectorField> CellSetValueExpressionDriver::makePositionField() const
{
    return getFromFieldInternal(this->mesh().C(),cellSet_());
}

tmp<scalarField> CellSetValueExpressionDriver::makeCellVolumeField() const
{
    return getFromFieldInternal(this->mesh().V(),cellSet_());
}


// tmp<vectorField> CellSetValueExpressionDriver::makePointField() const
// {
//     notImplemented("CellSetValueExpressionDriver::makePointField");
// }

tmp<scalarField> CellSetValueExpressionDriver::makeFaceAreaMagField() const
{
    FatalErrorIn("CellSetValueExpressionDriver::makeFaceAreaField()")
        << "cellSet knows nothing about faces"
            << endl
            << exit(FatalError);
    return tmp<scalarField>(new scalarField(0));
}

tmp<scalarField> CellSetValueExpressionDriver::makeFaceFlipField() const
{
    FatalErrorIn("CellSetValueExpressionDriver::makeFaceFlipField()")
        << "cellSet knows nothing about faces"
            << endl
            << exit(FatalError);
    return tmp<scalarField>(new scalarField(0));
}

tmp<vectorField> CellSetValueExpressionDriver::makeFaceNormalField() const
{
    FatalErrorIn("CellSetValueExpressionDriver::makeFaceNormalField()")
        << "cellSet knows nothing about faces"
            << endl
            << exit(FatalError);
    return tmp<vectorField>(new vectorField(0));
}

tmp<vectorField> CellSetValueExpressionDriver::makeFaceAreaField() const
{
    FatalErrorIn("CellSetValueExpressionDriver::makeFaceAreaField()")
        << "cellSet knows nothing about faces"
            << endl
            << exit(FatalError);
    return tmp<vectorField>(new vectorField(0));
}

bool CellSetValueExpressionDriver::update()
{
    if (debug) {
        Pout<< "CellSet: update " << cellSet_->name()
            << endl;
    }

    return updateSet(cellSet_,id_,origin_);
}

autoPtr<CommonPluginFunction> CellSetValueExpressionDriver::newPluginFunction(
    const word &name
) {
    return autoPtr<CommonPluginFunction>(
        CellSetValuePluginFunction::New(
            *this,
            name
        ).ptr()
    );
}

bool CellSetValueExpressionDriver::existsPluginFunction(
    const word &name
) {
    return CellSetValuePluginFunction::exists(
        *this,
        name
    );
}

tmp<scalarField> CellSetValueExpressionDriver::weightsNonPoint(
    label size
) const
{
    const label cellSize=this->size();
    bool isCell=(size==cellSize);
    reduce(isCell,andOp<bool>());

    if (!isCell) {
        Pout<< "Expected size: " << size
            << " Cell size: " << cellSize << endl;

        FatalErrorIn("CellSetValueExpressionDriver::weightsNonPoint")
            << "Can not construct weight field of the expected size. "
                << " For sizes on the processors see above"
                << endl
                << exit(FatalError);
    }

    return tmp<scalarField>(makeCellVolumeField());
}

// ************************************************************************* //

} // namespace
