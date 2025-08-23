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
    2010-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "swakExpressionFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swakExpressionFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        swakExpressionFunctionObject,
        dictionary
    );



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

swakExpressionFunctionObject::swakExpressionFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    timelineFunctionObject(name,t,dict),
    expression_(
        dict.lookup("expression"),
        dict
    ),
    maskExpression_(
        dict.lookupOrDefault<string>("mask",""),
        dict
    ),
    accumulations_(
        NumericAccumulationNamedEnum::readAccumulations(
            dict,"accumulations"
        )
    ),
    driver_(
        CommonValueExpressionDriver::New(
            dict,
            refCast<const fvMesh>(obr_)
        )
    )
{
    driver_->createWriterAndRead(name+"_"+type());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word swakExpressionFunctionObject::dirName()
{
    return typeName;
}

wordList swakExpressionFunctionObject::fileNames()
{
    return wordList(1,name());
}

stringList swakExpressionFunctionObject::columnNames()
{
    stringList result(accumulations_.size());

    forAll(accumulations_,i) {
        result[i]=NumericAccumulationNamedEnum::toString(
            accumulations_[i]
        );
    }

    return result;
}

void swakExpressionFunctionObject::writeSimple()
{

    if (verbose()) {
        Info<< "Expression " << name() << " : ";
    }

    driver_->clearVariables();
    driver_->parse(expression_);
    word rType=driver_->CommonValueExpressionDriver::getResultType();

    if (rType==pTraits<scalar>::typeName) {
        writeTheData<scalar>(driver_());
    } else if (rType==pTraits<vector>::typeName) {
        writeTheData<vector>(driver_());
    } else if (rType==pTraits<tensor>::typeName) {
        writeTheData<tensor>(driver_());
    } else if (rType==pTraits<symmTensor>::typeName) {
        writeTheData<symmTensor>(driver_());
    } else if (rType==pTraits<sphericalTensor>::typeName) {
        writeTheData<sphericalTensor>(driver_());
    } else {
        WarningIn("swakExpressionFunctionObject::writeSimple()")
            << "Don't know how to handle type " << rType
                << endl;
    }

    if (verbose()) {
        Info<< endl;
    }

    // make sure that the stored Variables are consistently written
    driver_->tryWrite();
}

} // namespace Foam

// ************************************************************************* //
