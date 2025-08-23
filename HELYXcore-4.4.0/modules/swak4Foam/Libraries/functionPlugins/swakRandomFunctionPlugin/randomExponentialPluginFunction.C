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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2017 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "randomExponentialPluginFunction.H"
#include "FieldValueExpressionDriver.H"
#include "plugins/FieldValuePluginFunction.H"
#include "primitives/random/Random/Random.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam {

typedef randomExponentialPluginFunction<FieldValuePluginFunction,FieldValueExpressionDriver> randomExponentialPluginFunctionField;
defineTemplateTypeNameAndDebug(randomExponentialPluginFunctionField,0);

addNamedToRunTimeSelectionTable(FieldValuePluginFunction, randomExponentialPluginFunctionField , name, randomExponential);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <typename FType,typename DType>
randomExponentialPluginFunction<FType,DType>::randomExponentialPluginFunction(
    const DType &parentDriver,
    const word &name
):
    FType(
        parentDriver,
        name,
        word("volScalarField"),
        string("seed primitive label,halfLife primitive scalar")
    ),
    halfLife_(1),
    seed_(666)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <typename FType,typename DType>
void randomExponentialPluginFunction<FType,DType>::doEvaluationInternal(
    scalarField &f
) {
    label seed=seed_;

    if (seed<=0) {
        seed=this->mesh().time().timeIndex()-seed;
    }

    Random rnd(seed);

    forAll(f,i) {
#ifdef FOAM_RANDOM_CLASS_NEW_INTERFACE
        f[i]=-log(1-rnd.sample01<scalar>())*halfLife_;
#else
        f[i]=-log(1-rnd.scalar01())*halfLife_;
#endif
    }
}

template <typename FType,typename DType>
void randomExponentialPluginFunction<FType,DType>::doEvaluation()
{
    autoPtr<volScalarField> pRandom(
        new volScalarField(
            IOobject(
                "exponentialRandom",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        )
    );

    doEvaluationInternal(
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        const_cast<scalarField&>(pRandom->internalField().field())
#else
        pRandom->internalField()
#endif
    );

    this->result().setObjectResult(pRandom);
}

template <typename FType,typename DType>
void randomExponentialPluginFunction<FType,DType>::setArgument(
    label index,
    const label &val
)
{
    assert(index==0);
    seed_=val;
}

template <typename FType,typename DType>
void randomExponentialPluginFunction<FType,DType>::setArgument(
    label index,
    const scalar &val
)
{
    assert(index==1);
    halfLife_=val;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
