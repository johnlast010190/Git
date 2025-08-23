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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "include/swak.H"
#include "primitives/contiguous/contiguous.H"
#include "findShiftedFieldPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshSearch/meshSearch.H"

namespace Foam {

    typedef findShiftedFieldPluginFunction<scalar> findShiftedFieldScalar;
    defineTemplateTypeNameAndDebug(findShiftedFieldScalar,0);
    addNamedToRunTimeSelectionTable(FieldValuePluginFunction, findShiftedFieldScalar , name, findShiftedScalarField);

    typedef findShiftedFieldPluginFunction<vector> findShiftedFieldVector;
    defineTemplateTypeNameAndDebug(findShiftedFieldVector,0);
    addNamedToRunTimeSelectionTable(FieldValuePluginFunction, findShiftedFieldVector , name, findShiftedVectorField);

    typedef findShiftedFieldPluginFunction<tensor> findShiftedFieldTensor;
    defineTemplateTypeNameAndDebug(findShiftedFieldTensor,0);
    addNamedToRunTimeSelectionTable(FieldValuePluginFunction, findShiftedFieldTensor , name, findShiftedTensorField);

    typedef findShiftedFieldPluginFunction<symmTensor> findShiftedFieldSymmTensor;
    defineTemplateTypeNameAndDebug(findShiftedFieldSymmTensor,0);
    addNamedToRunTimeSelectionTable(FieldValuePluginFunction, findShiftedFieldSymmTensor , name, findShiftedSymmTensorField);

    typedef findShiftedFieldPluginFunction<sphericalTensor> findShiftedFieldSphericalTensor;
    defineTemplateTypeNameAndDebug(findShiftedFieldSphericalTensor,0);
    addNamedToRunTimeSelectionTable(FieldValuePluginFunction, findShiftedFieldSphericalTensor , name, findShiftedSphericalTensorField);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
findShiftedFieldPluginFunction<Type>::findShiftedFieldPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        word(ResultType::typeName),
        "field internalField "+ResultType::typeName
        +",shift internalField volVectorField"
        +",default internalField "+ResultType::typeName
    ),
    defaultValue_(
        pTraits<Type>::zero
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class T>
struct FoundValue {
    label procNo_;
    label cell_;
    point position_;
    T value_;
    bool found_;

    FoundValue()
        : procNo_(-1),
          cell_(-1),
          position_(vector::zero),
          value_(pTraits<T>::zero),
          found_(false)
        {}

    FoundValue(label c,const point &p)
        : procNo_(Pstream::myProcNo()),
          cell_(c),
          position_(p),
          value_(pTraits<T>::zero),
          found_(false)
        {}

    void operator()(const T &v) {
        value_=v;
        found_=true;
    }
};

template<typename T>
constexpr bool is_contiguous<FoundValue<T>> = contiguous<T>();

template<class T>
bool operator!=(const FoundValue<T> &a,const FoundValue<T> &b)
{
    return !(
        a.procNo_==b.procNo_
        &&
        a.cell_==b.cell_
    );
}

template<class T>
Istream & operator>>(Istream &in,FoundValue<T> &v)
{
    in >> v.procNo_;
    in >> v.cell_;
    in >> v.position_;
    in >> v.value_;
    in >> v.found_;
    return in;
}

template<class T>
Ostream & operator<<(Ostream &out,const FoundValue<T> &v)
{
    out << v.procNo_ << tab;
    out << v.cell_ << tab;
    out << v.position_ << tab;
    out << v.value_ << tab;
    out << v.found_ << tab;
    return out;
}

template<class T>
struct KeepFound {
    FoundValue<T> operator()(const FoundValue<T> &a, const FoundValue<T> &b) const {
        FoundValue<T> out = a;
        if (
            !a.found_
            &&
            b.found_
            &&
            a.cell_==b.cell_
            &&
            a.procNo_==b.procNo_
        ) {
            out(b.value_);
        }
        return out;
    }
};

template<class Type>
void findShiftedFieldPluginFunction<Type>::doEvaluation()
{
    autoPtr<ResultType> mappedField(
        new ResultType(
            default_()
        )
    );

    DynamicList<FoundValue<Type>> notFound;

    meshSearch search(mesh());

    forAll(shift_(),cellI) {
        const vector &shift=shift_()[cellI];
        const point targetPlace=mesh().C()[cellI]+shift;
        label fromCell=search.findCell(targetPlace);
        if (fromCell>=0) {
            mappedField()[cellI]=field_()[fromCell];
        } else {
            notFound.append(
                FoundValue<Type>(
                    cellI,
                    targetPlace
                )
            );
        }
    }
    notFound.shrink();
    Pbug << notFound.size() << " cells not found" << endl;

    if (Pstream::parRun()) {
        List<List<FoundValue<Type>>> data(Pstream::nProcs());
        data[Pstream::myProcNo()]=notFound;
        Pstream::allGatherList(data);
        label cnt=0;
        forAll(data,procI) {
            if (procI!=Pstream::myProcNo()) {
                forAll(data[procI],i) {
                    FoundValue<Type> &d=data[procI][i];
                    label fromCell=search.findCell(d.position_);
                    if (fromCell>=0) {
                        d(field_()[fromCell]);
                        cnt++;
                    }
                }
            }
        }
        Pbug << "Found " << cnt << endl;
        forAll(data,procI) {
            forAll(data[procI],i) {
                // Yikes.
                reduce(data[procI][i], KeepFound<Type>());
            }
        }
        cnt=0;
        forAll(data[Pstream::myProcNo()],i) {
            FoundValue<Type> &d=data[Pstream::myProcNo()][i];
            if (d.found_) {
                cnt++;
                mappedField()[d.cell_]=d.value_;
            }
        }
        Pbug << "Found " << cnt << " of " << notFound.size() << endl;
    }

    mappedField->correctBoundaryConditions();

    result().setObjectResult(
        mappedField
    );
}

template<class Type>
void findShiftedFieldPluginFunction<Type>::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0 || index==1 || index==2);

    if (index==0) {
        Dbug<< "Setting field" << endl;
        field_.set(
            new ResultType(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<ResultType>()
            )
        );
    } else if (index==1) {
        Dbug<< "Setting shift" << endl;
        shift_.set(
            new volVectorField(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<volVectorField>()
            )
        );
    } else if (index==2) {
        Dbug<< "Setting default" << endl;
        default_.set(
            new ResultType(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<ResultType>()
            )
        );
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


} // namespace

// ************************************************************************* //
