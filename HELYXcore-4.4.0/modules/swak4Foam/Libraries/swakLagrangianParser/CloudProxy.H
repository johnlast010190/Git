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
    (c) 1991-2010 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::CloudProxy

Description
    Class that interfaces to the actual clouds

SourceFiles
    CloudProxy.C

Contributors/Copyright:
    2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CloudProxy_H
#define CloudProxy_H

#include "fields/cloud/cloud.H"
#include "fields/Fields/Field/Field.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class CloudProxy;

Ostream& operator<<(Ostream&, const CloudProxy&);

/*---------------------------------------------------------------------------*\
                         Class CloudProxy Declaration
\*---------------------------------------------------------------------------*/

class CloudProxy
{
    // Private data

        const cloud& cloud_;

        //- Types of the available functions
        HashTable<word, word> fieldTypes_;

        //- Descriptions of the functions
        HashTable<string, word> fieldDescriptions_;


protected:

    //- Add description of a "field"
    template<class Type>
    void addField(const word& name, const string& description);

    const cloud& theCloud() const
    {
        return cloud_;
    }


public:

    //- Runtime type information
    TypeName("CloudProxy");


    // Constructors

        //- Construct from components
        CloudProxy(const cloud&);

        //- Disallow default bitwise copy construction
        CloudProxy(const CloudProxy&) = delete;


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CloudProxy,
        cloud,
        (
            const cloud& aCloud
        ),
        (aCloud)
    );

    static autoPtr<CloudProxy> New
    (
        const cloud& c,
        const word& alternateType = ""
    );

    static autoPtr<CloudProxy> New
    (
        const cloud& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~CloudProxy();

    // Public Memeber Functions

    bool hasField(const word& name) const
    {
        return fieldTypes_.found(name);
    }

    word fieldType(const word& name) const
    {
        return fieldTypes_[name];
    }

    bool isScalarField(const word& name) const;

    bool isVectorField(const word& name) const;

    bool isTensorField(const word& name) const;

    bool isSymmTensorField(const word& name) const;

    bool isSphericalTensorField(const word& name) const;

    bool isBoolField(const word& name) const;

    virtual tmp<scalarField> getScalarField(const word& name) const = 0;

    virtual tmp<Field<vector>> getVectorField(const word& name) const = 0;

    virtual tmp<Field<tensor>> getTensorField(const word& name) const = 0;

    virtual tmp<Field<symmTensor>> getSymmTensorField
    (
        const word& name
    ) const = 0;

    virtual tmp<Field<sphericalTensor>> getSphericalTensorField
    (
        const word& name
    ) const = 0;

    virtual tmp<Field<bool>> getBoolField(const word& name) const = 0;

    virtual tmp<Field<vector>> getPositions() const = 0;

    virtual tmp<Field<label>> getCells() const = 0;

    virtual tmp<scalarField> weights() const = 0;

    virtual label size() const = 0;

    friend Ostream& operator<<(Ostream&, const CloudProxy&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CloudProxy&) = delete;
};


#define addCloudProxyToTable(ProxyType,CloudType)\
    CloudProxy::addcloudConstructorToTable                              \
    <ProxyType<CloudType>>                                             \
    add##ProxyType##CloudType##cloudConstructorToCloudProxyTable(       \
        Cloud<CloudType::particleType>::typeName)

#define addCloudProxyToTable2(ProxyType,CloudType,ParcelType)              \
    CloudProxy::addcloudConstructorToTable                              \
    <ProxyType<CloudType>>                                             \
    add##ProxyType##CloudType##cloudConstructorToCloudProxyTable(       \
        Cloud<ParcelType>::typeName)

    // This variation is used to add proxies that work in CloudFunctionObjects
#define addKinematicCloudProxyToTable(ProxyType,CloudType)      \
    CloudProxy::addcloudConstructorToTable                              \
    <ProxyType<CloudType::kinematicCloudType>>                         \
    addKinematic##ProxyType##CloudType##cloudConstructorToCloudProxyTable( \
        word("kinematic")+Cloud<CloudType::particleType>::typeName)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudProxyI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
