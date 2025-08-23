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
    (c) 2024 Engys Ltd.

Class
    Foam::ReaderParticleCloud

Description
    A Cloud of passive particles

SourceFiles
    ReaderParticleCloud.C

Contributors/Copyright:
    2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef ReaderParticleCloud_H
#define ReaderParticleCloud_H

#include "include/swakCloudTypes.H"

#ifdef FOAM_CLOUDTEMPLATE_IS_CLOUD_H
#include "CloudTemplate.H"
#else
#include "Cloud/Cloud.H"
#endif

#include "ReaderParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class CloudRepository;

/*---------------------------------------------------------------------------*\
                   Class ReaderParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class ReaderParticleCloud
:
    public Cloud<ReaderParticle>
{
    // Private Data

        HashTable<Field<label>,word> labelData_;
        HashTable<Field<scalar>,word> scalarData_;
        HashTable<Field<vector>,word> vectorData_;
        HashTable<Field<tensor>,word> tensorData_;
        HashTable<Field<symmTensor>,word> symmTensorData_;
        HashTable<Field<sphericalTensor>,word> sphericalTensorData_;


    // Private Member Functions

        template <typename T>
        void addDataInternal
        (
            const word& name,
            const Field<T>& inData,
            HashTable<Field<T>,word>& data,
            bool overwrite=false
        );


protected:

    void addScalarField(const word& name);
    void addLabelField(const word& name);
    void addVectorField(const word& name);
    void addTensorField(const word& name);
    void addSymmTensorField(const word& name);
    void addSphericalTensorField(const word& name);

    template <typename T>
    void writeData(const HashTable<Field<T>,word>& data) const;

    void clearData();


public:

    //- Type of parcel the cloud was instantiated for
    typedef ReaderParticle particleType;

    // Constructors

        //- Construct given mesh
        ReaderParticleCloud
        (
            const polyMesh&,
            const word& cloudName,
            bool readFields = true
        );

        //- Disallow default bitwise copy construction
        ReaderParticleCloud(const ReaderParticleCloud&) = delete;


    static autoPtr<ReaderParticleCloud> makeCloudFromVariables
    (
        const polyMesh& mesh,
        const word& cloudName,
        const wordList& globalNameSpacesToUse,
        const word& positionVar
    );

    friend class ReaderParticle;

    template <typename T>
    void setValues
    (
        const word& name,
        const Field<T>& values
    );

    template <typename T>
    const Field<T>& getValues(const word& name);

    wordList getScalarFieldNames() const
    {
        return scalarData_.toc();
    }

    const Field<scalar>& getScalarField(const word& name) const
    {
        return scalarData_[name];
    }

    bool hasScalar(const word& name) const
    {
        return scalarData_.found(name);
    }

    wordList getLabelFieldNames() const
    {
        return labelData_.toc();
    }

    const Field<label>& getLabelField(const word& name) const
    {
        return labelData_[name];
    }

    bool hasLabel(const word& name) const
    {
        return labelData_.found(name);
    }

    wordList getVectorFieldNames() const
    {
        return vectorData_.toc();
    }

    const Field<vector>& getVectorField(const word& name) const
    {
        return vectorData_[name];
    }

    bool hasVector(const word& name) const
    {
        return vectorData_.found(name);
    }

    wordList getTensorFieldNames() const
    {
        return tensorData_.toc();
    }

    const Field<tensor>& getTensorField(const word& name) const
    {
        return tensorData_[name];
    }

    bool hasTensor(const word& name) const
    {
        return tensorData_.found(name);
    }

    wordList getSymmTensorFieldNames() const
    {
        return symmTensorData_.toc();
    }

    const Field<symmTensor>& getSymmTensorField(const word& name) const
    {
        return symmTensorData_[name];
    }

    bool hasSymmTensor(const word& name) const
    {
        return symmTensorData_.found(name);
    }

    wordList getSphericalTensorFieldNames() const
    {
        return sphericalTensorData_.toc();
    }

    const Field<sphericalTensor>& getSphericalTensorField
    (
        const word& name
    ) const
    {
        return sphericalTensorData_[name];
    }

    bool hasSphericalTensor(const word& name) const
    {
        return sphericalTensorData_.found(name);
    }

    friend class CloudRepository;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ReaderParticleCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReaderParticleCloudI.H"
#endif

#endif

// ************************************************************************* //
