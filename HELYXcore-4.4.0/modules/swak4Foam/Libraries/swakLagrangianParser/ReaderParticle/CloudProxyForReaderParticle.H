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
    Foam::CloudProxyForReaderParticle

Description

SourceFiles
    CloudProxyForReaderParticle.C

Contributors/Copyright:
    2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CloudProxyForReaderParticle_H
#define CloudProxyForReaderParticle_H

#include "CloudProxyForParticle.H"

#include "ReaderParticleCloud.H"

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class CloudProxyForReaderParticle Declaration
\*---------------------------------------------------------------------------*/

class CloudProxyForReaderParticle
:
    public CloudProxyForParticle<ReaderParticleCloud>
{
    // Private Data

        const ReaderParticleCloud& readerCloud_;


    // Private Member Functions

        template <typename T>
        void internalAddFields(const wordList&);


public:

    //- Runtime type information
    TypeName("CloudProxyForReaderParticle");


    // Constructors

        //- Construct from components
        CloudProxyForReaderParticle(const cloud&);

        //- Disallow default bitwise copy construction
        CloudProxyForReaderParticle
        (
            const CloudProxyForReaderParticle&
        ) = delete;


    //- Destructor
    virtual ~CloudProxyForReaderParticle();


    // Public Member Functions

        virtual tmp<Field<scalar>> getScalarField(const word& name) const;

        virtual tmp<Field<vector>> getVectorField(const word& name) const;

        virtual tmp<Field<tensor>> getTensorField(const word& name) const;

        virtual tmp<Field<symmTensor>> getSymmTensorField
        (
            const word& name
        ) const;

        virtual tmp<Field<sphericalTensor>> getSphericalTensorField
        (
            const word& name
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CloudProxyForReaderParticle&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudProxyForReaderParticleI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
