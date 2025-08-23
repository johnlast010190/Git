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
    Foam::fv::matrixChangeBaseFvOption

Description
    This fvOption is the basis for two fvOptions that record the change of the
    matrix diagonal/source due to other fvOptions

SourceFiles
    matrixChangeBaseFvOption.C

Contributors/Copyright:
    2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef matrixChangeBaseFvOption_H
#define matrixChangeBaseFvOption_H

#include "memory/autoPtr/autoPtr.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "db/functionObjects/functionObjectList/functionObjectList.H"

#include "include/swakFvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class matrixChangeBaseFvOption Declaration
\*---------------------------------------------------------------------------*/

class matrixChangeBaseFvOption
:
    public option
{
    // Private data

        word fieldName_;
        word namePrefix_;
        bool doAtAddSup_;


protected:

    // Protected Member Functions

        template<class T>
        inline tmp<VolField<T>> calcResiduum(const fvMatrix<T> &) const;

        word fieldName() const
        {
            return namePrefix_ + fieldName_;
        }

        bool doAtAddSup()
        {
            return doAtAddSup_;
        }


public:

    //- Runtime type information
    TypeName("matrixChangeBase");


    // Constructors

        //- Construct from explicit source name and mesh
        matrixChangeBaseFvOption
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        matrixChangeBaseFvOption(const matrixChangeBaseFvOption&) = delete;


    //- Destructor
    ~matrixChangeBaseFvOption()
    {};

    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const matrixChangeBaseFvOption&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

#include "matrixChangeBaseFvOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
