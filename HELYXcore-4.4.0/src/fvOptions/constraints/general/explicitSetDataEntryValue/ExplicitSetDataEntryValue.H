/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017 Engys Ltd.

Class
    Foam::fv::explicitSetDataEntryValue

Description
    Set values field values explicity.

    Sources described by:

        \verbatim
        <Type>ExplicitSetDataEntryValueCoeffs
        {
            sourceData
            {
                p table
                (
                    (0 0)
                    (1 10)
                    (2 0)
                );
                epsilon  constant   1.5;
            }
        }
        \endverbatim

SeeAlso
    Foam::fvOption

SourceFiles
    explicitSetDataEntryValue.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitSetDataEntryValue_H
#define ExplicitSetDataEntryValue_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "primitives/Tuple2/Tuple2.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class explicitSetDataEntryValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ExplicitSetDataEntryValue
:
    public cellSetOption
{

protected:

    // Protected data

        //- Source value per field
        PtrList<Function1<Type>> sourceData_;


    // Protected functions

        //- Set the local field data
        void setFieldData(const dictionary& dict, wordList& fieldNames);


public:

    //- Runtime type information
    TypeName("explicitSetDataEntryValue");


    // Constructors

        //- Construct from components
        ExplicitSetDataEntryValue
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);


        // Evaluation

            //- Set value on field
            virtual void setValue(fvMatrix<Type>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ExplicitSetDataEntryValue.C"
#   include "ExplicitSetDataEntryValueIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
