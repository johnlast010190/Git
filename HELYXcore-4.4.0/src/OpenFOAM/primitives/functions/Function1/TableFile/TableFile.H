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
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::Function1Types::TableFile

Description
    Templated table container function where data is read from file.

    Usage:
    \verbatim
        <entryName> tableFile;
        <entryName>Coeffs
        {
            fileName            dataFile;    // name of data file
            outOfBounds         clamp;       // optional out-of-bounds handling
            interpolationScheme linear;      // optional interpolation method
        }
    \endverbatim

    Data is stored as a list of Tuple2's. First column is always stored as
    scalar entries.  Data is read in the form, e.g. for an entry \<entryName\>
    that is (scalar, vector):
    \verbatim
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        );
    \endverbatim

SourceFiles
    TableFile.C

\*---------------------------------------------------------------------------*/

#ifndef TableFile_H
#define TableFile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class TableFile Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableFile
:
    public TableBase<Type>
{
    // Private data

        //- File name for csv table (optional)
        fileName fName_;

        //- File name where restart data is read from
        fileName restartDictName_;

        //- Run-time modifiable
        bool runTimeModifiable_;
        label readInterval_;

        //- Store data from old table prior to re-read
        mutable Type integrationValueOld_;
        mutable scalar x2old_;
        mutable scalar x1cache_;
        mutable scalar x2cache_;

        mutable label counter_;


public:

    //- Runtime type information
    TypeName("tableFile");


    // Constructors

        //- Construct from entry name and Istream
        TableFile(const word& entryName, const dictionary& dict);

        //- Copy constructor
        TableFile(const TableFile<Type>& tbl);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const override
        {
            return tmp<Function1<Type>>(new TableFile<Type>(*this));
        }


    //- Destructor
    virtual ~TableFile();


    // I/O

        //- Re-read the dictionary at runtime
        List<Tuple2<scalar,Type>> readTableFile() const;
        void updateTable() const;
        void readRestartDict() const;
        void writeRestartDict() const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const override;


    // Member Functions

        //- Return Table value
        virtual Type value(const scalar x) const override;

        //- Integrate between two (scalar) values
        virtual Type integrate
        (
            const scalar x1,
            const scalar x2
        ) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TableFile<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/TableFile/TableFile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
