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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/TableFile/TableFile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableFile<Type>::TableFile
(
    const word& entryName,
    const dictionary& dict
)
:
    TableBase<Type>(entryName, dict),
    fName_(dict.lookup<fileName>("file")),
    restartDictName_(fName_.expand() + ".restart"),
    runTimeModifiable_(dict.lookupOrDefault<Switch>("runTimeModifiable", false)),
    readInterval_(dict.lookupOrDefault<label>("readInterval", 1)),
    x2old_(0),
    x1cache_(-1),
    x2cache_(-1)
{
    this->table_ = readTableFile();
    TableBase<Type>::check();

    readRestartDict();
}


template<class Type>
Foam::Function1Types::TableFile<Type>::TableFile(const TableFile<Type>& tbl)
:
    TableBase<Type>(tbl),
    fName_(tbl.fName_),
    restartDictName_(tbl.restartDictName_),
    runTimeModifiable_(tbl.runTimeModifiable_),
    readInterval_(tbl.readInterval_),
    x2old_(tbl.x2old_),
    x1cache_(tbl.x1cache_),
    x2cache_(tbl.x2cache_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableFile<Type>::~TableFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::TableFile<Type>::value(const scalar x) const
{
    updateTable();
    return TableBase<Type>::value(x);
}


template<class Type>
Type Foam::Function1Types::TableFile<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // only handles correctly integrate(0,t) and integrate(t-dt,t)
    if (runTimeModifiable_)
    {
        x1cache_ = x1;
        x2cache_ = x2;

        if (x1cache_ > 0)
        {
            return TableBase<Type>::integrate(x1,x2);
        }

        Type integrationValue
        (
            TableBase<Type>::integrate(x2old_, x2cache_) + integrationValueOld_
        );

        return integrationValue;
    }

    // else return default integration
    return TableBase<Type>::integrate(x1,x2);
}


template<class Type>
void Foam::Function1Types::TableFile<Type>::readRestartDict() const
{
    autoPtr<ISstream> isPtr
    (
        fileHandler().NewIFstream(restartDictName_)
    );
    ISstream& is = isPtr();

    if (is.good())
    {
        dictionary dict;
        is >> dict;
        integrationValueOld_ = dict.lookup<Type>("integrationValueOld");
        x2old_ = dict.lookup<scalar>("x2old");
    }
    else
    {
        integrationValueOld_ = TableBase<Type>::integrate(0,0);
        x2old_ = 0.0;
    }
}


template<class Type>
void Foam::Function1Types::TableFile<Type>::writeRestartDict()
const
{
    autoPtr<Ostream> osPtr(fileHandler().NewOFstream(restartDictName_));
    dictionary restartDict;
    restartDict.add("x2old",x2old_);
    restartDict.add("integrationValueOld",integrationValueOld_);
    restartDict.write(osPtr());
}


template<class Type>
Foam::List<Foam::Tuple2<Foam::scalar,Type>>
Foam::Function1Types::TableFile<Type>::readTableFile() const
{
    fileName expandedFile(fName_);
    autoPtr<ISstream> isPtr
    (
        fileHandler().NewIFstream(expandedFile.expand())
    );
    ISstream& is = isPtr();

    if (!is.good())
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open file." << exit(FatalIOError);
    }

    List<Tuple2<scalar,Type>> newTable;
    is  >> newTable;

    return newTable;
}


template<class Type>
void Foam::Function1Types::TableFile<Type>::updateTable() const
{
    if (runTimeModifiable_)
    {
        counter_++;
        if (counter_ % readInterval_ == 0)
        {
            List<Tuple2<scalar, Type>> newTable = readTableFile();

            if (newTable != this->table_)
            {
                // clean old table data to re-compute of demand
                this->tableSamplesPtr_.clear();
                this->interpolatorPtr_.clear();

                // cache old table integration value after
                // check if integrate has ever been called and starts from 0
                if (x1cache_ == 0)
                {
                    // add integration value from last update time to current
                    integrationValueOld_ +=
                        TableBase<Type>::integrate(x2old_,x2cache_);

                    // update lower integration bound (former upper)
                    x2old_ = x2cache_;

                    // write restart data to file
                    writeRestartDict();
                }

                this->table_ = newTable;
                TableBase<Type>::check();
            }
        }
    }
}


template<class Type>
void Foam::Function1Types::TableFile<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));

    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    os.writeEntry("file", fName_);
    if (runTimeModifiable_)
    {
        os.writeEntry("runTimeModifiable",runTimeModifiable_);
        os.writeEntry("readInterval",readInterval_);
    }

    os.endBlock() << flush;
}


// ************************************************************************* //
