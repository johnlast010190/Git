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
    timelineFunctionObject

Description
    Basis for function objects that generate timelines

SourceFiles
    timelineFunctionObject.C

Contributors/Copyright:
    2008-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2017 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef timelineFunctionObject_H
#define timelineFunctionObject_H

#include "simpleDataFunctionObject/simpleDataFunctionObject.H"
#include "primitives/enums/NamedEnum.H"
#include "primitives/strings/lists/stringList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class timelineFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class timelineFunctionObject
:
    public simpleDataFunctionObject
{
    // Private Data

        //- Factor to multiply results with (to account for symmetries etc)
        scalar factor_;


public:

    // Public Data

        enum outputFileMode
        {
            ofmFoam,
            ofmRaw,
            ofmCsv
        };


private:

    // Private Data

        //- Current open files
        HashPtrTable<OFstream> filePtrs_;


protected:

    // Protected Data

        //- Possible output-modes
        static const NamedEnum<outputFileMode, 3> outputFileModeNames_;
        outputFileMode outputFileMode_;

        //- The file extension that fits the file-mode
        word fileExtension_;

        //- Character that separates fields
        std::string separator_;


    // Protected Member Functions

        //- Name of the sub-dictionary of the case
        virtual word dirName() = 0;

        //- Names of the files
        virtual wordList fileNames() = 0;

        //- Text in the first line
        virtual stringList columnNames() = 0;

        //- Factor to multiply results with (for symmetries etc)
        scalar factor() const
        {
            return factor_;
        }

        //- Write similar data items
        template<class T>
        void writeData(const word& name, const UList<T>& data);

        //- Write the time
        template<class T>
        void writeTime(const word& name, const T& time);

        //- Ends a dataset (usually by writing an endline)
        void endData(const word& name);

        //- Close all files to enable a new slate
        void closeAllFiles();


private:

    // Private Data

        //- Should we write immediately at startup
        bool writeStartTime_;


public:

    //- Runtime type information
    TypeName("timelineFunctionObject");


    // Constructors

        //- Construct from components
        timelineFunctionObject
        (
            const word &name,
            const Time&,
            const dictionary&,
            bool writeStartTime = false
        );

        //- Disallow default bitwise copy construct
        timelineFunctionObject(const timelineFunctionObject&) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        virtual void flush();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const timelineFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "timelineFunctionObjectI.H"
#endif

#endif

// ************************************************************************* //
