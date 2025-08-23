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
    (c) 2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.
    (c) 2020-2024 Engys Ltd.

Class
    Foam::functionObjects::histogram

Group
    grpFieldFunctionObjects

Description
    Write the volume-weighted histogram of a volScalarField.

Usage
    Example:
    \verbatim
    histogram1
    {
        type            histogram;

        libs            ("libfieldFunctionObjects.so");

        field           p;
        nBins           100;
        min             -5;
        max             5;
        setFormat       gnuplot;
        sets
        (
            boxToCell {box (0 0 0)(1 1 1);}
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: histogram    | yes         |
        field        | Field to analyse        | yes         |
        nBins        | Number of bins for the histogram | yes|
        max          | Maximum value sampled   | no          | field max
        min          | minimum value sampled   | no          | 0
        setFormat    | Output format           | yes         |
        sets         | Output over cellSet     | no          |
    \endtable

Note
    If max is not provided it will use the field's min and max as the bin
    extremes. If max is provided but not min it will use 0. The set written
    contains two columns, the first the volume averaged values, the second
    the raw bin count.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    histogram.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_histogram_H
#define functionObjects_histogram_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "sampledSetWriters/writer.H"
#include "containers/Lists/PackedList/PackedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class histogram Declaration
\*---------------------------------------------------------------------------*/

class histogram
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Name of field
        word fieldName_;

        //- Maximum value
        scalar max_;

        //- Minimum value
        scalar min_;

        //- Mumber of bins
        label nBins_;

        //- Output formatter to write
        autoPtr<writer<scalar>> formatterPtr_;

        //- Cells to evaluate
        PackedList<1> reportCells_;


    // Private Member Functions

        void writeGraph
        (
            const coordSet& coords,
            const word& valueName,
            const scalarField& normalizedValues,
            const scalarField& absoluteValues
        ) const;


public:

    //- Runtime type information
    TypeName("histogram");


    // Constructors

        //- Construct from Time and dictionary
        histogram
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        histogram(const histogram&) = delete;


    //- Destructor
    virtual ~histogram();


    // Member Functions

        //- Read the histogram data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the histogram and write.
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const histogram&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
