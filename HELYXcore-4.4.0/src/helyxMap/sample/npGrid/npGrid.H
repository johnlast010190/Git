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
    (c) 2019-2025 Engys Ltd.

Class
    Foam::npGrid

Description
   Class designed to convert a Numpy array into a HELYX field for visualization
   in Paraview.

SourceFiles
    npGrid.C

\*---------------------------------------------------------------------------*/

#ifndef npGrid_H
#define npGrid_H

#include <map>
#include "sample/gridField/gridField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class npGrid Declaration
\*---------------------------------------------------------------------------*/

class npGrid
:
    public gridField
{
public:

    // Public Data

        std::map<word, nodeField> scalarFields_;

        std::map<word, nodeVectors> vectorFields_;

        std::map<word, nodeField> vectorComps_;


    // Static Member Functions

        //- Read content of a text file
        static void readTextFile(word& input, const word& fname);


    // Constructors

        //- Default construct
        npGrid();

        //- Construct from words
        npGrid(const word& blkdef, const word& npfiles);


    //- Destructor
    virtual ~npGrid();


    // Member Functions

        void getScalarField(const word& fname, DynamicList<scalar>& scfield);

        void combineVectorField();

        void writeMeshDict();

        void writeScalarField(label k = 0);

        void writeVectorField(label k = 0);

        void writeVisualization();

        word timeName() const
        {
            return "1";
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
