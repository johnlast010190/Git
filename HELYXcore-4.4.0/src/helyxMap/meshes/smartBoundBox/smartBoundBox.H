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
    (c) 2020-2025 Engys Ltd.

Class
    Foam::smartBoundBox

Description
    Smart bounding box with irregular boundary. The purpose of this class is
    to efficiently determine whether a given point (x,y,z) is within the
    boundBox. The determination is approximate, to the error of maximum cell
    dimension*ampCoef, where 'ampCoef' is the amplification coefficient.

SourceFiles
    smartBoundBox.C

\*---------------------------------------------------------------------------*/

#ifndef smartBoundBox_H
#define smartBoundBox_H

#include "containers/Array3d/Array3d.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class smartBoundBox Declaration
\*---------------------------------------------------------------------------*/

class smartBoundBox
:
    public boundBox
{
protected:

    // Protected Data

        //- Type of bounding box, simple = 0, smart = 1.
        label boxType_;

        //- Number of cells in a grid
        Array3d<bool> numCells_;

        //- Biggest cell dimension amplification coefficient
        scalar ampCoef_;

        //- Grid dimension in x
        label nx_;

        //- Grid dimension in y
        label ny_;

        //- Grid dimension in z
        label nz_;

        scalar deltaX_;
        scalar deltaY_;
        scalar deltaZ_;


public:

    // Constructors

        //- Default construct
        smartBoundBox();

        smartBoundBox
        (
            const fvMesh* mesh,
            const word type = "boundBox",
            scalar ampCoef = 1.5,
            label cellNum = 40000
        );

        smartBoundBox
        (
            const point& pmin,
            const point& pmax,
            const word type = "boundBox",
            scalar ampCoef = 1.5
        );


    // Member Functions

        void createBox(const fvMesh *mesh, label ncells = 40000);

        inline bool pointInBox(const point& pt) const;

        //- Get number of cells in each grid
        void getCellList(const fvMesh* mesh, label);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "smartBoundBoxI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
