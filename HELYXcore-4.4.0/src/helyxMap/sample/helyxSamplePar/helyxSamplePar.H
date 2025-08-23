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
    Foam::helyxSamplePar

Description
    Parallel version of the 'helyxSample' class.

SourceFile
    helyxSamplePar.C

\*---------------------------------------------------------------------------*/

#ifndef helyxSamplePar_H
#define helyxSamplePar_H

#include "sample/helyxSample/helyxSample.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class helyxSamplePar Declaration
\*---------------------------------------------------------------------------*/

class helyxSamplePar
:
    public helyxSample
{
public:

    // Public Data

        word sampleMethod_;

        DynamicList<point> bodySurfCenter_;


    // Constructors

        //- Default construct
        helyxSamplePar();

        //- Construct from Time and fvMesh
        helyxSamplePar(const Time* runTime, const fvMesh* mesh);


    //- Destructor
    virtual ~helyxSamplePar();


    // Member Functions

        void setOptions(const argList& args);

        void getDirNames(wordList& dnames, const word& pattern);

        void getMapTime();

        void combineFields();

        void combineFieldsBK();

        void constructKnn();

        void getBodySurface();

        void searchNodeField(gridField& grid, const word& funct = "field");
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Transfer a field to a specific processor
template<class T>
void interProcTrans
(
    DynamicList<T>& field,
    const DynamicList<T>& source,
    const word& name
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
