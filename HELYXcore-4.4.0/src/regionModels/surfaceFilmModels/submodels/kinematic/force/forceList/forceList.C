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
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/kinematic/force/forceList/forceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

forceList::forceList(surfaceFilmRegionModel& film)
:
    PtrList<force>()
{}


forceList::forceList
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    PtrList<force>()
{
    Info<< "    Selecting film force models" << endl;

    if (dict.isDict("forces"))
    {
        const dictionary& forcesDict = dict.subDict("forces");
        this->setSize(forcesDict.size());

        label i = 0;
        forAllConstIter(dictionary, forcesDict, iter)
        {
            set
            (
                i++,
                force::New
                (
                    film,
                    forcesDict.isDict(iter().keyword())
                  ? forcesDict.subDict(iter().keyword())
                  : dictionary::null,
                    iter().keyword()
                )
            );
        }
    }
    else if (dict.found("forces"))
    {
        const wordList models(dict.lookup("forces"));

        this->setSize(models.size());

        forAll(models, i)
        {
            set(i, force::New(film, dict, models[i]));
        }
    }

    if (!size())
    {
        Info<< "        none" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceList::~forceList()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> forceList::correct(volVectorField& U)
{
    tmp<fvVectorMatrix> tResult(new fvVectorMatrix(U, dimForce));
    fvVectorMatrix& result = tResult.ref();

    forAll(*this, i)
    {
        result += this->operator[](i).correct(U);
    }

    return tResult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
