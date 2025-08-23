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
    (c) 2014 Engys Ltd.

Class
    Foam::isotropy

Description
    Optimize without considering anisotropic cells

SourceFiles
    isotropy.C

\*---------------------------------------------------------------------------*/

#ifndef isotropy_H
#define isotropy_H

#include "hessianMeshOptimization/layerCellHandling/stencilBuilder/stencilBuilder.H"

namespace Foam
{

class isotropy : public stencilBuilder
{
    scalarField anisotropy_;

    vectorField anisotropicDirection_;

    isotropy(const isotropy&);

    void operator=(const isotropy&);

public:

    //- Runtime type information
    TypeName("isotropy");

    // Constructors

        isotropy
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    virtual ~isotropy();

    const scalarField& anisotropy() const;

    const vectorField& anisotropicDirection() const;
};

} // End namespace Foam

#endif
