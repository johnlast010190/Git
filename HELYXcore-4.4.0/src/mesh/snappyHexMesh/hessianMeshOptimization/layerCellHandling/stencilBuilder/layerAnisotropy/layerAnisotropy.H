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
    Optimize anisotropically for layer cells

SourceFiles
    layerAnisotropy.C

\*---------------------------------------------------------------------------*/

#ifndef layerAnisotropy_H
#define layerAnisotropy_H

#include "hessianMeshOptimization/layerCellHandling/stencilBuilder/stencilBuilder.H"
#include "hessianMeshOptimization/layerCellHandling/layerCellDetection.H"

namespace Foam
{

class layerAnisotropy : public stencilBuilder
{

    const fvMesh& mesh_;

    layerCellDetection layers_;

    scalarField anisotropy_;

    layerAnisotropy(const layerAnisotropy&);

    void operator=(const layerAnisotropy&);

    void calculateAnisotropy();

public:

    //- Runtime type information
    TypeName("layerAnisotropy");

    // Constructors

        layerAnisotropy
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    virtual ~layerAnisotropy();

    const scalarField& anisotropy() const;

    const vectorField& anisotropicDirection() const;

};
} //End namespace Foam

#endif
