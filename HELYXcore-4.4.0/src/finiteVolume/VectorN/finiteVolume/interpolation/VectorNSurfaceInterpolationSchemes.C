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
    (c) 2010 Ivor Clifford

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "VectorN/Fields/VectorNFieldTypes.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "interpolation/surfaceInterpolation/schemes/reverseLinear/reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBaseSurfaceInterpolationScheme(Type)                              \
                                                                              \
defineNamedTemplateTypeNameAndDebug(surfaceInterpolationScheme<Type>, 0);     \
                                                                              \
defineTemplateRunTimeSelectionTable                                           \
(                                                                             \
    surfaceInterpolationScheme<Type>,                                         \
    Mesh                                                                      \
);                                                                            \
                                                                              \
defineTemplateRunTimeSelectionTable                                           \
(                                                                             \
    surfaceInterpolationScheme<Type>,                                         \
    MeshFlux                                                                  \
);


#define doMakeSchemes(type, Type, args...)                                    \
    makeBaseSurfaceInterpolationScheme(type)                                  \
                                                                              \
    makeSurfaceInterpolationTypeScheme(linear, type)                          \
    makeSurfaceInterpolationTypeScheme(reverseLinear, type)

forAllVectorNTypes(doMakeSchemes)

forAllTensorNTypes(doMakeSchemes)

forAllDiagTensorNTypes(doMakeSchemes)

forAllSphericalTensorNTypes(doMakeSchemes)

#undef doMakeSchemes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
