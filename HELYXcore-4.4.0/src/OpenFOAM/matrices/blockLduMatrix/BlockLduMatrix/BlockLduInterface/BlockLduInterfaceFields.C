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
    (c) 2004-2006 Ivor Clifford

Class
    BlockLduInterfaceField

Description
    Abstract base class for interface fields with block coefficients

Author
    Ivor Clifford, ivor.clifford@gmail.com

\*----------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduInterface/BlockLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(BlockLduInterfaceField<scalar>, 0);
defineTemplateTypeNameAndDebug(BlockLduInterfaceField<vector>, 0);
defineTemplateTypeNameAndDebug(BlockLduInterfaceField<tensor>, 0);
defineTemplateTypeNameAndDebug(BlockLduInterfaceField<sphericalTensor>, 0);
defineTemplateTypeNameAndDebug(BlockLduInterfaceField<symmTensor>, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
