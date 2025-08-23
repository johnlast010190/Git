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
    (c) 2010-2024 Engys Ltd.
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2017 OpenFOAM Foundation

Class
    Foam::fv::effectivenessHeatExchangerSource

Group
    grpFvOptionsSources

Description
    Heat exchanger source model, in which the heat exchanger is defined as a
    selection of cells.

    The total heat exchange source is given by:
    \f[
        Q_t = e(\phi, \dot{m}_2) (T_2 - T_1) \phi c_p
    \f]

    where:
    \vartable
        Q_t       | total heat source
        e(\phi,\dot{m}_2) | effectivenes table
        \phi      | net mass flux entering heat exchanger [kg/s]
        \dot{m}_2 | secondary mass flow rate [kg/s]
        T_1       | primary inlet temperature [K]
        T_2       | secondary inlet temperature [K]
        c_p       | specific heat capacity [J/kg/K]
    \endvartable


    The distribution inside the heat exchanger is given by:
    \f[
        Q_c = \frac{V_c |U_c| (T_c - T_{ref})}{\sum(V_c |U_c| (T_c - T_{ref}))}
    \f]

    where:
    \vartable
        Q_c     | source for cell
        V_c     | volume of the cell [m3]
        U_c     | local cell velocity [m/s]
        T_c     | local call temperature [K]
        T_{ref} | min or max(T) in cell zone depending on the sign of Q_t [K]
    \endvartable

Usage
    Example usage:
    \verbatim
    effectivenessHeatExchangerSource1
    {
        type            effectivenessHeatExchangerSource;
        active          yes;

        selectionMode   cellZone;
        cellZone        porosity;

        secondaryMassFlowRate   1.0;
        secondaryInletT 336;
        secondaryCp  4.18; //units of J/kgK
        faceZone        facesZoneInletOriented;
        outOfBounds     clamp;
        file            "effTable";
    }
    \endverbatim

    The effectiveness table is described in terms of the primary and secondary
    mass flow rates.  For example, the table:

    \verbatim
                           secondary MFR
                       |  0.1   0.2   0.3
                  -----+-----------------
                  0.02 |   A     B     C
     primary MFR  0.04 |   D     E     F
                  0.06 |   G     H     I
    \endverbatim


    Is specified by the following:

    \verbatim
        (
            0.02
            (
                (0.1    A)
                (0.2    B)
                (0.3    C)
            ),
            0.04
            (
                (0.1    D)
                (0.2    E)
                (0.3    F)
            ),
            0.06
            (
                (0.1    G)
                (0.2    H)
                (0.3    I)
            )
        );
    \endverbatim

Note
- the table with name "file" should have the same units as the
  secondary mass flow rate and kg/s for phi
- faceZone is the faces at the inlet of the cellzone, it needs to be
  created with flip map flags. It is used to integrate the net mass flow
  rate into the heat exchanger


SourceFiles
    effectivenessHeatExchangerSource.C

\*---------------------------------------------------------------------------*/

#ifndef effectivenessHeatExchangerSource_H
#define effectivenessHeatExchangerSource_H

#include "cellSetOption/cellSetOption.H"
#include "memory/autoPtr/autoPtr.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class effectivenessHeatExchangerSource Declaration
\*---------------------------------------------------------------------------*/

class effectivenessHeatExchangerSource
:
    public cellSetOption
{
protected:

    // Protected data

        //- Secondary flow mass rate [kg/s]
        scalar secondaryMassFlowRate_;

        //- Inlet secondary temperature [K]
        scalar secondaryInletT_;

        //- Inlet secondary Cp [J/kg/K]
        scalar secondaryCp_;

        //- Primary air temperature at the heat exchanger inlet [K]
        mutable scalar primaryInletT_;

        //- Flag to use a user-specified primary inlet temperature
        bool userPrimaryInletT_;

        //- 2D look up table efficiency = function of primary and secondary
        // mass flow rates [kg/s]
        autoPtr<interpolation2DTable<scalar>> eTable_;

        //- Name of velocity field; default = U
        word UName_;

        //- Name of temperature field; default = T
        word TName_;

        //- Name of the flux
        word phiName_;

        //- Name of the faceZone at the heat exchange inlet
        word faceZoneName_;

        //- Face zone ID
        label zoneID_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- List of +1/-1 representing face flip map (1 use as is, -1 negate)
        labelList faceSign_;

        //- Inlet face area
        mutable scalar faceZoneArea_;

        //- Local list of cells upwind from inlet
        labelList upwindCells_;

        //- Volume of inlet cells for primary Inlet T averages
        scalar inletCellVolume_;

        //- Zone boundary faces for mass flow calc
        labelList zoneBoundaryFaces_;


private:

    // Private Member Functions

        //- Initialise heat exchanger source model
        void init();

        //- Calculate total area of faceZone accross processors
        void calculateTotalArea(scalar& area);

        //- cellZone boundary faces
        labelList zoneBoundaryFaces() const;

        //- determine cell ID upwind from inlet
        void determineUpwindCells();

        //- Return rho if rhoName is specified otherwise rhoRef
        tmp<volScalarField> getRho() const;

        //- return sum of mag of face flux
        scalar zoneMassFlux(const fvMesh& mesh, const labelList& faces);


public:

    //- Runtime type information
    TypeName("effectivenessHeatExchangerSource");


    // Constructors

        //- Construct from components
        effectivenessHeatExchangerSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        effectivenessHeatExchangerSource
        (
            const effectivenessHeatExchangerSource&
        ) = delete;


    //- Destructor
    virtual ~effectivenessHeatExchangerSource()
    {}


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;


        // Explicit and implicit source for ico equation

            //- Scalar
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;


        // Explicit and implicit source for compressible equation

            //- Scalar
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const effectivenessHeatExchangerSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
