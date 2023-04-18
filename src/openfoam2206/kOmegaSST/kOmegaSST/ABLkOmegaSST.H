/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::ABLkOmegaSST

Description
    Specialisation for RAS of the generic ABLkOmegaSSTBase base class.
    For more information, see Description of ABLkOmegaSSTBase.H

See also
    Foam::ABLkOmegaSST

SourceFiles
    ABLkOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef ABLkOmegaSST_H
#define ABLkOmegaSST_H

#include "ABLkOmegaSSTBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class ABLkOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class ABLkOmegaSST
:
    public Foam::ABLkOmegaSST
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
{
/*protected:

    // Protected data

        // Model coefficients

            dimensionedScalar betaStarMax_;

        // ABL model coefficients

        dimensionedScalar z0_;     // Aerodynamic roughness
        dimensionedScalar kappa_;  // von Karman constant
        dimensionedScalar ustar_;  // friction velocity
        dimensionedVector zDir_; // Normal of the wall
        dimensionedVector flowDir_;  // Direction of the reference flow
        volScalarField centres_;
        volScalarField uref_;
        volScalarField kref_;
        volScalarField oref_;

        volScalarField Uerr_;
        volScalarField Kerr_;
        volScalarField Oerr_;
        volScalarField Hybriderr_;
        volScalarField BlendingFun_;

        dimensionedScalar A_;
        dimensionedScalar B_;
        dimensionedScalar C_;
        dimensionedScalar D_;
        dimensionedScalar E_;
        dimensionedScalar blendCoeff_;
    
        dimensionedScalar thresholdU_;
        dimensionedScalar thresholdK_;
        dimensionedScalar thresholdO_;
        dimensionedScalar constantU_;
        dimensionedScalar constantK_;
        dimensionedScalar constantO_;
    
    
        dimensionedVector msel_;
    
        dimensionedScalar usmall;
        dimensionedScalar ksmall;
        dimensionedScalar osmall;
    
        volScalarField NLbetaStar_;
        volScalarField hombetaStar_;
        volScalarField blendbetaStar_;


        // Fields

        tmp<volScalarField> uError();
        tmp<volScalarField> kError();
        tmp<volScalarField> oError();
        tmp<volScalarField> hybridError();
        tmp<volScalarField> homogeneousbetaStar();
        tmp<volScalarField> nonlinearbetaStar();
        tmp<volScalarField> blendedbetaStar();

    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;

        tmp<volScalarField> blendingFunction();
        void errorsAndbetaStarComputation();
        void referencesComputation();*/


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("ABLkOmegaSST");


    // Constructors

        //- Construct from components
        ABLkOmegaSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~ABLkOmegaSST()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "ABLkOmegaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
