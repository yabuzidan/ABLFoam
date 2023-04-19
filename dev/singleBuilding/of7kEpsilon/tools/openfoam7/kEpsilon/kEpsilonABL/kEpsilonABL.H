/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::RASModels::kEpsilonABL

Description
    Standard k-epsilon turbulence model for incompressible and compressible
    flows including rapid distortion theory (RDT) based compression term.

    Reference:
    \verbatim
        Standard model:
            Launder, B. E., & Spalding, D. B. (1972).
            Lectures in mathematical models of turbulence.

            Launder, B. E., & Spalding, D. B. (1974).
            The numerical computation of turbulent flows.
            Computer methods in applied mechanics and engineering,
            3(2), 269-289.

        For the RDT-based compression term:
            El Tahry, S. H. (1983).
            k-epsilon equation for compressible reciprocating engine flows.
            Journal of Energy, 7(4), 345-353.
    \endverbatim

    The default model coefficients are
    \verbatim
        kEpsilonABLCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    kEpsilonABL.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonABL_H
#define kEpsilonABL_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "wallFvPatch.H"
#include "wallDist.H"
#include "wallDistData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonABL Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonABL
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
protected:

    // Protected data

        // Model coefficients

        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar C3_;
        dimensionedScalar sigmak_;
        dimensionedScalar sigmaEps_;
        dimensionedScalar CmuMax_;
    
        // ABL model coefficients

        dimensionedScalar z0_;     // Aerodynamic roughness
        dimensionedScalar kappa_;  // von Karman constant
        dimensionedScalar ustar_;  // friction velocity
        dimensionedVector zDir_; // Normal of the wall
        dimensionedVector flowDir_;  // Direction of the reference flow
        volScalarField centres_;
        volScalarField epsSource_;
        volScalarField kSource_;
    
        volScalarField uref_;
        volScalarField kref_;
        volScalarField eref_;
        volScalarField nref_;
    
        volScalarField Uerr_;
        volScalarField Kerr_;
        volScalarField Eerr_;
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
        dimensionedScalar thresholdE_;
        dimensionedScalar constantU_;
        dimensionedScalar constantK_;
        dimensionedScalar constantE_;
    
    
        dimensionedVector msel_;
    
        dimensionedScalar usmall;
        dimensionedScalar ksmall;
        dimensionedScalar esmall;
    
        volScalarField NLCmu_;
        volScalarField homCmu_;
        volScalarField blendCmu_;
        volScalarField wSwitch_;
    
        volScalarField Gref_;
        volScalarField eterm1_;
        volScalarField eterm2_;
        volScalarField kterm1_;
        volScalarField kterm2_;

        // Fields

        volScalarField k_;
        volScalarField epsilon_;
    

        tmp<volScalarField> uError();
        tmp<volScalarField> kError();
        tmp<volScalarField> eError();
        tmp<volScalarField> hybridError();
        tmp<volScalarField> homogeneousCmu();
        tmp<volScalarField> nonlinearCmu();
        tmp<volScalarField> blendedCmu();


    // Protected Member Functions

        virtual void correctNut();
    
        tmp<volScalarField> ABLepsilonSource();
        tmp<volScalarField> ABLkSource();
        tmp<volScalarField> blendingFunction();
        void errorsAndCmuComputation();
        void referencesComputation();




public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonABL");


    // Constructors

        //- Construct from components
        kEpsilonABL
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

        //- Disallow default bitwise copy construction
        kEpsilonABL(const kEpsilonABL&) = delete;


    //- Destructor
    virtual ~kEpsilonABL()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kEpsilonABL&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonABL.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
