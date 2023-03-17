/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "kEpsilonABL.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kEpsilonABL<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = blendCmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

}
    
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::ABLepsilonSource()
{
    bool debug=0;
    if (debug)
        Info << "Entered to ABLepsilonSource" << endl;
    
    return eterm1_ - eterm2_;
}
    
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::ABLkSource()
{
    bool debug=0;
    if (debug)
        Info << "Entered to ABLkSource" << endl;

    return -(kterm1_ + wSwitch_*kterm2_);
}

template<class BasicMomentumTransportModel>
void kEpsilonABL<BasicMomentumTransportModel>::referencesComputation()
{
    bool debug=0;
    if (debug)
        Info << "Entered to referencesComputation" << endl;
    
    // Set reference values for internal fields
    volScalarField coord = (centres_ + z0_)/z0_;
    volScalarField coord2 = (centres_ + z0_);
    coord2.dimensions().reset(dimless);
    
    uref_ = ustar_*log(coord)/kappa_;
    kref_ = A_*log(coord) + B_*sqr(coord) + C_*coord + D_ + E_*log(coord2);
    
    homCmu_ = homogeneousCmu();
    
    eref_ = sqrt(homCmu_)*kref_*ustar_/(kappa_*(centres_+z0_));
    nref_ = homCmu_*sqr(kref_)/eref_;
    
    Gref_   = pow(ustar_,3)/(kappa_*(centres_+z0_));
    eterm1_ = (pow(ustar_,4)/pow((centres_+z0_),2))*((C2_-C1_)*sqrt(homCmu_)/pow(kappa_,2));
    eterm2_ = (pow(ustar_,4)/pow((centres_+z0_),2))*(1/sigmaEps_);
    //eterm1_ = (pow(ustar_,4)/sigmaEps_)*(1/pow((centres_+z0_),2));
    //eterm2_ = C1_*Gref_*eref_/kref_ - C2_*eref_*eref_/kref_;
    kterm1_ = (ustar_*kappa_/(z0_*sigmak_))*(4*B_*(centres_+z0_)/z0_ + C_);
    kterm2_ = Gref_ - eref_;
    
    volScalarField mynu = this->nu();

    // Boundary fields
    const fvPatchList& patches = this->mesh_.boundary();
    volScalarField::Boundary& centresBf_ = this->centres_.boundaryFieldRef();
    volScalarField::Boundary& krefBf_ = this->kref_.boundaryFieldRef();
    volScalarField::Boundary& urefBf_ = this->uref_.boundaryFieldRef();
    volScalarField::Boundary& homCmuBf_ = this->homCmu_.boundaryFieldRef();
    volScalarField::Boundary& erefBf_ = this->eref_.boundaryFieldRef();
    volScalarField::Boundary& nrefBf_ = this->nref_.boundaryFieldRef();
    volScalarField::Boundary& GrefBf_ = this->Gref_.boundaryFieldRef();
    volScalarField::Boundary& eterm1Bf_ = this->eterm1_.boundaryFieldRef();
    volScalarField::Boundary& eterm2Bf_ = this->eterm2_.boundaryFieldRef();
    volScalarField::Boundary& kterm1Bf_ = this->kterm1_.boundaryFieldRef();
    volScalarField::Boundary& kterm2Bf_ = this->kterm2_.boundaryFieldRef();
    

    scalar coordt=0;
    scalar coord2t=0;
    
    forAll(patches, patchi)
    {
        fvPatchScalarField& pcentres_  = centresBf_[patchi];
        fvPatchScalarField& pkref_  = krefBf_[patchi];
        fvPatchScalarField& puref_  = urefBf_[patchi];
        fvPatchScalarField& phomCmu_  = homCmuBf_[patchi];
        fvPatchScalarField& peref_  = erefBf_[patchi];
        fvPatchScalarField& pnref_  = nrefBf_[patchi];
        fvPatchScalarField& pGref_  = GrefBf_[patchi];
        fvPatchScalarField& peterm1_  = eterm1Bf_[patchi];
        fvPatchScalarField& peterm2_  = eterm2Bf_[patchi];
        fvPatchScalarField& pkterm1_  = kterm1Bf_[patchi];
        fvPatchScalarField& pkterm2_  = kterm2Bf_[patchi];

        const fvPatch& curPatch = patches[patchi];
        
        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                wSwitch_[faceCelli] = 0.0;
            }
        }
        
        forAll(curPatch, facei)
        {

            coordt = (pcentres_[facei] + z0_.value())/z0_.value();
            coord2t = (pcentres_[facei] + z0_.value());
            
            puref_[facei] = ustar_.value()*log(coordt)/kappa_.value();
            pkref_[facei] = A_.value()*log(coordt) + B_.value()*sqr(coordt) + C_.value()*coordt + D_.value() +E_.value()*log(coord2t);
            
            peref_[facei] = sqrt(phomCmu_[facei])*pkref_[facei]*ustar_.value()/(kappa_.value()*(pcentres_[facei] + z0_.value()));

            pnref_[facei] = phomCmu_[facei]*sqr(pkref_[facei])/peref_[facei];

            pGref_[facei]   = pow(ustar_.value(),3)/(kappa_.value()*coord2t);
            
            peterm1_[facei] = (pow(ustar_.value(),4)/pow((coord2t),2))*((C2_.value()-C1_.value())*sqrt(phomCmu_[facei])/pow(kappa_.value(),2));
	        peterm2_[facei] = (pow(ustar_.value(),4)/pow((coord2t),2))*(1/sigmaEps_.value());
            //peterm1_[facei] = (pow(ustar_.value(),4)/sigmaEps_.value())*(1/pow(coord2t,2));
            //peterm2_[facei] = C1_.value()*pGref_[facei]*peref_[facei]/pkref_[facei] - C2_.value()*peref_[facei]*peref_[facei]/pkref_[facei];
            
            pkterm1_[facei] = (ustar_.value()*kappa_.value()/(z0_.value()*sigmak_.value()))*(4*B_.value()*coordt + C_.value());
            
            pkterm2_[facei] = pGref_[facei] - peref_[facei];
        }
    }
    
}
   
template<class BasicMomentumTransportModel>
void kEpsilonABL<BasicMomentumTransportModel>::errorsAndCmuComputation()
{
    bool debug=0;
    if (debug)
        Info << "Entered to computeErrorsAndCmu" << endl;
    Uerr_ = uError();
    Kerr_ = kError();
    Eerr_ = eError();
    Hybriderr_ = hybridError();
    BlendingFun_ = blendingFunction();
    NLCmu_ = nonlinearCmu();
    homCmu_ = homogeneousCmu();
    blendCmu_ = blendedCmu();

}

template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::uError()
{
    bool debug=0;
    if (debug)
        Info << "Entered to uError" << endl;
        
    volScalarField Udiff = mag((uref_ - this->U_.component(0))/max(uref_,usmall))/constantU_;
        
    forAll(Udiff,celli)
    {

        if(Udiff[celli] < thresholdU_.value())
        {
            Udiff[celli] = scalar(0);
        }
    }

    return min(Udiff,scalar(1));
}
 
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::kError()
{
    bool debug=0;
    if (debug)
        Info << "Entered to kError" << endl;
    volScalarField Kdiff = mag((kref_ - k_)/max(kref_,ksmall))/constantK_;
        
    forAll(Kdiff,celli)
    {

        if(Kdiff[celli] < thresholdK_.value())
        {
            Kdiff[celli] = scalar(0);
        }
    }
        
    return min(Kdiff,scalar(1));
}
    
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::eError()
{
    bool debug=0;
    if (debug)
        Info << "Entered to eError" << endl;

    volScalarField Ediff = mag((eref_ - epsilon_)/max(eref_,esmall))/constantE_;
        
    forAll(Ediff,celli)
    {

        if(Ediff[celli] < thresholdE_.value())
        {
            Ediff[celli] = scalar(0);
        }
    }

    return min(Ediff,scalar(1));
}
    
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::hybridError()
{
    bool debug=0;
    if (debug)
        Info << "Entered to hybridError" << endl;

    return max(max(Uerr_,Kerr_),Eerr_);
}
    
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::blendingFunction()
{
    bool debug=0;
    if (debug)
        Info << "Entered to blendingFunction" << endl;

    volScalarField marker = msel_.component(0)*Uerr_ + msel_.component(1)*Kerr_ + msel_.component(2)*Hybriderr_;

    volScalarField transErr = Foam::constant::mathematical::pi*max(marker- 0.5,-0.5);

    // Calculate and return blender value
    return pow(1.0 - 0.5*(1.0 + sin(transErr)),blendCoeff_);
}
    
/* Non-linear Eddy-viscosity model Cmu */
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::nonlinearCmu()
{
    bool debug=0;
    if (debug)
        Info << "Entered to nonlinearCmu" << endl;
    volSymmTensorField S = symm(fvc::grad(this->U_));
    volTensorField W = skew(fvc::grad(this->U_));

    volScalarField St = (this->k_/this->epsilon_)*sqrt(2.0)*mag(S);   //original
    volScalarField Ot = (this->k_/this->epsilon_)*sqrt(2.0)*mag(W);   //original

    return min((1/(0.9*pow(St,1.4)+0.4*pow(Ot,1.4)+3.5)),0.15);
}
    
// Height dependent blended Cmu
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::blendedCmu()
{
    bool debug=0;
    if (debug)
        Info << "Entered to blendedCmu" << endl;
            

    return NLCmu_ + (homCmu_ - NLCmu_)*blendingFunction();
}
    
// Smooth Cmu
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonABL<BasicMomentumTransportModel>::homogeneousCmu()
{
    bool debug=0;
    if (debug)
        Info << "Entered to homogeneousCmu" << endl;
    
    return min((pow(ustar_,4)/sqr(kref_)),CmuMax_);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsilonABL<BasicMomentumTransportModel>::kEpsilonABL
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    
    CmuMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuMax",
            this->coeffDict_,
            0.15
        )
    ),

    z0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "z0",
            this->coeffDict_,
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            0.01
        )
    ),
    
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.41
        )
    ),
    
    ustar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ustar",
            this->coeffDict_,
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            0.641
        )
    ),
    
    zDir_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "zDir",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            vector(0, 0, 1)
        )
    ),

    flowDir_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "flowDir",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            vector(1, 0, 0)
        )
    ),
    
    centres_
    (
        IOobject
        (
             "centres",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(zDir_ & this->mesh_.C())

    ),
    
    epsSource_
    (
        IOobject
        (
            "epsSource",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsSource", dimensionSet(0,2,-4,0,0,0,0),0.0)
    ),
    
    kSource_
    (
        IOobject
        (
            "kSource",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kSource", dimensionSet(0,2,-3,0,0,0,0),0.0)
    ),
    
    uref_
    (
        IOobject
        (
            "uref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("uref", dimensionSet(0,1,-1,0,0,0,0),0.0)
    ),

    kref_
    (
        IOobject
        (
            "kref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kref", dimensionSet(0,2,-2,0,0,0,0),0.0)
    ),

    eref_
    (
        IOobject
        (
            "eref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("eref", dimensionSet(0,2,-3,0,0,0,0),0.0)
    ),

    nref_
    (
        IOobject
        (
            "nref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("nref", dimensionSet(0,2,-1,0,0,0,0),0.0)
    ),
    
    Uerr_
    (
        IOobject
        (
            "Uerr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Uerr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    Kerr_
    (
        IOobject
        (
            "Kerr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Kerr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    Eerr_
    (
        IOobject
        (
            "Eerr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Eerr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    Hybriderr_
    (
        IOobject
        (
            "Hybriderr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Hybriderr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    BlendingFun_
    (
        IOobject
        (
            "BlendingFun",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("BlendingFun", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    A_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    B_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    C_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    D_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "D",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            1.3696
        )
    ),
    
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "E",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    blendCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "blendCoeff",
            this->coeffDict_,
            0
        )
    ),
    
    thresholdU_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "thresholdU",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.05
        )
    ),
    
    thresholdK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "thresholdK",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.3
        )
    ),
    
    thresholdE_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "thresholdE",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.1
        )
    ),
    
    constantU_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "constantU",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            1
        )
    ),
    
    constantK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "constantK",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            3
        )
    ),
    
    constantE_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "constantE",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            10
        )
    ),

    msel_
     (
         dimensioned<vector>::lookupOrAddToDict
         (
             "msel",
             this->coeffDict_,
             dimensionSet(0, 0, 0, 0, 0, 0, 0),
             vector(1, 0, 0)
         )
     ),

    usmall("usmall",dimensionSet(0,1,-1,0,0,0,0),scalar(1e-10)),
    ksmall("ksmall",dimensionSet(0,2,-2,0,0,0,0),scalar(1e-10)),
    esmall("esmall",dimensionSet(0,2,-3,0,0,0,0),scalar(1e-10)),
    
    NLCmu_
    (
        IOobject
        (
            "NLCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("NLCmu", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    homCmu_
    (
        IOobject
        (
            "homCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("homCmu", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    blendCmu_
    (
        IOobject
        (
            "blendCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("blendCmu", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    wSwitch_
    (
        IOobject
        (
            "wSwitch",
         this->runTime_.timeName(),
         this->mesh_,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
     this->mesh_,
        dimensionedScalar("wSwitch", dimensionSet(0,0,0,0,0,0,0),1.0)
    ),
    
    Gref_
    (
        IOobject
        (
            "Gref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Gref", dimensionSet(0,2,-3,0,0,0,0),0.0)
    ),
    
    eterm1_
    (
        IOobject
        (
            "eterm1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("eterm1", dimensionSet(0,2,-4,0,0,0,0),0.0)
    ),
    
    eterm2_
    (
        IOobject
        (
            "eterm2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("eterm2", dimensionSet(0,2,-4,0,0,0,0),0.0)
    ),
    
    kterm1_
    (
        IOobject
        (
            "kterm1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kterm1", dimensionSet(0,2,-3,0,0,0,0),0.0)
    ),
    
    kterm2_
    (
        IOobject
        (
            "kterm2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kterm2", dimensionSet(0,2,-3,0,0,0,0),0.0)
    ),
    
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    referencesComputation();
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    errorsAndCmuComputation();
    correctNut();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    
template<class BasicMomentumTransportModel>
bool kEpsilonABL<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        CmuMax_.readIfPresent(this->coeffDict());
        
        z0_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        ustar_.readIfPresent(this->coeffDict());
        zDir_.readIfPresent(this->coeffDict());
        flowDir_.readIfPresent(this->coeffDict());
        
        A_.readIfPresent(this->coeffDict());
        B_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());
        D_.readIfPresent(this->coeffDict());
        E_.readIfPresent(this->coeffDict());
        blendCoeff_.readIfPresent(this->coeffDict());
        
        thresholdU_.readIfPresent(this->coeffDict());
        thresholdK_.readIfPresent(this->coeffDict());
        thresholdE_.readIfPresent(this->coeffDict());
        constantU_.readIfPresent(this->coeffDict());
        constantK_.readIfPresent(this->coeffDict());
        constantE_.readIfPresent(this->coeffDict());
        
        msel_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kEpsilonABL<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    // const alphaField& alpha = this->alpha_;
    // const rhoField& rho = this->rho_;
    // const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    // fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();
    
    if (this->mesh_.changing())
    {
        referencesComputation();
    }

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    
    tgradU.clear();

    // Source terms computation
    epsSource_ = ABLepsilonSource()*blendingFunction();
    
    if(B_.value() + C_.value() != scalar(0))
        kSource_ = ABLkSource()*blendingFunction();
    else
        kSource_ = 0.0*ABLkSource();
    
    
    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
     (
         fvm::ddt(epsilon_)
       + fvm::div(phi, epsilon_)
       //- fvm::Sp(fvc::div(phi),epsilon_)
       - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
         C1_*G*epsilon_()/k_()
       - fvm::Sp(C2_*epsilon_()/k_(), epsilon_)
       + epsSource_
     );
    

    epsEqn.ref().relax();
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    bound(epsilon_, this->epsilonMin_);
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi, k_)
     - fvm::Sp(fvc::div(phi),k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_()/k_(), k_)
      + kSource_
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

    errorsAndCmuComputation();
    correctNut();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
