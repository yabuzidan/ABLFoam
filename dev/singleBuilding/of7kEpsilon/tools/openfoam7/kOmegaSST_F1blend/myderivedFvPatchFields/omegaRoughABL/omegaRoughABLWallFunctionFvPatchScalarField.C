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

\*---------------------------------------------------------------------------*/

#include "omegaRoughABLWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaRoughABLWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::omegaRoughABLWallFunctionFvPatchScalarField::blendingFunction
(
    const scalar BIA,
    const scalar blendCoeff
) const
{
    scalar transErr = Foam::constant::mathematical::pi*max(BIA-0.5,-0.5);
    /* Calculate and return blender value */
   return pow(1.0 - 0.5*(1.0 + sin(transErr)),blendCoeff);
}

void omegaRoughABLWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaRoughABLWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaRoughABLWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void omegaRoughABLWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaRoughABLWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                label celli = faceCells[i];
                weights[celli]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchi = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


omegaRoughABLWallFunctionFvPatchScalarField&
omegaRoughABLWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const omegaRoughABLWallFunctionFvPatchScalarField& opf =
        refCast<const omegaRoughABLWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaRoughABLWallFunctionFvPatchScalarField&>(opf);
}


void omegaRoughABLWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaRoughABLWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaRoughABLWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void omegaRoughABLWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> te = turbModel.epsilon();
    const volScalarField& epsilon = te();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const FieldType& G =
        db().lookupObject<FieldType>(turbModel.GName());

    /*const volVectorField& Ut = turbModel.U();
    vectorField Uwt = Uw;
    forAll(nutw, facei)
    {
      label faceCelli = patch.faceCells()[facei];
      Uwt[facei] = Ut[faceCelli];
    }*/
    // Get first cell centroid
    const scalarField ccentre = zDir_ & patch.Cn();
    /*forAll(nutw, facei)
    {
        if(roughWall_[facei] > 0)
        {
            ccentre = y[facei];
        }
    }*/

    const scalarField coord = (ccentre + z0_)/z0_;
    const scalarField coord2 = (ccentre + z0_);
    
    // Calculate references
    const scalarField uref = ustar_*log(coord)/kappaWF_;   
    const scalarField kref = A_*log(coord) + B_*sqr(coord) + C_*(coord)+ D_ +Etke_*log(coord2);
    const scalarField hombetaStar = min(pow(ustar_,4)/sqr(kref),betaStarMax_);
    const scalarField oref = (ustar_/sqrt(hombetaStar))/(kappaWF_*coord2);

    // Calculate errors
    const volVectorField& Ut = turbModel.U();
    vectorField Uwt = Uw;
    forAll(nutw, facei)
    {
      label faceCelli = patch.faceCells()[facei];
      Uwt[facei] = Ut[faceCelli];
    }
    
    scalarField Uerr  = mag((uref - Uwt.component(0))/max(uref,uref*1e-14))/constantU_;
    scalarField kerr = mag((kref - k)/max(kref,kref*1e-14))/constantK_;
    scalarField Oerr = mag((oref - (epsilon/(hombetaStar*k)))/max(oref,oref*1e-14))/constantO_;
    
    forAll(Uerr, facei)
  
    {
        if(Uerr[facei] < thresholdU_[facei])
        {
            Uerr[facei] = scalar(0.0);
        }
        else if(Uerr[facei] > scalar(1))
        {
            Uerr[facei] = scalar(1.0);;
        }
    }

    forAll(kerr, facei)
  
    {
        if(kerr[facei] < thresholdK_[facei])
        {
            kerr[facei] = scalar(0);
        }
        if(kerr[facei] > scalar(1))
        {
            kerr[facei] = scalar(1);
        }
    }
    
    forAll(Oerr, facei)
  
    {
        if(Oerr[facei] < thresholdO_[facei])
        {
            Oerr[facei] = scalar(0);
        }
        if(Oerr[facei] > scalar(1))
        {
            Oerr[facei] = scalar(1);
        }
    }

    scalarField herr = max(max(Uerr,kerr),Oerr);

    const labelUList& faceCells = patch.faceCells();


    // Set omega and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];


    	const scalar betaStar25 = pow025(hombetaStar[facei]);
    	const scalar betaStar5 = sqrt(hombetaStar[facei]);

        const scalar w = cornerWeights[facei];

        const scalar Rey = y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar yPlus = betaStar25*Rey;
        const scalar uPlus = (1/nutw.kappa())*log(nutw.E()*yPlus);

        if (blended_)
        {
            const scalar lamFrac = exp(-Rey/11);
            const scalar turbFrac = 1 - lamFrac;

            const scalar uStar = sqrt
            (
                lamFrac*nuw[facei]*magGradUw[facei] + turbFrac*betaStar5*k[celli]
            );

            const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));
            const scalar omegaLog = uStar/(betaStar5*nutw.kappa()*y[facei]);

            omega0[celli] += w*(lamFrac*omegaVis + turbFrac*omegaLog);

            G0[celli] +=
                w
               *(
                   lamFrac*G[celli]

                 + turbFrac
                  *sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                  /(nuw[facei]*nutw.kappa()*yPlus)
               );
        }
        else
        {
            /*if (yPlus < nutw.yPlusLam())
            {
                const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));

                omega0[celli] += w*omegaVis;

                G0[celli] += w*G[celli];
            }
            else
            {
                const scalar uStar = sqrt(betaStar5*k[celli]);
                const scalar omegaLog = uStar/(betaStar5*nutw.kappa()*y[facei]);*/

                //omega0[celli] += w*omegaLog;
		omega0[celli] += w*(1/(nutw.kappa()*ustar_[facei]*(y[facei] + z0_[facei])))
				  *k[celli];

                G0[celli] +=
            	w
           	*(nutw[facei] + nuw[facei])
           	*magGradUw[facei]
           	*betaStar25*sqrt(k[celli])
           /(nutw.kappa()*(y[facei] + z0_[facei]));
                    //w*
                    //sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                   ///(nuw[facei]*nutw.kappa()*yPlus);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaRoughABLWallFunctionFvPatchScalarField::omegaRoughABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    beta1_(0.075),
    blended_(false),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_(p.size(), 0.0),
    ustar_(p.size(), 0.0),
    kappaWF_(p.size(), 0.0),
    A_(p.size(), 0.0),
    B_(p.size(), 0.0),
    C_(p.size(), 0.0),
    D_(p.size(), 1.37),
    Etke_(p.size(), 0),
    WFbetaStar_(p.size(),0.15),
    betaStarMax_(p.size(),0.15),
    blendCoeff_(p.size(), 0.0),
    zDir_(p.size(), vector(0, 0, 1)),
    roughWall_(p.size(), 1.0),
    thresholdU_(p.size(), 0.3),
    thresholdO_(p.size(), 0.5),
    thresholdK_(p.size(), 0.5),
    constantU_(p.size(), 1),
    constantO_(p.size(), 3),
    constantK_(p.size(), 10)
{}


omegaRoughABLWallFunctionFvPatchScalarField::omegaRoughABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    blended_(dict.lookupOrDefault<Switch>("blended", false)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_("z0", dict, p.size()),
    ustar_("ustar", dict, p.size()),
    kappaWF_("kappaWF", dict, p.size()),
    A_("A", dict, p.size()),
    B_("B", dict, p.size()),
    C_("C", dict, p.size()),
    D_("D", dict, p.size()),
    Etke_("Etke", dict, p.size()),
    WFbetaStar_("WFbetaStar", dict, p.size()),
    betaStarMax_("betaStarMax", dict, p.size()),
    blendCoeff_("blendCoeff", dict, p.size()),
    zDir_("zDir", dict, p.size()),
    roughWall_("roughWall", dict, p.size()),
    thresholdU_("thresholdU", dict, p.size()),
    thresholdO_("thresholdO", dict, p.size()),
    thresholdK_("thresholdK", dict, p.size()),
    constantU_("constantU", dict, p.size()),
    constantO_("constantO", dict, p.size()),
    constantK_("constantK", dict, p.size())
{
    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


omegaRoughABLWallFunctionFvPatchScalarField::omegaRoughABLWallFunctionFvPatchScalarField
(
    const omegaRoughABLWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    beta1_(ptf.beta1_),
    blended_(ptf.blended_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_(mapper(ptf.z0_)),
    ustar_(mapper(ptf.ustar_)),
    kappaWF_(mapper(ptf.kappaWF_)),
    A_(mapper(ptf.A_)),
    B_(mapper(ptf.B_)),
    C_(mapper(ptf.C_)),
    D_(mapper(ptf.D_)),
    Etke_(mapper(ptf.Etke_)),
    WFbetaStar_(mapper(ptf.WFbetaStar_)),
    betaStarMax_(mapper(ptf.betaStarMax_)),
    blendCoeff_(mapper(ptf.blendCoeff_)),
    zDir_(mapper(ptf.zDir_)),
    roughWall_(mapper(ptf.roughWall_)),
    thresholdU_(mapper(ptf.thresholdU_)),
    thresholdO_(mapper(ptf.thresholdO_)),
    thresholdK_(mapper(ptf.thresholdK_)),
    constantU_(mapper(ptf.constantU_)),
    constantO_(mapper(ptf.constantO_)),
    constantK_(mapper(ptf.constantK_))
{}


omegaRoughABLWallFunctionFvPatchScalarField::omegaRoughABLWallFunctionFvPatchScalarField
(
    const omegaRoughABLWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_(owfpsf.z0_),
    ustar_(owfpsf.ustar_),
    kappaWF_(owfpsf.kappaWF_),
    A_(owfpsf.A_),
    B_(owfpsf.B_),
    C_(owfpsf.C_),
    D_(owfpsf.D_),
    Etke_(owfpsf.Etke_),
    WFbetaStar_(owfpsf.WFbetaStar_),
    betaStarMax_(owfpsf.betaStarMax_),
    blendCoeff_(owfpsf.blendCoeff_),
    zDir_(owfpsf.zDir_),
    roughWall_(owfpsf.roughWall_),
    thresholdU_(owfpsf.thresholdU_),
    thresholdO_(owfpsf.thresholdO_),
    thresholdK_(owfpsf.thresholdK_),
    constantU_(owfpsf.constantU_),
    constantO_(owfpsf.constantO_),
    constantK_(owfpsf.constantK_)

{}


omegaRoughABLWallFunctionFvPatchScalarField::omegaRoughABLWallFunctionFvPatchScalarField
(
    const omegaRoughABLWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_(owfpsf.z0_),
    ustar_(owfpsf.ustar_),
    kappaWF_(owfpsf.kappaWF_),
    A_(owfpsf.A_),
    B_(owfpsf.B_),
    C_(owfpsf.C_),
    D_(owfpsf.D_),
    WFbetaStar_(owfpsf.WFbetaStar_),
    betaStarMax_(owfpsf.betaStarMax_),
    blendCoeff_(owfpsf.blendCoeff_),
    zDir_(owfpsf.zDir_),
    roughWall_(owfpsf.roughWall_),
    thresholdU_(owfpsf.thresholdU_),
    thresholdO_(owfpsf.thresholdO_),
    thresholdK_(owfpsf.thresholdK_),
    constantU_(owfpsf.constantU_),
    constantO_(owfpsf.constantO_),
    constantK_(owfpsf.constantK_)
{}

void Foam::omegaRoughABLWallFunctionFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
    m(z0_,z0_);
    m(ustar_,ustar_);
    m(kappaWF_,kappaWF_);
    m(A_,A_);
    m(B_,B_);
    m(C_,C_);
    m(D_,D_);
    m(Etke_,Etke_);
    m(WFbetaStar_,WFbetaStar_);
    m(betaStarMax_,betaStarMax_);
    m(blendCoeff_,blendCoeff_);
    m(zDir_,zDir_);
    m(roughWall_,roughWall_);
    m(thresholdU_,thresholdU_);
    m(thresholdO_,thresholdO_);
    m(thresholdK_,thresholdK_);
    m(constantU_,constantU_);
    m(constantO_,constantO_);
    m(constantK_,constantK_);
}

void Foam::omegaRoughABLWallFunctionFvPatchScalarField::rmap
(
    const omegaRoughABLWallFunctionFvPatchScalarField& blptf,
    const labelList& addr
)
{
    z0_.rmap(blptf.z0_,addr);
    ustar_.rmap(blptf.ustar_,addr);
    kappaWF_.rmap(blptf.kappaWF_,addr);
    A_.rmap(blptf.A_, addr);
    B_.rmap(blptf.B_, addr);
    C_.rmap(blptf.C_, addr);
    D_.rmap(blptf.D_, addr);
    Etke_.rmap(blptf.Etke_, addr);
    WFbetaStar_.rmap(blptf.WFbetaStar_, addr);
    betaStarMax_.rmap(blptf.betaStarMax_, addr);
    blendCoeff_.rmap(blptf.blendCoeff_, addr);
    zDir_.rmap(blptf.zDir_,addr);
    roughWall_.rmap(blptf.roughWall_, addr);
    thresholdU_.rmap(blptf.thresholdU_, addr);
    thresholdO_.rmap(blptf.thresholdO_, addr);
    thresholdK_.rmap(blptf.thresholdK_, addr);
    constantU_.rmap(blptf.constantU_, addr);
    constantO_.rmap(blptf.constantO_, addr);
    constantK_.rmap(blptf.constantK_, addr);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& omegaRoughABLWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalarField& omegaRoughABLWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void omegaRoughABLWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        omega[celli] = omega0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaRoughABLWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaRoughABLWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaRoughABLWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintomega(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintomega.append(omega[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintomega)
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaRoughABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeEntry(os,"z0", z0_);
    writeEntry(os,"ustar", ustar_);
    writeEntry(os,"kappaWF", kappaWF_);
    writeEntry(os,"A", A_);
    writeEntry(os,"B", B_);
    writeEntry(os,"C", C_);
    writeEntry(os,"D", D_);
    writeEntry(os,"Etke", Etke_);
    writeEntry(os,"WFbetaStar", WFbetaStar_);
    writeEntry(os,"betaStarMax", betaStarMax_);
    writeEntry(os,"blendCoeff", blendCoeff_);
    writeEntry(os,"zDir", zDir_);
    writeEntry(os,"roughWall", roughWall_);
    writeEntry(os,"thresholdU", thresholdU_);
    writeEntry(os,"thresholdO", thresholdO_);
    writeEntry(os,"thresholdK", thresholdK_);
    writeEntry(os,"constantU", constantU_);
    writeEntry(os,"constantO", constantO_);
    writeEntry(os,"constantK", constantK_);   
    writeEntry(os, "beta1", beta1_);
    writeEntry(os, "blended", blended_);
    //fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaRoughABLWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
