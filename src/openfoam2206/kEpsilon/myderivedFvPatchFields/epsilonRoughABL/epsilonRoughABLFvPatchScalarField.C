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

#include "epsilonRoughABLFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonRoughABLFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
Foam::scalar Foam::epsilonRoughABLFvPatchScalarField::blendingFunction
(
    const scalar BIA,
    const scalar blendCoeff
) const
{
    scalar transErr = Foam::constant::mathematical::pi*max(BIA-0.5,-0.5);
    /* Calculate and return blender value */
   return pow(1.0 - 0.5*(1.0 + sin(transErr)),blendCoeff);
}

void Foam::epsilonRoughABLFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonRoughABLFvPatchScalarField>(bf[patchi]))
        {
            epsilonRoughABLFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonRoughABLFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

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

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonRoughABLFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchi = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    epsilon_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::epsilonRoughABLFvPatchScalarField&
Foam::epsilonRoughABLFvPatchScalarField::epsilonPatch(const label patchi)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const epsilonRoughABLFvPatchScalarField& epf =
        refCast<const epsilonRoughABLFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonRoughABLFvPatchScalarField&>(epf);
}


void Foam::epsilonRoughABLFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonRoughABLFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonRoughABLFvPatchScalarField& epf = epsilonPatch(patchi);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void Foam::epsilonRoughABLFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    
    const tmp<volScalarField> te = turbModel.epsilon();
    const volScalarField& epsilon = te();


    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Get first cell centroid
    scalarField ccentre = zDir_ & patch.Cn();
    forAll(nutw, facei)
    {
        if(roughWall_[facei] > 0)
        {
            ccentre = y[facei];
        }
    }

    const scalarField coord = (ccentre + z0_)/z0_;
    const scalarField coord2 = (ccentre + z0_);
    
    // Calculate references
    const scalarField uref = ustar_*log(coord)/kappaWF_;
    const scalarField kref = A_*log(coord) + B_*sqr(coord) + C_*(coord)+ D_ +Etke_*log(coord2);
    const scalarField homCmu = min(pow(ustar_,4)/sqr(kref),CmuMax_);
    const scalarField eref = (sqrt(homCmu)*ustar_*kref)/(kappaWF_*coord2);
    
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
    scalarField Eerr = mag((eref - epsilon)/max(eref,eref*1e-14))/constantE_;
    
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
    
    forAll(Eerr, facei)
  
    {
        if(Eerr[facei] < thresholdE_[facei])
        {
            Eerr[facei] = scalar(0);
        }
        if(Eerr[facei] > scalar(1))
        {
            Eerr[facei] = scalar(1);
        }
    }

    scalarField herr = max(max(Uerr,kerr),Eerr);
    
    // Set epsilon and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];
        
        scalar bCmu = WFCmu_[facei] + (homCmu[facei] - WFCmu_[facei])*blendingFunction(herr[facei],blendCoeff_[facei]);
        
        const scalar Cmu25 = pow025(bCmu);
        const scalar Cmu05 = pow(bCmu,0.5);
        const scalar Cmu75 = pow(bCmu,0.75);

        const scalar w = cornerWeights[facei];
        
        if(roughWall_[facei] > scalar(0)) // Rough case
        {
            /*epsilon0[celli] +=
            w*Cmu05*k[celli]*ustar_[facei]/(kappaWF_[facei]*(y[facei]+z0_[facei]));*/
            epsilon0[celli] +=
            w*Cmu75*pow(k[celli],1.5)/(kappaWF_[facei]*(y[facei]+z0_[facei]));

            // Production of k - Original implementation (faulty?)
            // G0[celli] +=
            //     w
            //    *(0.52)*log((2*y[facei]+z0_[facei])/z0_[facei])*(nutw[facei] + nuw[facei])
            //    *magGradUw[facei]
            //    *Cmu25*sqrt(k[celli])
            //   /(kappaWF_[facei]*(y[facei] + z0_[facei]));
        
            // Production of k (Parante 2011 Eq 24)
            G0[celli] +=
                w
                *pow((nutw[celli] + nuw[facei])
                *magGradUw[celli],2)
                /(Cmu25*sqrt(k[celli])
                *kappaWF_[facei]*(y[facei] + z0_[facei]));
        }
        else
        {
            epsilon0[celli] +=
            w*Cmu05*k[celli]*ustar_[facei]/(kappaWF_[facei]*(y[facei]));

            G0[celli] +=
                w
               *(nutw[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(k[celli])
              /(kappaWF_[facei]*(y[facei]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonRoughABLFvPatchScalarField::
epsilonRoughABLFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    G_(),
    epsilon_(),
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
    WFCmu_(p.size(),0.15),
    CmuMax_(p.size(),0.15),
    blendCoeff_(p.size(), 0.0),
    zDir_(p.size(), vector(0, 0, 1)),
    roughWall_(p.size(), 1.0),
    thresholdU_(p.size(), 0.3),
    thresholdE_(p.size(), 0.5),
    thresholdK_(p.size(), 0.5),
    constantU_(p.size(), 1),
    constantE_(p.size(), 3),
    constantK_(p.size(), 10)
{}


Foam::epsilonRoughABLFvPatchScalarField::
epsilonRoughABLFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    G_(),
    epsilon_(),
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
    WFCmu_("WFCmu", dict, p.size()),
    CmuMax_("CmuMax", dict, p.size()),
    blendCoeff_("blendCoeff", dict, p.size()),
    zDir_("zDir", dict, p.size()),
    roughWall_("roughWall", dict, p.size()),
    thresholdU_("thresholdU", dict, p.size()),
    thresholdE_("thresholdE", dict, p.size()),
    thresholdK_("thresholdK", dict, p.size()),
    constantU_("constantU", dict, p.size()),
    constantE_("constantE", dict, p.size()),
    constantK_("constantK", dict, p.size())
{
    // Apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::epsilonRoughABLFvPatchScalarField::
epsilonRoughABLFvPatchScalarField
(
    const epsilonRoughABLFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    G_(),
    epsilon_(),
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
    WFCmu_(mapper(ptf.WFCmu_)),
    CmuMax_(mapper(ptf.CmuMax_)),
    blendCoeff_(mapper(ptf.blendCoeff_)),
    zDir_(mapper(ptf.zDir_)),
    roughWall_(mapper(ptf.roughWall_)),
    thresholdU_(mapper(ptf.thresholdU_)),
    thresholdE_(mapper(ptf.thresholdE_)),
    thresholdK_(mapper(ptf.thresholdK_)),
    constantU_(mapper(ptf.constantU_)),
    constantE_(mapper(ptf.constantE_)),
    constantK_(mapper(ptf.constantK_))
{}


Foam::epsilonRoughABLFvPatchScalarField::
epsilonRoughABLFvPatchScalarField
(
    const epsilonRoughABLFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_(ewfpsf.z0_),
    ustar_(ewfpsf.ustar_),
    kappaWF_(ewfpsf.kappaWF_),
    A_(ewfpsf.A_),
    B_(ewfpsf.B_),
    C_(ewfpsf.C_),
    D_(ewfpsf.D_),
    Etke_(ewfpsf.Etke_),
    WFCmu_(ewfpsf.WFCmu_),
    CmuMax_(ewfpsf.CmuMax_),
    blendCoeff_(ewfpsf.blendCoeff_),
    zDir_(ewfpsf.zDir_),
    roughWall_(ewfpsf.roughWall_),
    thresholdU_(ewfpsf.thresholdU_),
    thresholdE_(ewfpsf.thresholdE_),
    thresholdK_(ewfpsf.thresholdK_),
    constantU_(ewfpsf.constantU_),
    constantE_(ewfpsf.constantE_),
    constantK_(ewfpsf.constantK_)

{}


Foam::epsilonRoughABLFvPatchScalarField::
epsilonRoughABLFvPatchScalarField
(
    const epsilonRoughABLFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    z0_(ewfpsf.z0_),
    ustar_(ewfpsf.ustar_),
    kappaWF_(ewfpsf.kappaWF_),
    A_(ewfpsf.A_),
    B_(ewfpsf.B_),
    C_(ewfpsf.C_),
    D_(ewfpsf.D_),
    WFCmu_(ewfpsf.WFCmu_),
    CmuMax_(ewfpsf.CmuMax_),
    blendCoeff_(ewfpsf.blendCoeff_),
    zDir_(ewfpsf.zDir_),
    roughWall_(ewfpsf.roughWall_),
    thresholdU_(ewfpsf.thresholdU_),
    thresholdE_(ewfpsf.thresholdE_),
    thresholdK_(ewfpsf.thresholdK_),
    constantU_(ewfpsf.constantU_),
    constantE_(ewfpsf.constantE_),
    constantK_(ewfpsf.constantK_)
{}

void Foam::epsilonRoughABLFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
    z0_.autoMap(m);
    ustar_.autoMap(m);
    kappaWF_.autoMap(m);
    A_.autoMap(m);
    B_.autoMap(m);
    C_.autoMap(m);
    D_.autoMap(m);
    Etke_.autoMap(m);
    WFCmu_.autoMap(m);
    CmuMax_.autoMap(m);
    blendCoeff_.autoMap(m);
    zDir_.autoMap(m);
    roughWall_.autoMap(m);
    thresholdU_.autoMap(m);
    thresholdE_.autoMap(m);
    thresholdK_.autoMap(m);
    constantU_.autoMap(m);
    constantE_.autoMap(m);
    constantK_.autoMap(m);
}

void Foam::epsilonRoughABLFvPatchScalarField::rmap
(
    const epsilonRoughABLFvPatchScalarField& blptf,
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
    WFCmu_.rmap(blptf.WFCmu_, addr);
    CmuMax_.rmap(blptf.CmuMax_, addr);
    blendCoeff_.rmap(blptf.blendCoeff_, addr);
    zDir_.rmap(blptf.zDir_,addr);
    roughWall_.rmap(blptf.roughWall_, addr);
    thresholdU_.rmap(blptf.thresholdU_, addr);
    thresholdE_.rmap(blptf.thresholdE_, addr);
    thresholdK_.rmap(blptf.thresholdK_, addr);
    constantU_.rmap(blptf.constantU_, addr);
    constantE_.rmap(blptf.constantE_, addr);
    constantK_.rmap(blptf.constantK_, addr);
}
void Foam::epsilonRoughABLFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    os.writeEntry("z0", z0_);
    os.writeEntry("ustar", ustar_);
    os.writeEntry("kappaWF", kappaWF_);
    os.writeEntry("A", A_);
    os.writeEntry("B", B_);
    os.writeEntry("C", C_);
    os.writeEntry("D", D_);
    os.writeEntry("Etke", Etke_);
    os.writeEntry("WFCmu", WFCmu_);
    os.writeEntry("CmuMax", CmuMax_);
    os.writeEntry("blendCoeff", blendCoeff_);
    os.writeEntry("zDir", zDir_);
    os.writeEntry("roughWall", roughWall_);
    os.writeEntry("thresholdU", thresholdU_);
    os.writeEntry("thresholdE", thresholdE_);
    os.writeEntry("thresholdK", thresholdK_);
    os.writeEntry("constantU", constantU_);
    os.writeEntry("constantE", constantE_);
    os.writeEntry("constantK", constantK_);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::epsilonRoughABLFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


Foam::scalarField& Foam::epsilonRoughABLFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::epsilonRoughABLFvPatchScalarField::updateCoeffs()
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
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        epsilon[celli] = epsilon0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonRoughABLFvPatchScalarField::updateWeightedCoeffs
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
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    scalarField& epsilonf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            epsilon[celli] = (1.0 - w)*epsilon[celli] + w*epsilon0[celli];
            epsilonf[facei] = epsilon[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonRoughABLFvPatchScalarField::manipulateMatrix
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


void Foam::epsilonRoughABLFvPatchScalarField::manipulateMatrix
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
    DynamicList<scalar> constraintEpsilon(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilon
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintEpsilon.append(epsilon[celli]);
        }
    }

    bool debug =1;
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
        scalarField(constraintEpsilon)
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonRoughABLFvPatchScalarField
    );
}


// ************************************************************************* //
