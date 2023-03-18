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
    Foam::EpsilonABL

Description
    This boundary condition specifies an inlet value for the turbulence
    dissipation, \f$\epsilon\f$, appropriate for atmospheric boundary layers.

    See Foam::mainABL for details.

    Example of the boundary condition specification:
    \verbatim
    ground
    {
        type            mainABLInletEpsilon;
        z               (0 0 1);
        Uref            10.0;
        Zref            20.0;
        z0              uniform 0.1;
        zGround         uniform 0.0;
    }
    \endverbatim

See also
    Foam::mainABL,
    Foam::mainABLInletVelocityFvPatchVectorField,
    Foam::mainABLInletKFvPatchScalarField

SourceFiles
    EpsilonABL.C

\*---------------------------------------------------------------------------*/

#ifndef EpsilonABL_H
#define EpsilonABL_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "mainABL.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class EpsilonABL Declaration
\*---------------------------------------------------------------------------*/

class EpsilonABL
:
    public inletOutletFvPatchScalarField,
    public mainABL
{

public:

    //- Runtime type information
    TypeName("EpsilonABL");


    // Constructors

        //- Construct from patch and internal field
        EpsilonABL
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        EpsilonABL
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  EpsilonABL onto a new patch
        EpsilonABL
        (
            const EpsilonABL&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new EpsilonABL(*this)
            );
        }

        //- Construct as copy setting internal field reference
        EpsilonABL
        (
            const EpsilonABL&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new EpsilonABL(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
