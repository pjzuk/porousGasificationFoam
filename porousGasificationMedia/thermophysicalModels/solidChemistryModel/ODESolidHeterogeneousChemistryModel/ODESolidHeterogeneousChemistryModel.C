/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "ODESolidHeterogeneousChemistryModel.H"
#include "solidChemistrySolver.H"
#include "reactingSolidHeterogeneousMixture.H"
#include "psiChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::
ODESolidHeterogeneousChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& compTypeName,
    const word& solidThermoName,
    PtrList<volScalarField>& gasPhaseGases
)
:
    CompType(mesh, obj, compTypeName, solidThermoName, gasPhaseGases),
    ODE(),
    gasPhaseGases_(gasPhaseGases),
    Ys_(this->solidThermo().composition().Y()),
    pyrolisisGases_
    (
        mesh.lookupObject<dictionary>
            ("chemistryProperties").lookup("species")
    ),
    reactions_
    (
        static_cast<const reactingSolidHeterogeneousMixture<SolidThermo>& >
            (this->solidThermo().composition())
    ),
    solidThermo_
    (
        static_cast<const reactingSolidHeterogeneousMixture<SolidThermo>& >
            (this->solidThermo().composition()).solidData()
    ),
    gasThermo_(
	dynamic_cast<const reactingMixture<GasThermo>&>
	(mesh.lookupObject <psiChemistryModel> ("chemistryProperties").thermo()).speciesData()
    ),
    nGases_(pyrolisisGases_.size()),
    nSpecie_(Ys_.size() + nGases_),
    nSolids_(Ys_.size()),
    nReaction_(reactions_.size()),
    solver_
    (
        solidChemistrySolver<CompType, SolidThermo, GasThermo>::New
        (
            *this,
            compTypeName,
	    solidThermoName
        )
    ),
    RRs_(nSolids_),
    RRg_(nGases_),
    shReactionHeat_
    (
	IOobject
        (
                "shReactionHeat",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    coeffs_(nSpecie_ + 2),
    Ys0_(nSolids_),
    cellCounter_(0),
    reactingCells_(mesh.nCells(), true),
    V_
    (
        IOobject
        (
            "V_",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("zero",dimVolume,0.0)
    ),
    solidReactionEnergyFromEnthalpy_
    (
        mesh.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("solidReactionEnergyFromEnthalpy",true)
    ),
    stoichiometricReactions_
    (
        mesh.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("stoichiometricReactions",false)
    ),
    diffusionLimitedReactions_
    (
        mesh.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("diffusionLimitedReactions",false)
    ),
    solidReactionDeltaEnergy_(0.0),
    showRRR_
    (
        mesh.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("showRelativeReactionRates",false)
    ),
    rhoG_
    (
        mesh.lookupObject<volScalarField>("rho")
    ),
    porosityF_
    (
        mesh.lookupObject<volScalarField>("porosityF")
    ),
    ST_
    (
        mesh.lookupObject<volScalarField>("STvol")
    )
{
    // create the fields for the chemistry sources

    Info << "gases in gas phase: " << gasPhaseGases_.size() << " \n" << endl;
    
    forAll(gasPhaseGases_,i)
    {
    	Info << gasPhaseGases_[i].name() << " \n";
    }
    
    Info << "gases from pyrolysis: " << endl;
    
    forAll(pyrolisisGases_,i)
    {
        Info << pyrolisisGases_[i] << " \n";
    }

    forAll(RRs_, fieldI)
    {
        RRs_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );


        IOobject header
        (
            Ys_[fieldI].name() + "0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            Ys0_.set
            (
                fieldI,
                new volScalarField
                (
                    IOobject
                    (
                        Ys_[fieldI].name() + "0",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField YsDefault
            (
                IOobject
                (
                    "YsDefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Ys0_.set
            (
                fieldI,
                new volScalarField
                (
                    IOobject
                    (
                        Ys_[fieldI].name() + "0",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    YsDefault
                )
            );

            // Calculate inital values of Ysi0 = rho*delta*Yi
            Ys0_[fieldI].internalField() =
                this->solidThermo().rho()
               *Ys_[fieldI]*mesh.V();
        }
    }

    V_.internalField() = mesh.V();

    forAll(RRg_, fieldI)
    {
        RRg_.set(fieldI, new scalarField(mesh.nCells(), 0.0));
    }

    Info<< "ODESolidHeterogeneousChemistryModel: Number of solids = " << nSolids_
        << " and reactions = " << nReaction_ << endl;

    Info<< "Number of gases from pyrolysis = " << nGases_ << endl;

    forAll(reactions_, i)
    {
        Info<< indent << "Reaction " << i << nl << reactions_[i] << nl;
    }

    Info << "solidReactionEnergyFromEnthalpy " << solidReactionEnergyFromEnthalpy_ << nl;
    Info << "stoichiometricReactions " << stoichiometricReactions_ << nl;

    gasDictionary_.resize(nGases_);
    forAll(gasDictionary_,gasI)
    {
        forAll(gasPhaseGases_,gasJ)
        {
            if (gasPhaseGases_[gasJ].name() == pyrolisisGases_[gasI])
            {
                gasDictionary_[gasI] = gasJ;
            }
        }
    }

    Info << "diffusionLimitedReactions " << diffusionLimitedReactions_ << nl;

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::
~ODESolidHeterogeneousChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
Foam::scalarField Foam::
ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& rR,
    const bool updateC0
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const label cellI = cellCounter_;
    label nEq = nEqns();

    if (not solidReactionEnergyFromEnthalpy_)
    {
	nEq = nEqns()+1;
    }

    scalarField om(nEq,0.0);

    forAll(reactions_, i)
    {
        const solidHeterogeneousReaction& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, 0.0, pf, cf, lRef, pr, cr, rRef
        );

 	scalar massCoefficient = 1;

	scalar substrates = 0;
	scalar products = 0;

	scalar solidSubstrates = 0;
	scalar solidProducts = 0;

        scalar massStream = 0;

        if (stoichiometricReactions_)
        {
            forAll(R.grhs(), g)
            {
                label gi = gasDictionary_[R.grhs()[g]];
                products += gasThermo_[gi].W()*R.grhsSto()[g];
            }
    
            forAll(R.glhs(), g)
            {
                label gi = R.glhs()[g];
                substrates += gasThermo_[gi].W()*R.glhsSto()[g];
            }
    
            forAll(R.slhs(), s)
            {
                solidSubstrates += R.slhsSto()[s];
            }
    
            forAll(R.srhs(), s)
            {
       	        solidProducts += R.srhsSto()[s];
            }

            if ((substrates + solidSubstrates  == 0) || (products + solidProducts == 0))
            {
                FatalErrorIn("omega")
                    << "Reaction:\n" << R
                    << "\nis not really a reaction" << exit(FatalError);
            }
   
            if (solidSubstrates > solidProducts)
            {
                scalar sr = solidProducts/solidSubstrates;
                massCoefficient = 1./(products-substrates); 

                forAll(R.slhs(), s)
                {
                    label si = R.slhs()[s];
                    om[si] -= omegai*R.slhsSto()[s]/solidSubstrates;
                    massStream -= omegai*R.slhsSto()[s]/solidSubstrates;
                    if (updateC0)
                    {
        		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                    } 
                }
        
                forAll(R.srhs(), s)
                {
                    label si = R.srhs()[s];
                    om[si] += sr*omegai*R.srhsSto()[s]/solidProducts;
                    if (updateC0)
                    {
        		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                    }
                }
       
                forAll(R.grhs(), g)
                {
                    om[R.grhs()[g] + nSolids_] +=  (1.0 - sr)*omegai*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g];
                }
        
                forAll(R.glhs(), g)
                {
                    label gi = R.glhs()[g];
                    om[gi + nSolids_] -=  (1.0 - sr)*massCoefficient*omegai*gasThermo_[gi].W()*R.glhsSto()[g];
                    massStream -= (1.0 - sr)*massCoefficient*omegai*gasThermo_[gi].W()*R.glhsSto()[g];
                }
 
                if (not solidReactionEnergyFromEnthalpy_)
                {
                    om[nEqns()] -= omegai*R.heatReact();
                }
            }
            else if (solidSubstrates < solidProducts)
            {
                if (solidSubstrates > 0) 
                {
                    scalar sr = solidProducts/solidSubstrates;
                    massCoefficient = 1./(products-substrates); 
        
                    forAll(R.slhs(), s)
                    {
                        label si = R.slhs()[s];
                        om[si] -= omegai*R.slhsSto()[s]/solidSubstrates;
                        massStream -= omegai*R.slhsSto()[s]/solidSubstrates;
                        if (updateC0)
                        {
            		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                        }
                    }
            
                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        om[si] += omegai*R.srhsSto()[s]/solidSubstrates;
            
                        if (updateC0)
                        {
            		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                        }
                    }
                   
                    forAll(R.grhs(), g)
                    {
                        om[R.grhs()[g] + nSolids_] +=  (1.0 - sr)*omegai*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g];
                    }
                
                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        om[gi + nSolids_] -=  (1.0 - sr)*omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g];
                        massStream -= (1.0 - sr)*omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g];
                    }

                    if (not solidReactionEnergyFromEnthalpy_)
                    {
                        om[nEqns()] -= omegai*R.heatReact();
                    }
                }
                else
                {
                    scalar sr = products/substrates;
                    forAll(R.grhs(), g)
                    {
                        om[R.grhs()[g] + nSolids_] +=  sr*omegai*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g]/products;
                    }
             
                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        om[gi + nSolids_] -=  omegai*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                        massStream -= omegai*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                    }
         
                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        om[si] += omegai*(1-sr)*R.srhsSto()[s]/solidProducts;
            
                        if (updateC0)
                        {
            		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                        }
                    }
                    if (not solidReactionEnergyFromEnthalpy_)
                    {
                        om[nEqns()] -= omegai*R.heatReact();
                    }
                }
            }
            else
            {
                if (products == substrates)
                {

                    forAll(R.slhs(), s)
                    {
                        label si = R.slhs()[s];
                        om[si] -= omegai*R.slhsSto()[s]/solidSubstrates;
                        massStream -= omegai*R.slhsSto()[s]/solidSubstrates;
                        if (updateC0)
                        {
            		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                        }
                    }
            
                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        om[si] += omegai*R.srhsSto()[s]/solidProducts;
                        if (updateC0)
                        {
            		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                        }
                    }
                    if (products > 0)
                    { 
                        forAll(R.grhs(), g)
                        {
                            om[R.grhs()[g] + nSolids_] +=  omegai*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g]/products;
                        }
                
                        forAll(R.glhs(), g)
                        {
                            label gi = R.glhs()[g];
                            om[gi + nSolids_] -=  omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                            massStream -= omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                        }
                    }     

                    if (not solidReactionEnergyFromEnthalpy_)
                    {
                        om[nEqns()] -= omegai*R.heatReact();
                    }

                }
                else
                {
                    FatalErrorIn("omega")
                        << "Reaction:\n" << R
                        << "\ntype is not implemented" << exit(FatalError);
                }
            } 
        }
        else
        {
            forAll(R.grhs(), g)
            {
                products += R.grhsSto()[g];
            }
    
            forAll(R.glhs(), g)
            {
                substrates += R.glhsSto()[g];
            }
    
            forAll(R.slhs(), s)
            {
                solidSubstrates += R.slhsSto()[s];
            }
    
            forAll(R.srhs(), s)
            {
       	        solidProducts += R.srhsSto()[s];
            }

            if ((substrates + solidSubstrates  == 0) || (products + solidProducts == 0))
            {
                FatalErrorIn("omega")
                    << "Reaction:\n" << R
                    << "\nis not really a reaction" << exit(FatalError);
            }

            scalar totalSubstrates = substrates + solidSubstrates;
  
            if ((totalSubstrates > 0) and (mag(totalSubstrates - (products + solidProducts)) < SMALL))
            {
                forAll(R.slhs(), s)
                {
                    label si = R.slhs()[s];
                    om[si] -= omegai*R.slhsSto()[s]/totalSubstrates;
                    massStream -= omegai*R.slhsSto()[s]/totalSubstrates;
                    if (updateC0)
                    {
        		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                    }
                }
        
                forAll(R.srhs(), s)
                {
                    label si = R.srhs()[s];
                    om[si] += omegai*R.srhsSto()[s]/totalSubstrates;
                    if (updateC0)
                    {
        		        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI]; 
                    }
                }
             
                forAll(R.grhs(), g)
                {
                    om[R.grhs()[g] + nSolids_] +=  omegai*R.grhsSto()[g]/totalSubstrates;
                }
        
                forAll(R.glhs(), g)
                {
                    label gi = R.glhs()[g];
                    om[gi + nSolids_] -= omegai*R.glhsSto()[g]/totalSubstrates;
                    massStream -= omegai*R.glhsSto()[g]/totalSubstrates;
                }
                if (not solidReactionEnergyFromEnthalpy_)
                {
                    om[nEqns()] -= omegai*R.heatReact();
                }
            }
            else
            {
                    FatalErrorIn("omega")
                        << "Reaction:\n" << R
                        << "\nviolates mass conservation" << exit(FatalError);
            }
        }
        rR[i] = mag(massStream);
    }
    return om;
}

template<class CompType, class SolidThermo, class GasThermo>
Foam::scalar
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::omega
(
    const solidHeterogeneousReaction& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
// eqZx2uHGn003
{
    scalarField c1(nSpecie_, 0.0);

    label cellI = cellCounter_;

    for (label i=0; i<nSpecie_; i++)
    {
        c1[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(T, 0.0, c1);

    const label Nl = R.slhs().size();
    if (Nl > 0)
    {
        if ( R.glhs().size() > 0 )
        {
            for (label s=0; s < Nl; s++)
            {
        	label si = R.slhs()[s];
                kf *= pow(Ys_[si][cellI],R.nReact()[s]);
            }
            forAll(R.glhs(),i)
            {
                kf *= pow(gasPhaseGases_[R.glhs()[i]].internalField()[cellI],R.nReact()[Nl+i]);
            }
        }
        else
        {
            for (label s=0; s<Nl; s++)
            {
                label si = R.slhs()[s];
                kf *= pow(Ys_[si][cellI],R.nReact()[s]);
            }
        }
        kf *= this->solidThermo().rho()[cellI];
    }
    else
    {
        forAll(R.glhs(),i)
        { 
            kf *= pow(gasPhaseGases_[R.glhs()[i]].internalField()[cellI],R.nReact()[i]);
        }
        kf *= rhoG_[cellI];
    }

    if (diffusionLimitedReactions_ and (R.glhs().size() > 0 ) and (kf != 0))
    {
        scalar avKf = 1./kf;
        forAll(R.glhs(),i)
        {
            scalar addAvKf = (ST_[cellI]*gasThermo_[R.glhs()[i]].alpha(T)*pow(gasPhaseGases_[R.glhs()[i]].internalField()[cellI],R.nReact()[Nl+i]));
            if (addAvKf != 0)
            {
                avKf += 1./addAvKf;
            }
            else
            {
                kf = 0;
            }
        }
        if (kf != 0)
        {
            kf = 1./avKf;
        }
    }

    return kf;
}


template<class CompType, class SolidThermo, class GasThermo>
void Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{

    scalar T = c[nSpecie_];

    dcdt = 0.0;
    scalarField rR(nReaction_,0.0);

    label celli = cellCounter_;

    scalarField omegaPreq(omega(c,T,0,rR)*(1.-porosityF_[celli]));
    if (solidReactionEnergyFromEnthalpy_)
    {
	    dcdt = omegaPreq;
    }
    else
    {
	    forAll(dcdt,i)
	    {
	        dcdt[i] = omegaPreq[i];
	    }
    }

    //Total mass concentration
    scalar cTot = 0.0;
    for (label i=0; i<nSolids_; i++)
    {
        cTot += c[i];
    }

    scalar newCp = 0.0;
    scalar newhi = 0.0;

    if (solidReactionEnergyFromEnthalpy_)
    {
	    for (label i=0; i<nSolids_; i++)
	    {
                scalar dYidt = dcdt[i];
                scalar Yi = c[i];
                newCp += Yi*solidThermo_[i].Cp(T);
                newhi -= dYidt*solidThermo_[i].hf();
	    }
    }
    else
    {
	    for (label i=0; i<nSolids_; i++)
	    {
                scalar Yi = c[i];
                newCp += Yi*solidThermo_[i].Cp(T);
	    }
        newhi += omegaPreq[nEqns()];
    }


    scalar dTdt = newhi/newCp;
    scalar dtMag = min(500.0, mag(dTdt));
    dcdt[nSpecie_] = dTdt*dtMag/(mag(dTdt) + 1.0e-10);

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class CompType, class SolidThermo, class GasThermo>
void Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{

    label cellI = cellCounter_;

    scalar T = c[nSpecie_];

    scalarField c2(nSpecie_, 0.0);
    scalarField rR(nReaction_, 0.0);

    for (label i=0; i<nSolids_; i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    scalarField omegaPreq(omega(c2,T,0.0,rR)*(1.-porosityF_[cellI]));
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdt = omegaPreq;
    }
    else
    {
        forAll(dcdt,i)
        {
            dcdt[i] = omegaPreq[i];
        }
    }

    for (label ri=0; ri<reactions_.size(); ri++)
    {
        const solidHeterogeneousReaction& R = reactions_[ri];

        scalar kf0 = R.kf(T, 0.0, c2);

        const label Ns = R.slhs().size();
        const label Ng = R.glhs().size();

        for (label rSj=0; rSj < Ns + Ng; rSj++)
        {
            label sj;
            if (rSj < Ns)
            {
                sj = R.slhs()[rSj];
            }
            else
            {
                sj = R.glhs()[rSj-Ns] +  nSolids_;
            }
            scalar kf = kf0;

            for (label rSi=0; rSi < Ns + Ng; rSi++)
            {
                label si;
                if (rSi < Ns)
                {
                    si = R.slhs()[rSi];
                }
                else
                {
                    si = R.glhs()[rSi-Ns] +  nSolids_;
                }

                scalar el = R.nReact()[rSi];
                if (rSi == rSj)
                {
                    if (el < 1.0)
                    {
                        if (c2[si]>SMALL)
                        {
                            kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c2[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c2[si], el);
                }
            }

            for (label rSi=0; rSi < Ns + Ng; rSi++)
            {
                label si;
                if (rSi < Ns)
                {
                    si = R.slhs()[rSi];
                }
                else
                {
                    si = R.glhs()[rSi-Ns] +  nSolids_;
                }
                dfdc[si][sj] -= kf;
            }

            forAll(R.srhs(), i)
            {
                label si = R.srhs()[i];
                dfdc[si][sj] += kf;
            }
            forAll(R.grhs(), i)
            {
                label gi = R.grhs()[i];
                dfdc[gi+nSolids_][sj] += kf;
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;

    scalarField dcdT0(dcdt);
    omegaPreq = omega(c2,T - delta ,0,rR)*(1.-porosityF_[cellI]);
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdT0 = omegaPreq;
    }
    else
    {
        forAll(dcdT0,i)
        {
            dcdT0[i] = omegaPreq[i];
        }
    }


    scalarField dcdT1(dcdt);
    omegaPreq = omega(c2,T + delta,0,rR)*(1.-porosityF_[cellI]);
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdT1 = omegaPreq;
    }
    else
    {
        forAll(dcdT1,i)
        {
            dcdT1[i] = omegaPreq[i];
        }
    }

    for (label i=0; i<nEqns(); i++)
    {
        dfdc[i][nSpecie_] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }

}


template<class CompType, class SolidThermo, class GasThermo>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::tc() const
{
    notImplemented
    (
        "ODESolidHeterogeneousChemistryModel::tc()"
    );

    return volScalarField::null();
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::Sh() const
// eqZx2uHGn004
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

	    if (solidReactionEnergyFromEnthalpy_) // eqZx2uHGn004
	    {
	        forAll(Ys_, i)
	        {
	    	    forAll(Sh, cellI)
	    	    {
	    	        scalar hf = solidThermo_[i].hf();
	    	        Sh[cellI] -= hf*RRs_[i][cellI];
	    	    }
	        }
	        forAll(Sh, cellI)
	        {
	    	    forAll(pyrolisisGases_, i)
	    	    {
	    	        scalar Hc = gasThermo_[gasDictionary_[i]].Hc();
	    	        Sh[cellI] -= Hc*RRg_[i][cellI];
	    	    }
	        }
	    }
	    else  // eqZx2uHGn017
	    {
	        forAll(Sh,cellI)
	        {
	    	    Sh[cellI] = shReactionHeat_[cellI];
	        }	
	    }
    }
    return tSh;
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::RRpor(const volScalarField T) const
{
    tmp<volScalarField> tRRpor
    (
        new volScalarField
        (
            IOobject
            (
                "RRpor",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& RRpor = tRRpor();

        forAll(Ys_, i)
        {
            forAll(RRpor, cellI)
            {
                scalar rho = solidThermo_[i].rho(T[cellI]);
                RRpor[cellI] -= RRs_[i][cellI]/rho;
            }
        }

    }
    return tRRpor;
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::label Foam::
ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::nEqns() const
{
    return (nSpecie_ + 2);
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::scalar
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::solve
(
    const scalar t0,
    const scalar deltaT
)
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    volScalarField newDeltaTMin
    (
        IOobject
        (
            "newDeltaTMin",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, GREAT)
    );


    if (this->mesh().changing())
    {
        forAll(RRs_, i)
        {
            RRs_[i].setSize(rho.size());
        }
        forAll(RRg_, i)
        {
            RRg_[i].setSize(rho.size());
        }
    }

    forAll(RRs_, i)
    {
        RRs_[i] = 0.0;
    }
    forAll(RRg_, i)
    {
        RRg_[i] = 0.0;
    }

    shReactionHeat_ *= 0.0;

    if (!this->chemistry_)
    {
        return GREAT;
    }
    else
    {
        scalar deltaTMin = GREAT;
        newDeltaTMin = GREAT;

        forAll(rho, celli)
        {
            if (reactingCells_[celli])
            {
                cellCounter_ = celli;

                scalar rhoi = rho[celli]*(1.-porosityF_[celli]);
                scalar rhoiG = rhoG_[celli]*porosityF_[celli];

                scalar Ti = this->solidThermo().T()[celli];

                scalarField c(nSpecie_, 0.0);
                scalarField c0(nSpecie_, 0.0);
                scalarField dc(nSpecie_, 0.0);
                scalarField rR(nReaction_, 0.0);

                scalarField omegaPreq(omega(c0,Ti,0.0,rR)*(1.-porosityF_[celli]));
                if (showRRR_ && gSum(rR) > 0)
                {
                    Info << "relative reaction rates: " << rR/gSum(rR) << endl;
                }

                for (label i=0; i<nSolids_; i++)
                {
                    c[i] = rhoi*Ys_[i][celli];
                }
                for (label i=0; i < nGases_; i++)
                {
                    c[nSolids_ + i] = rhoiG*gasPhaseGases_[i][celli];
                }
                c0 = c;

                scalar t = t0;
                scalar tauC = this->deltaTChem_[celli];
                scalar dt = min(deltaT, tauC);
                scalar timeLeft = deltaT;

                // calculate the source terms
                while (timeLeft > SMALL)
                {
                    tauC = solver().solve(c, Ti, 0.0, t, dt);
            	    t += dt;
                        // update the temperature
                        scalar cTot = 0.0;

                        //Total mass density
                        for (label i=0; i<nSolids_; i++)
                        {
                            cTot += c[i];
                        }

                        scalar newCp = 0.0;
                        scalar newhi = 0.0;
                        scalar invRho = 0.0;
                        scalarList dcdt = (c - c0)/dt;

            	    if (solidReactionEnergyFromEnthalpy_)
            	    {
                        for (label i=0; i<nSolids_; i++)
                        {
                    	scalar dYi = dcdt[i];
                    	scalar Yi = c[i];
                    	newhi -= dYi*solidThermo_[i].hf();
                    	newCp += Yi*solidThermo_[i].Cp(Ti);
                    	invRho += Yi/solidThermo_[i].rho(Ti);
                        }
            	    }
            	    else
            	    {
                        for (label i=0; i<nSolids_; i++)
                        {
                            scalar Yi = c[i];
                            newCp += Yi*solidThermo_[i].Cp(Ti);
                            invRho += Yi/solidThermo_[i].rho(Ti);
                        }
            	        newhi += omegaPreq[nEqns()];
            	    }

                    scalar dTi = (newhi/(newCp*rhoi))*dt;

                    Ti += dTi;
                    

                    timeLeft -= dt;
                    this->deltaTChem_[celli] = tauC;
                    dt = min(timeLeft, tauC);
                    dt = max(dt, SMALL);
                }
                deltaTMin = min(tauC, deltaTMin);

                dc = c - c0;

                // the integration inside time step tends to loose accuracy
                // for conserving the mass we correct for it
                // to normalize to calculated solid loss
                scalar sourceSolid = 0.0;
                scalar sourceGas = 0.0;
                scalar sourceCorrect = 1.0;
                for (label i=0; i<nSolids_; i++)
                {
                    sourceSolid += dc[i]/deltaT;
                }
                for (label i=0; i < nGases_; i++)
                {
                    sourceGas += dc[nSolids_ + i]/deltaT;
                }

                if (sourceGas != 0) 
                {
                    sourceCorrect = mag(sourceSolid/sourceGas);
                }

                scalar sumC = 0.;

                forAll(RRs_, i)
                {
                    RRs_[i][celli] = dc[i]/deltaT;
                    sumC += c[i];

            	    if ((rhoi*RRs_[i][celli] != 0) and (mag(RRs_[i][celli]) > ROOTVSMALL))
                    {
                        if (1. < c[i]/rhoi) 
                        {
                            newDeltaTMin[celli] = min(newDeltaTMin[celli], (1.01-Ys_[i][celli])/RRs_[i][celli]*rhoi);
                            deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                            if (1.02 < c[i]/rhoi)
                            {
                                Info << indent << indent << " too much   " << c[i]/rhoi << " of " << Ys_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << nl;
                            }
                        }
                        if (0. > c[i]/rhoi) 
                        {
                            newDeltaTMin[celli] = min(newDeltaTMin[celli], (-.01-Ys_[i][celli])/RRs_[i][celli]*rhoi);
                            deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                            if (-0.02 > c[i]/rhoi)
                            {
                                Info << indent << indent << " too little " << c[i]/rhoi << " of " << Ys_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << nl;
                            }
                        }
                        if (deltaTMin < 0) Info << dt << " " << deltaTMin << " " << tauC << " an error occured: negative deltaT from solid chemistry" << endl;
            	    } 
                }

                sumC=0;

                forAll(RRg_, i)
                {
                    RRg_[i][celli] = dc[nSolids_ + i]/deltaT*sourceCorrect;
                    sumC += c[nSolids_ + i];
                    if ((dc[nSolids_ + i]*rhoiG != 0) and (mag(RRg_[i][celli]) > ROOTVSMALL))
            	    {
                        scalar dtm = deltaTMin;
            	        if (1. < c[nSolids_ + i]/rhoiG)
                        {
                            newDeltaTMin[celli] = min(newDeltaTMin[celli], (1.01-gasPhaseGases_[i][celli])/RRg_[i][celli]*rhoiG);
                            deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                            if (1.02 < c[nSolids_ + i]/rhoiG)
                            {
                                Info << indent << indent << " too much   " << c[nSolids_ + i]/rhoiG << " of " << gasPhaseGases_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << " " << newDeltaTMin[celli] << nl;
                            }
                        }
            	        if (0. > c[nSolids_ + i]/rhoiG)
                        {
                            newDeltaTMin[celli] = min(newDeltaTMin[celli], (-.01-gasPhaseGases_[i][celli])/RRg_[i][celli]*rhoiG);
                            deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                            if (-0.02 > c[nSolids_ + i]/rhoiG)
                            {
                                Info << indent << indent << " too little " << c[nSolids_ + i]/rhoiG << " of " << gasPhaseGases_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << nl;
                            }
                        }
                        if (deltaTMin < 0) Info << dt  << " " << dtm << " " << deltaTMin << " " << tauC << " " << (1.-gasPhaseGases_[i][celli]) << " " << c[nSolids_ + i]/rhoiG << " " << gasPhaseGases_[i][celli]  << " an error occured: negative deltaT from solid chemistry" << endl;
                    }
                }

                if (not solidReactionEnergyFromEnthalpy_)
                { 
                    shReactionHeat_[celli] = omegaPreq[nEqns()];
                }
            
            }
        }
        // Don't allow the time-step to rise
        deltaTMin = min(deltaTMin, deltaT);
        reduce(deltaTMin,minOp<scalar>());
        scalar output = gMin(newDeltaTMin);
        return output;
    }
}


template<class CompType, class SolidThermo,class GasThermo>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::gasHs
(
    const volScalarField& T,
    const label index
) const
{

    tmp<volScalarField> tHs
    (
        new volScalarField
        (
            IOobject
            (
                "Hs_" + pyrolisisGases_[index],
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimMass, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& gasHs = tHs();

    const GasThermo& mixture = gasThermo_[index];

    forAll(gasHs.internalField(), cellI)
    {
        gasHs[cellI] = mixture.Hs(T[cellI]);
    }

    return tHs;
}


template<class CompType, class SolidThermo,class GasThermo>
Foam::scalar
Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    notImplemented
    (
        "ODESolidHeterogeneousChemistryModel::solve"
        "("
            "scalarField&, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar"
        ")"
    );
    return (0);
}


template<class CompType, class SolidThermo,class GasThermo>
void Foam::ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>::
setCellReacting(const label cellI, const bool active)
{
    reactingCells_[cellI] = active;
}


// ************************************************************************* //
