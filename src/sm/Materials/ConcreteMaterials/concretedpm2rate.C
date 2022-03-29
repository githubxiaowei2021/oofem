/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "concretedpm2rate.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "timestep.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"
#include <limits>



namespace oofem {
REGISTER_Material(ConcreteDPM2Rate);
  
  //****************************************
  //Status of ConcreteDPM2Rate
  //*************************************
  
ConcreteDPM2RateStatus::ConcreteDPM2RateStatus(GaussPoint *gp) : ConcreteDPM2Status(gp)
{
}


void
ConcreteDPM2RateStatus::initTempStatus()
{
    ConcreteDPM2Status::initTempStatus();
    this->tempReducedStrain = this->reducedStrain;
    this->tempKappaOne = this->kappaOne;
    this->tempKappaTwo = this->kappaTwo;
    this->tempBeta = this->beta;
    this->tempStrainRate = this->strainRate;
    this->tempRateFactor = this->rateFactor;
    this->tempKappaDTension = this->kappaDTension;
    this->tempKappaDTensionOne = this->kappaDTensionOne;
    this->tempKappaDTensionTwo = this->kappaDTensionTwo;
    this->tempDamageTension = this->damageTension;
}


void
ConcreteDPM2RateStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    // Call corresponding function of the parent class to print
    ConcreteDPM2Status::printOutputAt(file, tStep);
    fprintf(file, "kappaOne %f kappaTwo %f beta %f strainRate %f rateFactor %f\n", this->kappaOne, this->kappaTwo, this->beta, this->strainRate, this->rateFactor);
}

void
ConcreteDPM2RateStatus::updateYourself(TimeStep *tStep)
{
    ConcreteDPM2Status::updateYourself(tStep);
    this->tempReducedStrain = this->reducedStrain;
    this->kappaOne = this->tempKappaOne;
    this->kappaTwo = this->tempKappaTwo;
    this->beta = this->tempBeta;
    this->strainRate = this->tempStrainRate;
    this->rateFactor = this->tempRateFactor;
    this->kappaDTension = this->tempKappaDTension;
    this->kappaDTensionOne = this->tempKappaDTensionOne;
    this->kappaDTensionTwo = this->tempKappaDTensionTwo;
    this->damageTension = this->tempDamageTension;

}


void
ConcreteDPM2RateStatus::saveContext(DataStream &stream, ContextMode mode)
{
    ConcreteDPM2Status::saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( iores = reducedStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(kappaOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    
    if ( !stream.write(strainRate) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactor) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damageTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

void
ConcreteDPM2RateStatus::restoreContext(DataStream &stream, ContextMode mode)
{
    ConcreteDPM2Status::restoreContext(stream, mode);

    contextIOResultType iores;

    if ( ( iores = reducedStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(kappaOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(strainRate) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateFactor) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damageTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    
}


  
  //***************************
  //ConcreteDPM2Rate Class
  //***************************



FloatArrayF< 2 >
ConcreteDPM2Rate::computeDamage(const FloatArrayF< 6 > &strain,
                            const FloatMatrixF< 6, 6 > &D,
                            double deltaTime,
                            GaussPoint *gp,
                            TimeStep *tStep,
                            double tempAlpha,
                            const FloatArrayF< 6 > &effectiveStress) const
{
    auto status = static_cast< ConcreteDPM2RateStatus * >( this->giveStatus(gp) );

    //   const auto &strain = status->giveTempReducedStrain();
    //    auto principalStrain = StructuralMaterial::computePrincipalValues(from_voigt_strain(strain) );
    
    double tempEquivStrain;
    double deltaPlasticStrainNorm;
    double tempDamageTension = 0.0;

    double tempKappaDTension = 0.0;
    double tempKappaDTensionOne = 0.0, tempKappaDTensionTwo = 0.0;
    double tempKappaOne = 0.0, tempKappaTwo = 0.0;

    double tempStrainRate = 0., tempRateFactor = 1.;
    double omega = 0.0;
    double tempBeta = 0.;


    //Determine max and min value;
    double maxStrain = -1.e20, minStrain = 1.e20;
    for ( int k = 1; k <= principalStrain.giveSize(); k++ ) {
        //maximum
        if ( principalStrain.at(k) > maxStrain ) {
            maxStrain = principalStrain.at(k);
        }

    }

    //Evaluate the equivalent strains
    double oldRateStrain = status->giveRateStrain();

    tempStrainRate = ( maxStrain - oldRateStrain ) / deltaTime;
    status->letTempRateStrainBe(maxStrain);
    tempRateFactor = computeRateFactor(tempStrainRate, gp, tStep);

    double minEquivStrain = 0.;

    double sig, rho, theta;
    
    //Calculate coordinates
    computeCoordinates(effectiveStress, sig, rho, theta);

    int unAndReloadingFlag = checkForUnAndReloading(tempEquivStrain, minEquivStrain, D, gp);

    //Compute equivalent strains for  tension and compression
    double tempEquivStrainTension = 0.;
    tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;

    status->letTempRateFactorBe(rateFactor);

    double fTension = tempEquivStrainTension - status->giveKappaDTension();

    //Normalize the fs
    fTension = fTension / e0;

    double rateFactor;


    double ductilityMeasure = computeDuctilityMeasureDamage(gp, sig, rho);
    double deltaPlasticStrainNormTension, deltaPlasticStrainNormCompression;
    if ( fTension < -yieldTolDamage  ) {
        //Neither tension nor compression is active

        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        tempDamageTension = status->giveDamageTension();

        tempKappaOne = status->giveKappaOne();
        tempKappaTwo = status->giveKappaTwo();

    } else if ( fTension >= -yieldTolDamage  ) {   //Only tension is active
        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNorm = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNorm / ductilityMeasure / rateFactor;
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        omega = computeDamageParameter(tempKappaOne, tempKappaTwo, gp);
        tempKappaOne = status->giveKappaOne() + tempKappaDTension / tempRateFactor;
        tempKappaTwo = status->giveKappaTwo() + status->giveLe()( tempKappaDTensionOne + omega * tempKappaDTensionTwo ) / tempRateFactor;

        tempBeta = status->giveStrainRate()/(status->giveLe()*
            ( maxStrain - oldRateStrain ) / deltaTime);
        tempStrainRate = tempBeta*status->giveLe()*
            ( maxStrain - oldRateStrain ) / deltaTime;

        tempRateFactor = computeRateFactor(tempStrainRate, gp, tStep);

        //Initialise damage with tensile history variable
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), rateFactor);

    } else if ( fTension < -yieldTolDamage  ) {
        //Only compression is active

        //Nothing changes for the history variables in tension
        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        tempKappaOne = status->giveKappaOne();
        tempKappaTwo = status->giveKappaTwo();


        //Determine damage parameters
        tempDamageTension = status->giveDamageTension();

    } else if ( fTension >= -yieldTolDamage  ) {
        //Both tension and compression is active

        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ( ductilityMeasure * rateFactor );
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        omega = computeDamageParameter(tempKappaOne, tempKappaTwo, gp);
        tempKappaOne = status->giveKappaOne() + tempKappaDTension / tempRateFactor;
        tempKappaTwo = status->giveKappaTwo() + status->giveLe()( tempKappaDTensionOne + omega * tempKappaDTensionTwo ) / tempRateFactor;

        tempBeta = status->giveStrainRate()/(status->giveLe()*
            ( maxStrain - oldRateStrain ) / deltaTime);
        tempStrainRate = tempBeta*status->giveLe()*
            ( maxStrain - oldRateStrain ) / deltaTime;

        tempRateFactor = computeRateFactor(tempStrainRate, gp, tStep);

        //Determine the damage parameters
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), rateFactor);

    }

    //Write all temp history variables to the status
    status->letTempEquivStrainBe(tempEquivStrain);

    //Tension
    status->letTempReducedStrainBe(reducedTotalStrainVector);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->letTempKappaDTensionBe(tempKappaDTension);
    status->letTempKappaDTensionOneBe(tempKappaDTensionOne);
    status->letTempKappaDTensionTwoBe(tempKappaDTensionTwo);
    status->letTempDamageTensionBe(tempDamageTension);
    status->setTempKappaOne(tempKappaOne);
    status->setTempKappaTwo(tempKappaTwo);
    status->setTempDamage(omega);
    status->setTempRateFactor(tempRateFactor);
    status->setTempBeta(tempBeta);
    status->setTempStrainRate(tempStrainRate);



    return {
        tempDamageTension
    };
}




double
ConcreteDPM2Rate::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const
{
    double omega = 0.;

    //So that damage does not turn out to be negative if function is entered for equivstrains smaller thatn e0.
    double ftTemp = this->ft * ( 1. - yieldTolDamage );

    double wfMod = this->wf;
    double wfOneMod = this->wfOne;

    if ( this->strengthRateType > 0 ) {
        if ( this->energyRateType == 0 ) {
            wfMod /= pow(rateFactor, 2.);
            wfOneMod /= pow(rateFactor, 2.);
        } else if ( this->energyRateType == 1 ) {
            wfMod /= rateFactor;
            wfOneMod /= rateFactor;
        }
    }

    double help;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        if ( softeningType == 0 ) { //linear
            omega = ( this->eM * equivStrain * wfMod - ftTemp * wfMod + ftTemp * kappaOne * le ) /
                ( this->eM * equivStrain * wfMod - ftTemp * le * kappaTwo );
        } else if ( softeningType == 1 ) { //bilinear: Calculate damage parameter for both parts of bilinear curve  and check which fulfils limits.
            omega = ( this->eM * equivStrain * wfOneMod - ftTemp * wfOneMod - ( this->ftOne - ftTemp ) * kappaOne * le ) /
                ( this->eM * equivStrain * wfOneMod + ( this->ftOne - ftTemp ) * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help >= 0. && help < wfOneMod ) {
                return omega;
            }

            omega = ( this->eM * equivStrain * ( wfMod - wfOneMod ) - this->ftOne * ( wfMod - wfOneMod ) +
                this->ftOne * kappaOne * le  - this->ftOne * wfOneMod ) /
                    ( this->eM * equivStrain * ( wfMod - wfOneMod )  - this->ftOne * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help > wfOneMod && help < wfMod ) {
                return omega;
            }
        } else if ( softeningType == 2 ) { //exponential: Iterative solution
            omega = 1.; //Initial guess
            double residual = 0.;
            double dResidualDOmega = 0.;
            int nite = 0;

            do {
                nite++;

                residual  = ( 1 - omega ) * this->eM * equivStrain - ftTemp * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);
                dResidualDOmega = -this->eM * equivStrain + ftTemp * le * kappaTwo / wfMod * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);

                omega -= residual / dResidualDOmega;
                if ( nite > newtonIter ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(residual / ftTemp) >= 1.e-8 );
        }
    } else {
        omega = 0.;
    }


    if ( omega > 1. ) {
        omega = 1.;
    }

    if ( omega < 0. || omega < omegaOld ) {
        omega = omegaOld;
    }


    return omega;
}





double
ConcreteDPM2Rate::computeRateFactor(double strainRate,
                            GaussPoint *gp,
                            TimeStep *tStep) const
{
  double rateFactor = 1.;
  double strainRateRatio = strainRate / 1.e-6;


    if ( strainRate < 1.e-6 ) {
      rateFactor = 1.;
    } else if ( 1.e-6 <= strainRate && strainRate < 10 ) {
      rateFactor = pow(strainRateRatio, 0.018);
    } else {
      rateFactor = 0.0062 * pow(strainRateRatio, 1. / 3.);
    }
    
    
    return rateFactor;
}


int
ConcreteDPM2Rate :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    ConcreteDPM2RateStatus *status = static_cast< ConcreteDPM2RateStatus * >( this->giveStatus(gp) );
    if ( type == IST_RateFactor ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveRateFactor();
        return 1;
    } else {
        return ConcreteDPM2 :: giveIPValue(answer, gp, type, tStep);
    }

    return 1; // to make the compiler happy
}


  

MaterialStatus *
ConcreteDPM2Rate::CreateStatus(GaussPoint *gp) const
{
    return new ConcreteDPM2RateStatus(gp);
}

MaterialStatus *
ConcreteDPM2Rate::giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == nullptr ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status ) {
            gp->setMaterialStatus(status);
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}

ConcreteDPM2Rate::ConcreteDPM2RateStatus(GaussPoint *g) : ConcreteDPM2Material1Status(g), reducedStrain(), tempReducedStrain()
{
    int rsize = ConcreteDPM2Rate::giveSizeOfVoigtSymVector(gp->giveMaterialMode() );
    reducedStrain.resize(rsize);

    // reset temp vars.
    tempReducedStrain = reducedStrain;
}

}// end namespace oofem
