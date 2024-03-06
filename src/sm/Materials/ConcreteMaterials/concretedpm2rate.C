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
    this->tempBeta = this->beta;
    this->tempStrainRateTension = this->strainRateTension;
    this->tempStrainRateCompression = this->strainRateCompression;
    this->tempRateFactorTension = this->rateFactorTension;
    this->tempRateFactorCompression = this->rateFactorCompression;    
}


void
ConcreteDPM2RateStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    // Call corresponding function of the parent class to print
    ConcreteDPM2Status::printOutputAt(file, tStep);
    fprintf(file, "beta %f strainRateTension %f strainRateCompression %f rateFactorTension %f rateFactorCompression %f\n", this->beta, this->strainRateTension, this->strainRateCompression, this->rateFactorTension, this->rateFactorCompression);
}

void
ConcreteDPM2RateStatus::updateYourself(TimeStep *tStep)
{
    ConcreteDPM2Status::updateYourself(tStep);
    this->reducedStrain = this->tempReducedStrain;
    this->beta = this->tempBeta;
    this->strainRateTension = this->tempStrainRateTension;
    this->strainRateCompression = this->tempStrainRateCompression;
    this->rateFactorTension = this->tempRateFactorTension;
    this->rateFactorCompression = this->tempRateFactorCompression;
}


void
ConcreteDPM2RateStatus::saveContext(DataStream &stream, ContextMode mode)
{
    ConcreteDPM2Status::saveContext(stream, mode);

    //    contextIOResultType iores;
    contextIOResultType iores;


    if ( ( iores = reducedStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    
    if ( !stream.write(strainRateTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(strainRateCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactorTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactorCompression) ) {
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
    //    contextIOResultType iores;

    if ( !stream.read(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(strainRateTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(strainRateCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    
    if ( !stream.read(rateFactorTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateFactorCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }
  
}
  
  //***************************
  //ConcreteDPM2Rate Class
  //***************************

#define IDM_ITERATION_LIMIT 1.e-8
ConcreteDPM2Rate::ConcreteDPM2Rate(int n, Domain *d) :
    ConcreteDPM2(n, d)
{}

  
void
ConcreteDPM2Rate::initializeFrom(InputRecord &ir)
{

    ConcreteDPM2::initializeFrom(ir);
    
    this->atOne = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atOne, _IFT_ConcreteDPM2Rate_atOne);
    this->atTwo = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atTwo, _IFT_ConcreteDPM2Rate_atTwo);
    this->atThree = 0.018;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atThree, _IFT_ConcreteDPM2Rate_atThree);
    this->atFour = 0.0062;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atFour, _IFT_ConcreteDPM2Rate_atFour);
    this->atFive = 0.333;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atFive, _IFT_ConcreteDPM2Rate_atFive);

    this->acOne = 30.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acOne, _IFT_ConcreteDPM2Rate_acOne);
    this->acTwo = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acTwo, _IFT_ConcreteDPM2Rate_acTwo);
    this->acThree = 0.014;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acThree, _IFT_ConcreteDPM2Rate_acThree);
    this->acFour = 0.012;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acFour, _IFT_ConcreteDPM2Rate_acFour);
    this->acFive = 0.333;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acFive, _IFT_ConcreteDPM2Rate_acFive);  
}

  
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

    double tempEquivStrain;
    double deltaPlasticStrainNorm;
    double tempDamageTension = 0.0;
    double tempDamageCompression = 0.0;

    double tempKappaDTension = 0.0, tempKappaDCompression = 0.0;
    double tempKappaDTensionOne = 0.0, tempKappaDTensionTwo = 0.0;
    double tempKappaDCompressionOne = 0.0, tempKappaDCompressionTwo = 0.0;

    double minEquivStrain = 0.;    

    double sig, rho, theta;
    //Calculate coordinates
    computeCoordinates(effectiveStress, sig, rho, theta);

    int unAndReloadingFlag = checkForUnAndReloading(tempEquivStrain, minEquivStrain, D, gp);

    FloatArrayF<2> tempRateFactor;

    //For tempStrainRateCompression should this be (tempEquivStrain - minEquivStrain)?
    double tempStrainRateCompression = ( tempEquivStrain - status->giveEquivStrain() ) / deltaTime * this->fc / this->ft;
    double tempRateFactorCompression = computeRateFactorCompression(tempStrainRateCompression, gp, tStep);

    //Rate factor in tension has to be made mesh independent once damage has started, because the model is based on the crack band approach
    double tempStrainRateTension = 0.,tempRateFactorTension=1., tempBeta = 0.;
    if(status->giveTempDamageTension() == 0){
      //Damage is zero
      tempStrainRateTension = ( tempEquivStrain - status->giveEquivStrain() ) / deltaTime;
      tempRateFactorTension = computeRateFactorTension(tempStrainRateTension,gp, tStep);
    }
    else{
      //Damage in previous step is not zero
      tempBeta = status->giveTempBeta();
      if(tempBeta == 0){
	//Calculate tempBeta only once 
	tempBeta = status->giveStrainRateTension()/(status->giveLe()*
						    ( tempEquivStrain - status->giveEquivStrain() ) / deltaTime);
      }
      
      tempStrainRateTension = tempBeta*status->giveLe()*
	( tempEquivStrain - status->giveEquivStrain() ) / deltaTime;
      
      tempRateFactorTension = computeRateFactorTension(tempStrainRateTension, gp, tStep);	  
    }

    //Update the status here.
    status->setTempStrainRateTension(tempStrainRateTension);
    status->setTempRateFactorTension(tempRateFactorTension);
    status->setTempBeta(tempBeta);
    status->setTempStrainRateCompression(tempStrainRateCompression);
    status->setTempRateFactorCompression(tempRateFactorCompression);
   
    //Compute equivalent strains for  tension and compression
    double tempEquivStrainTension = 0.;
    double tempEquivStrainCompression = 0.;

    tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / tempRateFactorTension;

    if ( unAndReloadingFlag == 0 ) { //Standard way
      tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / tempRateFactorCompression;
    } else {
      tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / tempRateFactorCompression + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / tempRateFactorCompression;
    }

    double fTension = tempEquivStrainTension - status->giveKappaDTension();
    double fCompression = tempEquivStrainCompression - status->giveKappaDCompression();

    //Normalize the fs
    fTension = fTension / e0;
    fCompression = fCompression / e0;

    double ductilityMeasure = computeDuctilityMeasureDamage(gp, sig, rho);
    double deltaPlasticStrainNormTension, deltaPlasticStrainNormCompression;

    if ( fTension < -yieldTolDamage && fCompression < -yieldTolDamage ) {
        //Neither tension nor compression is active
        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension >= -yieldTolDamage && fCompression < -yieldTolDamage ) {   //Only tension is active
        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNorm = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNorm / ductilityMeasure * tempRateFactorTension;
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure * pow(tempRateFactorTension,2.);

        //Nothing changes for compression history variables
        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        //Initialise damage with tensile history variable
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension());

        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension < -yieldTolDamage && fCompression >= -yieldTolDamage ) {
      //Only compression is active

        //Nothing changes for the history variables in tension
        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        //Update compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression = computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression /  ductilityMeasure * tempRateFactorCompression;
	
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() + ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure * pow(tempRateFactorCompression,2.);

        //Determine damage parameters
        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression());
    } else if ( fTension >= -yieldTolDamage && fCompression >= -yieldTolDamage ) {
        //Both tension and compression is active

        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ductilityMeasure * tempRateFactorTension;
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure * pow(tempRateFactorTension,2.);

        //Update the compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression =
	  computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ductilityMeasure * tempRateFactorCompression;
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() +
	  ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure * pow(tempRateFactorCompression,2.);

        //Determine the damage parameters
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension());

        tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression());
    }

    //Write all temp history variables to the status
    status->letTempEquivStrainBe(tempEquivStrain);

    //Tension
    status->letTempEquivStrainTensionBe(tempEquivStrainTension);
    status->letTempKappaDTensionBe(tempKappaDTension);
    status->letTempKappaDTensionOneBe(tempKappaDTensionOne);
    status->letTempKappaDTensionTwoBe(tempKappaDTensionTwo);
    status->letTempDamageTensionBe(tempDamageTension);

    //Compression
    status->letTempEquivStrainCompressionBe(tempEquivStrainCompression);
    status->letTempKappaDCompressionBe(tempKappaDCompression);
    status->letTempKappaDCompressionOneBe(tempKappaDCompressionOne);
    status->letTempKappaDCompressionTwoBe(tempKappaDCompressionTwo);
    status->letTempDamageCompressionBe(tempDamageCompression);

    return {
        tempDamageTension, tempDamageCompression
    };
}

double
ConcreteDPM2Rate::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld) const
{
    double omega = 0.;

    //So that damage does not turn out to be negative if function is entered for equivstrains smaller thatn e0.
    double ftTemp = this->ft * ( 1. - yieldTolDamage );

    double help;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        if ( softeningType == 0 ) { //linear
            omega = ( this->eM * equivStrain * wf - ftTemp * wf + ftTemp * kappaOne * le ) /
                ( this->eM * equivStrain * wf - ftTemp * le * kappaTwo );
        } else if ( softeningType == 1 ) { //bilinear: Calculate damage parameter for both parts of bilinear curve  and check which fulfils limits.
            omega = ( this->eM * equivStrain * wfOne - ftTemp * wfOne - ( this->ftOne - ftTemp ) * kappaOne * le ) /
                ( this->eM * equivStrain * wfOne + ( this->ftOne - ftTemp ) * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help >= 0. && help < wfOne ) {
                return omega;
            }

            omega = ( this->eM * equivStrain * ( wf - wfOne ) - this->ftOne * ( wf - wfOne ) +
                this->ftOne * kappaOne * le  - this->ftOne * wfOne ) /
                    ( this->eM * equivStrain * ( wf - wfOne )  - this->ftOne * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help > wfOne && help < wf ) {
                return omega;
            }
        } else if ( softeningType == 2 ) { //exponential: Iterative solution
            omega = 1.; //Initial guess
            double residual = 0.;
            double dResidualDOmega = 0.;
            int nite = 0;

            do {
                nite++;

                residual  = ( 1 - omega ) * this->eM * equivStrain - ftTemp * exp(-le * ( omega * kappaTwo + kappaOne ) / wf);
                dResidualDOmega = -this->eM * equivStrain + ftTemp * le * kappaTwo / wf * exp(-le * ( omega * kappaTwo + kappaOne ) / wf);

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
ConcreteDPM2Rate::computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld) const
{
    if ( this->damageFlag == 3 ) {
        return 0.;
    }

    double ftTemp = this->ft * ( 1. - yieldTolDamage );

    double omega = 1.;
    int nite = 0;
    double residual = 0.;
    double dResidualDOmega = 0.;

    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        do {
            nite++;

            residual = ( 1. - omega ) * this->eM * equivStrain - ftTemp * exp(-( kappaOne + omega * kappaTwo ) / efCompression);
            dResidualDOmega = -this->eM * equivStrain + ftTemp * kappaTwo / efCompression * exp(-( kappaOne + omega * kappaTwo ) / efCompression);

            omega -= residual / dResidualDOmega;
            if ( nite > newtonIter ) {
                OOFEM_ERROR("algorithm not converging");
            }
        } while ( fabs(residual / ft) >= 1.e-8 );
    } else {
        omega = 0.;
    }

    if ( omega > 1. ) {
        omega = 1.;
    }
    if ( omega < omegaOld || omega < 0. ) {
        omega = omegaOld;
    }

    return omega;
}

double
ConcreteDPM2Rate::computeRateFactorTension(double strainRate,
					   GaussPoint *gp,
					   TimeStep *tStep) const
{

    //For tension according to Model Code 2010
    double rateFactorTension = 1.;
    double strainRateRatioTension = strainRate / this->atOne;
    
    if (strainRate >this->atOne && strainRate < this->atTwo ) {
      rateFactorTension = pow(strainRateRatioTension, this->atThree);
    } else {
      rateFactorTension =  this->atFour * pow(strainRateRatioTension, this->atFive);
    }
    
    return rateFactorTension;
}


double
ConcreteDPM2Rate::computeRateFactorCompression(double strainRate,
					   GaussPoint *gp,
					   TimeStep *tStep) const
{
    //For compression according to Model Code 2010
    double rateFactorCompression = 1.;
    double strainRateRatioCompression = strainRate / this->acOne ;
    if ( this->acOne < strainRate && strainRate < this->acTwo ) {
      rateFactorCompression = pow(strainRateRatioCompression, this->acThree);
    } else if ( this->acTwo < strainRate ) {
      rateFactorCompression =  this->acFour * pow(strainRateRatioCompression, this->acFive);
    }    

    return rateFactorCompression;
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

}// end namespace oofem
