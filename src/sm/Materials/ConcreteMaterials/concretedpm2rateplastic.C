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

#include "concretedpm2rateplastic.h"
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
    REGISTER_Material(ConcreteDPM2RatePlastic);

    //****************************************
    //Status of ConcreteDPM2RatePlastic
    //*************************************

    ConcreteDPM2RatePlasticStatus::ConcreteDPM2RatePlasticStatus(GaussPoint *gp) : ConcreteDPM2Status(gp)
    {}


    void
    ConcreteDPM2RatePlasticStatus::initTempStatus()
    {
        ConcreteDPM2Status::initTempStatus();
        this->tempKappaRate = this->KappaRate;
        this->tempBeta = this->beta;


    }


    void
    ConcreteDPM2RatePlasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
    {
        // Call corresponding function of the parent class to print
        ConcreteDPM2Status::printOutputAt(file, tStep);
    }

    void
    ConcreteDPM2RatePlasticStatus::updateYourself(TimeStep *tStep)
    {
        ConcreteDPM2Status::updateYourself(tStep);
        this->KappaRate = this->tempKappaRate;
        this->beta = this->tempBeta;
    }


    void
    ConcreteDPM2RatePlasticStatus::saveContext(DataStream &stream, ContextMode mode)
    {
        ConcreteDPM2Status::saveContext(stream, mode);

        contextIOResultType iores;

        if ( !stream.write(KappaRate) ) {
            THROW_CIOERR(CIO_IOERR);
        }


        if ( !stream.write(beta) ) {
            THROW_CIOERR(CIO_IOERR);
        }


    }

    void
    ConcreteDPM2RatePlasticStatus::restoreContext(DataStream &stream, ContextMode mode)
    {
        ConcreteDPM2Status::restoreContext(stream, mode);
        contextIOResultType iores;

        if ( !stream.write(KappaRate) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(beta) ) {
            THROW_CIOERR(CIO_IOERR);
        }

    }

    //***************************
    //ConcreteDPM2Rate Class
    //***************************

#define IDM_ITERATION_LIMIT 1.e-8
    ConcreteDPM2RatePlastic::ConcreteDPM2RatePlastic(int n, Domain *d) :
        ConcreteDPM2(n, d)
    {}


    void
    ConcreteDPM2RatePlastic::initializeFrom(InputRecord &ir)
    {
        ConcreteDPM2::initializeFrom(ir);

        this->cCompression = 1.e-6;
        IR_GIVE_OPTIONAL_FIELD(ir, this->cCompression, _IFT_ConcreteDPM2RatePlastic_cCompression);
        this->kappaRate0Compression = 10;
        IR_GIVE_OPTIONAL_FIELD(ir, this->kappaRate0Compression, _IFT_ConcreteDPM2RatePlastic_kappaRate0Compression);

        //this->cTension = 0.01935;
        this->cTension = 0.01935;
        IR_GIVE_OPTIONAL_FIELD(ir, this->cTension, _IFT_ConcreteDPM2RatePlastic_cTension);
        //this->Kapparate0tension = 0.00000434165;
        this->kappaRate0Tension = 0.00000434165;
        IR_GIVE_OPTIONAL_FIELD(ir, this->kappaRate0Tension, _IFT_ConcreteDPM2RatePlastic_Kapparate0tension);
    }

    FloatArrayF < 2 >
    ConcreteDPM2RatePlastic::computeDamage(const FloatArrayF < 6 > & strain,
                                           const FloatMatrixF < 6, 6 > & D,
                                           double deltaTime,
                                           GaussPoint * gp,
                                           TimeStep * tStep,
                                           double tempAlpha,
                                           const FloatArrayF < 6 > & effectiveStress) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

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

        int unAndReloadingFlag = checkForUnAndReloadingP(tempEquivStrain, minEquivStrain, D, gp);

        //        double rateFactor;
        //        if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) {
        //            rateFactor = computeRateFactor(tempAlpha);
        //        } else {
        //            rateFactor = status->giveRateFactor();
        //        }


        //Compute equivalent strains for  tension and compression
        double tempEquivStrainTension = 0.;
        double tempEquivStrainCompression = 0.;

        tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() );

        if ( unAndReloadingFlag == 0 ) { //Standard way
            tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) );
        } else {
            tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) );
        }


        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double ftYield = this->ft * ( 1 + this->cTension * log(1. + tempKappaRate / this->kappaRate0Tension) );
        double e01 = ftYield / eM;

        //Rate factor in tension has to be made mesh independent once damage has started, because the model is based on the crack band approach
        double tempBeta = 0.;
        if(status->giveTempDamageTension() == 0){
            //Damage is zero
            tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        }
        else{
            //Damage in previous step is not zero
            tempBeta = status->giveTempBeta();
            if(tempBeta == 0){
                //Calculate tempBeta only once
                tempBeta = status->giveTempKappaRate()/(status->giveLe()*
                               ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime);
            }

            tempKappaRate = tempBeta*status->giveLe()*
                ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;

        }
        //Update the status here.

        status->setTempKappaRate(tempKappaRate);
        status->setTempBeta(tempBeta);

        //        //If damage threshold is exceeded determine the rate factor from the previous step
        //        if ( ( tempEquivStrainTension > e01 || tempEquivStrainCompression > e01 ) &&
        //             ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) && !tStep->isTheFirstStep() ) {
        //            //Rate factor from last step
        //            rateFactor = status->giveRateFactor();
        //
        //            tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;
        //            if ( unAndReloadingFlag == 0 ) { //Standard way
        //                tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
        //            } else {
        //                tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
        //            }
        //        }

        //        status->letTempRateFactorBe(rateFactor);

        double fTension = tempEquivStrainTension - status->giveKappaDTension();
        double fCompression = tempEquivStrainCompression - status->giveKappaDCompression();

        //Normalize the fs
        fTension = fTension / e01;
        fCompression = fCompression / e01;

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
        } else if ( fTension >= -yieldTolDamage && fCompression < -yieldTolDamage ) { //Only tension is active
            //Update tension history variables
            tempKappaDTension = tempEquivStrainTension;
            deltaPlasticStrainNorm = computeDeltaPlasticStrainNormTensionP(tempKappaDTension, status->giveKappaDTension(), gp);
            tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNorm / ductilityMeasure;
            tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

            //Nothing changes for compression history variables
            tempKappaDCompression = status->giveKappaDCompression();
            tempKappaDCompressionOne = status->giveKappaDCompressionOne();
            tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

            //Initialise damage with tensile history variable
            this->initDamagedP(tempKappaDTension, strain, gp);

            tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), gp);

            tempDamageCompression = status->giveDamageCompression();
        } else if ( fTension < -yieldTolDamage && fCompression >= -yieldTolDamage ) {
            //Only compression is active

            //Nothing changes for the history variables in tension
            tempKappaDTension = status->giveKappaDTension();
            tempKappaDTensionOne = status->giveKappaDTensionOne();
            tempKappaDTensionTwo = status->giveKappaDTensionTwo();

            //Update compression history variables
            tempKappaDCompression = tempEquivStrainCompression;
            deltaPlasticStrainNormCompression = computeDeltaPlasticStrainNormCompressionP(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
            tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure );
            tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() + ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

            //Determine damage parameters
            tempDamageTension = status->giveDamageTension();
            tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression(), gp);
        } else if ( fTension >= -yieldTolDamage && fCompression >= -yieldTolDamage ) {
            //Both tension and compression is active

            //Update tension history variables
            tempKappaDTension = tempEquivStrainTension;
            deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTensionP(tempKappaDTension, status->giveKappaDTension(), gp);
            tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ( ductilityMeasure );
            tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

            //Update the compression history variables
            tempKappaDCompression = tempEquivStrainCompression;
            deltaPlasticStrainNormCompression =
                computeDeltaPlasticStrainNormCompressionP(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
            tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure );
            tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() +
                                       ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

            //Determine the damage parameters
            this->initDamagedP(tempKappaDTension, strain, gp);

            tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), gp);

            tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression(), gp);
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

    //    double
    //    ConcreteDPM2RatePlastic::computeRateFactor(double alpha) const
    //    {
    //        if ( this->strengthRateType == 0 ) {
    //            return 1;
    //        }
    ///*
    //        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );
    //
    //        const auto &strain = status->giveTempReducedStrain();
    //
    //        //Determine the principal values of the strain
    //
    //        auto principalStrain = StructuralMaterial::computePrincipalValues( from_voigt_strain(strain) );  ///@todo CHECK
    //
    //        //Determine max and min value;
    //        double maxStrain = -1.e20, minStrain = 1.e20;
    //        for ( int k = 1; k <= principalStrain.giveSize(); k++ ) {
    //            //maximum
    //            if ( principalStrain.at(k) > maxStrain ) {
    //                maxStrain = principalStrain.at(k);
    //            }
    //
    //            //minimum
    //            if ( principalStrain.at(k) < minStrain ) {
    //                minStrain = principalStrain.at(k);
    //            }
    //        }
    //
    //        //Evaluate the equivalent strains
    //        double strainRate;
    //        double oldRateStrain = status->giveRateStrain();
    //        if ( 1. - alpha > CDPM2_TOL ) { //Tension
    //            strainRate = ( maxStrain - oldRateStrain ) / deltaTime;
    //            status->letTempRateStrainBe(maxStrain);
    //        } else { //Compression
    //            strainRate = ( minStrain - oldRateStrain ) / deltaTime;
    //            status->letTempRateStrainBe(minStrain);
    //        }
    //
    //        //Tension
    //        //For tension according to Model Code 2010
    //        double rateFactorTension = 1.;
    //        double strainRateRatioTension = strainRate / 1.e-6;
    //
    //        if ( this->strengthRateType == 1 ) {
    //            if ( strainRate < 1.e-6 ) {
    //                rateFactorTension = 1.;
    //            } else if ( 1.e-6 < strainRate ) {
    //                rateFactorTension = pow(strainRateRatioTension, 0.018);
    //            }
    //        } else if ( this->strengthRateType == 2 ) {
    //            if ( strainRate < 1.e-6 ) {
    //                rateFactorTension = 1.;
    //            } else if ( 1.e-6 < strainRate && strainRate < 10 ) {
    //                rateFactorTension = pow(strainRateRatioTension, 0.018);
    //            } else {
    //                rateFactorTension =  0.0062 * pow(strainRateRatioTension, 1. / 3.);
    //            }
    //        }
    //
    //        //For compression according to Model Code 2010
    //        double rateFactorCompression = 1.;
    //        double strainRateRatioCompression = strainRate / ( -30.e-6 );
    //        if ( this->strengthRateType == 1 ) {
    //            if ( strainRate > -30.e-6 ) {
    //                rateFactorCompression = 1.;
    //            } else if ( -30.e-6 > strainRate ) {
    //                rateFactorCompression = pow(strainRateRatioCompression, 0.014);
    //            }
    //        } else if ( this->strengthRateType == 2 ) {
    //            if ( strainRate > -30.e-6 ) {
    //                rateFactorCompression = 1.;
    //            } else if ( -30.e-6 > strainRate && strainRate > -30 ) {
    //                rateFactorCompression = pow(strainRateRatioCompression, 0.014);
    //            } else if ( -30 > strainRate && strengthRateType == 2 ) {
    //                rateFactorCompression =  0.012 * pow(strainRateRatioCompression, 0.333);
    //            }
    //        }
    //*/
    //        //double rateFactor = ( 1. - alpha ) * rateFactorTension + alpha * rateFactorCompression;
    //        double rateFactor = 1;
    //
    //        return rateFactor;
    //    }






    double
    ConcreteDPM2RatePlastic::computeDeltaPlasticStrainNormTensionP(double tempKappaD, double kappaD, GaussPoint *gp) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        const auto &tempPlasticStrain = status->giveTempPlasticStrain();
        const auto &plasticStrain = status->givePlasticStrain();

        auto deltaPlasticStrain = tempPlasticStrain - plasticStrain;

        double deltaPlasticStrainNorm = 0;

        //Distinguish pre-peak, peak and post-peak

        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );
        double e01 = ftYield / eM;

        double factor = 0.;
        if ( tempKappaD < e01 * ( 1. - yieldTolDamage ) ) {
            deltaPlasticStrainNorm = 0.;
        } else if ( tempKappaD > e01 * ( 1. - yieldTolDamage ) && kappaD < e01  * ( 1. - yieldTolDamage ) ) {
            factor = ( 1. - ( e01 - kappaD ) / ( tempKappaD - kappaD ) );
            deltaPlasticStrain *= factor;
            deltaPlasticStrainNorm = norm(deltaPlasticStrain);
        } else {
            deltaPlasticStrainNorm = norm(deltaPlasticStrain);
        }

        return deltaPlasticStrainNorm;
    }

    double
    ConcreteDPM2RatePlastic::computeDeltaPlasticStrainNormCompressionP(double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp, const double rho) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        const auto &tempPlasticStrain = status->giveTempPlasticStrain();
        const auto &plasticStrain = status->givePlasticStrain();

        auto deltaPlasticStrain = tempAlpha * ( tempPlasticStrain - plasticStrain );

        double deltaPlasticStrainNorm = 0;


        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );
        double e01 = ftYield / eM;

        //Distinguish pre-peak, peak and post-peak
        if ( tempKappaD < e01 * ( 1. - yieldTolDamage ) ) {
            deltaPlasticStrainNorm = 0.;
        } else if ( tempKappaD > e01 * ( 1. - yieldTolDamage ) && kappaD < e01 * ( 1. - yieldTolDamage ) ) {
            double factor = ( 1. - ( e01 - kappaD ) / ( tempKappaD - kappaD ) );
            deltaPlasticStrain *= factor;
            deltaPlasticStrainNorm = norm(deltaPlasticStrain);
        } else {
            deltaPlasticStrainNorm = norm(deltaPlasticStrain);
        }

        double tempKappaP = status->giveTempKappaP();
        double yieldHardTwo = computeHardeningTwo(tempKappaP);
        double extraFactor;
        if ( rho < 1.e-16 ) {
            extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / 1.e-16 / sqrt(1. + 2. * pow(this->dilationConst, 2.) );
        } else {
            extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / rho / sqrt(1. + 2. * pow(this->dilationConst, 2.) );
        }

        return deltaPlasticStrainNorm * extraFactor;
    }

    int
    ConcreteDPM2RatePlastic::checkForUnAndReloadingP(double &tempEquivStrain, double &minEquivStrain, const FloatMatrixF < 6, 6 > &D, GaussPoint *gp) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        //Access old and new strains
        const auto &oldStrain = status->giveReducedStrain();
        const auto &strain = status->giveTempReducedStrain();

        //Compute the temp equivalent strain
        auto tempElasticStrain = strain - status->giveTempPlasticStrain();
        auto tempEffectiveStress = dot(D, tempElasticStrain);

        double sigEffective, rhoEffective, thetaEffective;
        computeCoordinates(tempEffectiveStress, sigEffective, rhoEffective, thetaEffective);
        tempEquivStrain = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective, gp);
        //Get the equivalent strain from the status
        double equivStrain = status->giveEquivStrain();

        //Compute the increment of effective stress
        auto elasticStrain = oldStrain - status->givePlasticStrain();
        auto effectiveStress = dot(D, elasticStrain);

        auto deltaEffectiveStress = tempEffectiveStress - effectiveStress;

        //Compute equivalent strains for stress state slightly greater than the effective stress and smaller than the temp effective stress
        //For slightly more than effective stress
        auto intermediateEffectiveStressPlus = effectiveStress + 0.01 * deltaEffectiveStress;
        computeCoordinates(intermediateEffectiveStressPlus, sigEffective, rhoEffective, thetaEffective);
        double equivStrainPlus = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective, gp);

        //For slightly less than temp effective stress
        auto intermediateEffectiveStressMinus = effectiveStress + 0.99 * deltaEffectiveStress;
        computeCoordinates(intermediateEffectiveStressMinus, sigEffective, rhoEffective, thetaEffective);
        double tempEquivStrainMinus = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective, gp);

        //Check for unloading and reloading in the same step
        int unloadingFlag = 0;
        minEquivStrain = equivStrain;

        if ( ( equivStrain > equivStrainPlus && tempEquivStrain > tempEquivStrainMinus ) &&
             ( fabs(equivStrainPlus - equivStrain) > yieldTolDamage / 100. && fabs(tempEquivStrainMinus - tempEquivStrain) > yieldTolDamage / 100. ) ) {
            unloadingFlag = 1;
            //Unloading and reloading takes place. Find the minimum equivalent strain by subincrementing the effective stress increment
            for ( double k = 1.0; k <= 100.0; k = k + 1.0 ) {
                auto intermediateEffectiveStress = effectiveStress + k / 100. * deltaEffectiveStress;
                computeCoordinates(intermediateEffectiveStress, sigEffective, rhoEffective, thetaEffective);
                double midEquivStrain = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective, gp);

                if ( midEquivStrain <= minEquivStrain ) {
                    minEquivStrain = midEquivStrain;
                } else {
                    return unloadingFlag;
                }
            }
        }
        return unloadingFlag;
    }



    double
    ConcreteDPM2RatePlastic::computeEquivalentStrainP(double sig, double rho, double theta, GaussPoint *gp) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );
        //Xiaowei: This seems to be wrong! You need to calculate the equivalent strain based on ftYield and fcYield.initDamaged
        // This will affect fc but also m. The equivalent strain is based on the yield function.

        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double fcYield = this->fc * ( 1 + cCompression * log(1. + tempKappaRate / kappaRate0Compression) );
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );

        double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );



        double rFunction = ( 4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(theta), 2.) + pow(2. * this->ecc - 1., 2.) ) / ( 2. * ( 1. - pow(this->ecc, 2.) ) * cos(theta) + ( 2. * this->ecc - 1. ) * sqrt(4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(theta), 2.) + 5. * pow(this->ecc, 2.) - 4. * this->ecc) );

        double pHelp = -myield * ( rho * rFunction / ( sqrt(6.) * fcYield ) + sig / fcYield );

        double qHelp = -3. / 2. * pow(rho, 2.) / pow(fcYield, 2.);

        double help = -0.5 * pHelp + sqrt(pow(pHelp, 2.) / 4. - qHelp);

        double e01 = ftYield / eM;


        double tempEquivStrain = 0.;
        if ( help > 0 ) {
            tempEquivStrain = help * e01;
        }
        return tempEquivStrain;
    }

    double
    ConcreteDPM2RatePlastic::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, GaussPoint *gp) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        double omega = 0.;

        //So that damage does not turn out to be negative if function is entered for equivstrains smaller thatn e0.


        double wfMod = this->wf;
        double wfOneMod = this->wfOne;

        //        if ( this->strengthRateType > 0 ) {
        //            if ( this->energyRateType == 0 ) {
        //                wfMod /= pow(rateFactor, 2.);
        //                wfOneMod /= pow(rateFactor, 2.);
        //            } else if ( this->energyRateType == 1 ) {
        //                wfMod /= rateFactor;
        //                wfOneMod /= rateFactor;
        //            }
        //        }


        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );
        double e01 = ftYield / this->eM;

        double ftTemp = ftYield * ( 1. - yieldTolDamage );


        double help;
        if ( equivStrain > e01 * ( 1. - yieldTolDamage ) ) {
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
    ConcreteDPM2RatePlastic::computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld, GaussPoint *gp) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        if ( this->damageFlag == 3 ) {
            return 0.;
        }

        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );
        double e01 = ftYield / eM;


        double ftTemp = ftYield * ( 1. - yieldTolDamage );
        double efCompressionMod = this->efCompression;

        //        if ( this->strengthRateType > 0 ) {
        //            if ( this->energyRateType == 0 ) {
        //                efCompressionMod /= pow(rateFactor, 2.);
        //            } else if ( this->energyRateType == 1 ) {
        //                efCompressionMod /= rateFactor;
        //            }
        //        }

        double omega = 1.;
        int nite = 0;
        double residual = 0.;
        double dResidualDOmega = 0.;



        if ( equivStrain > e01 * ( 1. - yieldTolDamage ) ) {
            do {
                nite++;

                residual = ( 1. - omega ) * this->eM * equivStrain - ftTemp * exp(-( kappaOne + omega * kappaTwo ) / efCompressionMod);
                dResidualDOmega = -this->eM * equivStrain + ftTemp * kappaTwo / efCompressionMod * exp(-( kappaOne + omega * kappaTwo ) / efCompressionMod);

                omega -= residual / dResidualDOmega;
                if ( nite > newtonIter ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(residual / ftYield) >= 1.e-8 );
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


    void
    ConcreteDPM2RatePlastic::initDamagedP(double kappaD, const FloatArrayF < 6 > &strain, GaussPoint *gp) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );


        double tempKappaRate = ( status->giveTempKappaP() - status->giveKappaP() ) / deltaTime;
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );
        double e01 = ftYield / eM;



        if ( kappaD <= e01 * ( 1. - yieldTolDamage ) ) {
            return;
        }


        if ( helem > 0. ) {
            status->setLe(helem);
        } else if ( status->giveDamageTension() == 0. && status->giveDamageCompression() == 0. ) {
            //auto [principalStrains, principalDir] = computePrincipalValDir(from_voigt_strain(strain)); // c++17
            auto tmp = computePrincipalValDir(from_voigt_strain(strain) );
            auto principalStrains = tmp.first;
            auto principalDir = tmp.second;

            // find index of max positive principal strain
            int indx = 1;
            for ( int i = 2; i <= 3; i++ ) {
                if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                    indx = i;
                }
            }

            FloatArray crackPlaneNormal(3);
            for ( int i = 1; i <= 3; i++ ) {
                crackPlaneNormal.at(i) = principalDir.at(i, indx);
            }

            // evaluate the projected element size
            double le = gp->giveElement()->giveCharacteristicLength(crackPlaneNormal);
            if ( le == 0. ) {
                le = gp->giveElement()->computeMeanSize();
            }

            // store le in the corresponding status
            status->setLe(le);
        } else if ( status->giveLe() == 0. ) {
            // this happens if the status is initialized from a file
            // with nonzero damage
            // le determined as square root of element area or cube root of el. volume
            double le = gp->giveElement()->computeMeanSize();
            status->setLe(le);
        }
    }


    FloatArrayF < 6 >
    ConcreteDPM2RatePlastic::performPlasticityReturn(GaussPoint * gp, const FloatMatrixF < 6, 6 > & D, const FloatArrayF < 6 > & strain, const double deltaTime) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        ConcreteDPM2_ReturnResult returnResult = RR_Unknown;
        ConcreteDPM2_ReturnType returnType = RT_Unknown;

        //get plastic strain and kappa
        auto tempPlasticStrain = status->givePlasticStrain();
        double tempKappaP = status->giveKappaP();

        //this theta computed here should stay constant for the rest of procedure.
        const auto &oldStrain = status->giveReducedStrain();

        // introduce a strange subincrementation flag
        int subIncrementFlag = 0;

        double apexStress = 0.;
        int subincrementcounter = 0;
        //Attempt to implement subincrementation
        // initialize variables
        subIncrementFlag = 0;
        auto convergedStrain = oldStrain;
        auto tempStrain = strain;
        auto deltaStrain = strain - oldStrain;

        FloatArrayF < 6 > effectiveStress;

        //To get into the loop
        returnResult = RR_NotConverged;
        while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
            auto elasticStrain = tempStrain - tempPlasticStrain;

            effectiveStress = dot(D, elasticStrain);

            double sig, rho, theta;
            computeCoordinates(effectiveStress, sig, rho, theta);
            double yieldValue = computeYieldValue(sig, rho, theta, tempKappaP, deltaTime, gp);

            apexStress = 0.;

            if ( yieldValue > 0. ) {
                checkForVertexCase(apexStress, returnType, sig, tempKappaP, gp);
                if ( returnType == RT_Tension || returnType == RT_Compression ) {
                    tempKappaP = performVertexReturn(effectiveStress, returnResult, returnType, apexStress, tempKappaP, gp, deltaTime);
                    status->letTempKappaPBe(tempKappaP);
                    if ( returnType == RT_Tension ) {
                        status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_VertexTension);
                    } else if ( returnType == RT_Compression ) {
                        status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_VertexCompression);
                    }
                }
                if ( returnType == RT_Regular ) {
                    tempKappaP = performRegularReturn(effectiveStress, returnResult, returnType, tempKappaP, gp, theta, deltaTime);
                    status->letTempKappaPBe(tempKappaP);
                }
            } else {
                returnResult = RR_Converged;
                tempPlasticStrain = status->givePlasticStrain();
                status->letTempPlasticStrainBe(tempPlasticStrain);
                status->letTempKappaPBe(tempKappaP);
                break;
            }

            if ( returnResult == RR_NotConverged ) {
                subincrementcounter++;
                if ( subincrementcounter > 10 ) {
                    OOFEM_LOG_INFO("Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                    OOFEM_LOG_INFO("Old strain vector %g %g %g %g %g %g  \n", oldStrain.at(1), oldStrain.at(2), oldStrain.at(3), oldStrain.at(4), oldStrain.at(5), oldStrain.at(6) );

                    const auto &help = status->giveTempPlasticStrain();
                    OOFEM_LOG_INFO("Old plastic strain vector %g %g %g %g %g %g  \n", help.at(1), help.at(2), help.at(3), help.at(4), help.at(5), help.at(6) );
                    OOFEM_LOG_INFO("New strain vector %g %g %g %g %g %g  \n", strain.at(1), strain.at(2), strain.at(3), strain.at(4), strain.at(5), strain.at(6) );

                    computeCoordinates(effectiveStress, sig, rho, theta);
                    double sig1, rho1, theta1;
                    auto help1 = dot(D, oldStrain - help);
                    computeCoordinates(help1, sig1, rho1, theta1);
                    yieldValue = computeYieldValue(sig, rho, theta, tempKappaP, deltaTime, gp);
                    OOFEM_LOG_INFO("OLD Sig %g rho %g theta %g  \n", sig1, rho1, theta1);
                    OOFEM_LOG_INFO("NEW Sig %g rho %g theta %g  \n", sig, rho, theta);
                    if ( returnType == RT_Tension || returnType == RT_Compression ) {
                        OOFEM_LOG_INFO("Vertex case apexstress %g\n", apexStress);
                    } else {
                        OOFEM_LOG_INFO("Regular case %g \n", 15.18);
                    }
                    OOFEM_LOG_INFO("KappaP old %g new %g yieldfun %g\n", status->giveTempKappaP(), tempKappaP, yieldValue);
                    OOFEM_WARNING("ConcreteDamagePlasticity2:: performPlasticityReturn: Could not reach convergence with small deltaStrain, giving up. Delete Element number %d", gp->giveElement()->giveNumber() );
                    //OOFEM_ERROR("Could not reach convergence with small deltaStrain, giving up.");
                    status->setTempDeletionFlag(1);
                    for (int k = 0; k < 6; k++) {
                        effectiveStress.at(k + 1) = 0.;
                    }
                    return effectiveStress;
                } else if ( subincrementcounter > 9 && tempKappaP < 1. ) {
                    tempKappaP = 1.;
                    status->letTempKappaPBe(tempKappaP);
                }

                subIncrementFlag = 1;
                deltaStrain *= 0.5;
                tempStrain = convergedStrain + deltaStrain;
            } else if ( returnResult == RR_Converged && subIncrementFlag == 0 ) {
                auto C = inv(D); // compliance
                elasticStrain = dot(C, effectiveStress);
                tempPlasticStrain = strain - elasticStrain;
                status->letTempPlasticStrainBe(tempPlasticStrain);
            } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
                subincrementcounter = 0;
                auto C = inv(D); // compliance
                elasticStrain = dot(C, effectiveStress);
                tempPlasticStrain = tempStrain - elasticStrain;
                status->letTempPlasticStrainBe(tempPlasticStrain);

                subIncrementFlag = 0;
                returnResult = RR_NotConverged;
                convergedStrain = tempStrain;
                deltaStrain = strain - convergedStrain;
                tempStrain = strain;
            }
        }

        return effectiveStress;
    }

    double
    ConcreteDPM2RatePlastic::performVertexReturn(FloatArrayF < 6 > &effectiveStress,
                                                 ConcreteDPM2_ReturnResult &returnResult,
                                                 ConcreteDPM2_ReturnType &returnType,
                                                 double apexStress, double tempKappaP,
                                                 GaussPoint *gp, double deltaTime) const
    {
        //auto [deviatoricStressTrial, sigTrial] = computeDeviatoricVolumetricSplit(effectiveStress); // c++17
        auto tmp = computeDeviatoricVolumetricSplit(effectiveStress);
        auto deviatoricStressTrial = tmp.first;
        auto sigTrial = tmp.second;

        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        double rhoTrial = computeSecondCoordinate(deviatoricStressTrial);

        double kappaInitial = tempKappaP;

        double sig2 = apexStress;

        tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial);

        double yieldValue = computeYieldValue(sigTrial, 0., 0., tempKappaP, deltaTime, gp);

        tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2);

        double yieldValueMid = computeYieldValue(sig2, 0., 0., tempKappaP, deltaTime, gp);

        if ( yieldValue * yieldValueMid >= 0. ) {
            returnType = RT_Regular;
            returnResult = RR_NotConverged;
            return kappaInitial;
        }

        double dSig, sigAnswer;
        if ( yieldValue < 0.0 ) {
            dSig = sig2 - sigTrial;
            sigAnswer = sig2;
        } else {
            dSig = sigTrial - sig2;
            sigAnswer = sig2;
        }

        for ( int j = 0; j < 250; j++ ) {
            dSig = 0.5 * dSig;

            double sigMid = sigAnswer + dSig;


            tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigMid);

            yieldValueMid = computeYieldValue(sigMid, 0., 0., tempKappaP, deltaTime, gp);

            if ( yieldValueMid <= 0. ) {
                sigAnswer = sigMid;
            }

            if ( fabs(yieldValueMid) < yieldTol && yieldValueMid <= 0. ) {
                double ratioPotential =
                    computeRatioPotential(sigAnswer,  0, tempKappaP);


                double ratioTrial = rhoTrial / ( sigTrial - sigAnswer );

                if ( ( ( ( ratioPotential >= ratioTrial ) && returnType == RT_Tension ) ) ||
                     ( ( ratioPotential <= ratioTrial ) && returnType == RT_Compression ) ) {
                    for ( int i = 0; i < 3; i++ ) {
                        effectiveStress.at(i + 1) = sigAnswer;
                    }

                    for ( int i = 3; i < 6; i++ ) {
                        effectiveStress.at(i + 1) = 0.;
                    }
                    returnResult = RR_Converged;
                    return tempKappaP;
                } else {
                    returnType = RT_Regular;
                    returnResult = RR_NotConverged;
                    return kappaInitial;
                }
            }
        }

        for ( int i = 0; i < 3; i++ ) {
            effectiveStress.at(i + 1) = sigAnswer;
        }

        for ( int i = 3; i < 6; i++ ) {
            effectiveStress.at(i + 1) = 0.;
        }
        returnResult = RR_Converged;

        OOFEM_WARNING("Perform vertex return not converged!\n");

        return tempKappaP;
    }


    double
    ConcreteDPM2RatePlastic::performRegularReturn(FloatArrayF < 6 > &effectiveStress,
                                                  ConcreteDPM2_ReturnResult &returnResult,
                                                  ConcreteDPM2_ReturnType &returnType,
                                                  double kappaP,
                                                  GaussPoint *gp,
                                                  double theta,
                                                  const double deltaTime) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        //Define stressVariables
        double trialSig, trialRho;

        auto trialStress = effectiveStress;

        //compute invariants from stress state
        //auto [deviatoricTrialStress, trialSig] = computeDeviatoricVolumetricSplit(trialStress); // c++17
        auto tmp = computeDeviatoricVolumetricSplit(trialStress);
        auto deviatoricTrialStress = tmp.first;
        trialSig = tmp.second;
        trialRho = computeSecondCoordinate(deviatoricTrialStress);

        double sig = trialSig;
        double rho = trialRho;

        // Starting guess:
        double tempKappaP = kappaP;

        // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
        double yieldValue = computeYieldValue(sig, rho, theta, tempKappaP, deltaTime, gp);

        //initialise unknowns
        FloatArray unknowns;
        FloatArray residuals;
        residuals.resize(4);
        residuals.at(4) = yieldValue; //store in the last element of the array
        unknowns.resize(4);
        unknowns.at(1) = trialSig;
        unknowns.at(2) = trialRho;
        unknowns.at(3) = tempKappaP;
        unknowns.at(4) = 0.;

        double deltaLambda = 0.;
        double normOfResiduals  = 1.;//just to get into the loop

        // N.R. iteration for finding the correct plastic return which is found when the norm of the residuals are equal to zero

        int iterationCount = 0;
        while ( normOfResiduals > yieldTol   ) {
            iterationCount++;
            if ( iterationCount == newtonIter ) {
                returnResult = RR_NotConverged;

                return kappaP;
            }

            auto residualsNorm = residuals;
            //Normalize residuals. Think about it more.
            residualsNorm.at(1) /= this->kM;
            residualsNorm.at(2) /= 2. * this->gM;

            normOfResiduals = norm(residualsNorm);

            printf("normOfResiduals = %e\n", normOfResiduals);

            if ( std::isnan(normOfResiduals) ) {
                returnResult = RR_NotConverged;
                return kappaP;
            }

            if ( normOfResiduals > yieldTol  ) {
                // Test to run newton iteration using inverse of Jacobian
                auto jacobian = computeJacobian(sig, rho, theta, tempKappaP, deltaLambda, gp, deltaTime);

                try {
                    auto deltaIncrement = solve( jacobian, FloatArrayF < 4 > ( residuals ) );
                    unknowns -= deltaIncrement;
                } catch(...) {
                    returnResult = RR_NotConverged;
                    return kappaP;
                }

                unknowns.at(2) = max(unknowns.at(2), 0.); //Keep rho greater than zero!
                unknowns.at(3) = max(unknowns.at(3), kappaP); //Keep deltaKappa greater than zero!
                unknowns.at(4) = max(unknowns.at(4), 0.); //Keep deltaLambda greater than zero!

                //compute residuals
                sig = unknowns.at(1);
                rho = unknowns.at(2);
                tempKappaP = unknowns.at(3);
                deltaLambda = unknowns.at(4);

                /* Compute the mVector holding the derivatives of the g function and the hardening function*/
                auto dGDInv = computeDGDInv(sig, rho, tempKappaP);
                double dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, theta, tempKappaP);

                residuals.at(1) = sig - trialSig + this->kM * deltaLambda * dGDInv.at(1);
                residuals.at(2) = rho - trialRho + ( 2. * this->gM ) * deltaLambda * dGDInv.at(2);
                residuals.at(3) = -tempKappaP + kappaP + deltaLambda * dKappaDDeltaLambda;
                residuals.at(4) = computeYieldValue(sig, rho, theta, tempKappaP, deltaTime, gp);
            }
        }


        //compute the principal directions of the stress
        //auto [helpStress, stressPrincipalDir] = StructuralMaterial :: computePrincipalValDir(from_voigt_stress(trialStress)); // c++17
        auto tmpEig = StructuralMaterial::computePrincipalValDir( from_voigt_stress(trialStress) );
        auto stressPrincipalDir = tmpEig.second;

        FloatArrayF < 6 > stressPrincipal;
        stressPrincipal [ 0 ] = sig + sqrt(2. / 3.) * rho * cos(theta);
        stressPrincipal [ 1 ] = sig + sqrt(2. / 3.) * rho * cos(theta - 2. * M_PI / 3.);
        stressPrincipal [ 2 ] = sig + sqrt(2. / 3.) * rho * cos(theta + 2. * M_PI / 3.);
        effectiveStress = transformStressVectorTo(stressPrincipalDir, stressPrincipal, 1);
        returnResult = RR_Converged;

        //Store deltaLambda in status
        status->letDeltaLambdaBe(deltaLambda);

        return tempKappaP;
    }




    FloatMatrixF < 4, 4 >
    ConcreteDPM2RatePlastic::computeJacobian(double sig,
                                             double rho,
                                             double theta,
                                             double kappa,
                                             double deltaLambda,
                                             GaussPoint * gp,
                                             const double deltaTime) const
    {
        auto dFDInv = computeDFDInv(sig, rho, theta, kappa, gp);
        auto dGDInv = computeDGDInv(sig, rho, kappa);
        auto dDGDDInv = computeDDGDDInv(sig, rho, kappa);


        //Calculate time increment



        double dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, theta, kappa);
        double dFDKappa = computeDFDKappa(sig, rho, theta, kappa, deltaTime, gp);

        auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, kappa);

        double dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, theta, kappa);
        auto dDKappaDDeltaLambdaDInv = computeDDKappaDDeltaLambdaDInv(sig, rho, theta, kappa);

        FloatMatrixF < 4, 4 > answer;
        /* Compute matrix*/
        answer.at(1, 1) = 1. + this->kM * deltaLambda * dDGDDInv.at(1, 1);
        answer.at(1, 2) = this->kM * deltaLambda * dDGDDInv.at(1, 2);
        answer.at(1, 3) = this->kM * deltaLambda * dDGDInvDKappa.at(1);
        answer.at(1, 4) = this->kM * dGDInv.at(1);
        /**/
        answer.at(2, 1) = 2. * this->gM * deltaLambda * dDGDDInv.at(2, 1);
        answer.at(2, 2) = 1. + 2. * this->gM * deltaLambda * dDGDDInv.at(2, 2);
        answer.at(2, 3) = 2. * this->gM * deltaLambda * dDGDInvDKappa.at(2);
        answer.at(2, 4) = 2. * this->gM * dGDInv.at(2);
        /**/
        answer.at(3, 1) = deltaLambda * dDKappaDDeltaLambdaDInv.at(1);
        answer.at(3, 2) = deltaLambda * dDKappaDDeltaLambdaDInv.at(2);
        answer.at(3, 3) = deltaLambda * dDKappaDDeltaLambdaDKappa - 1.;
        answer.at(3, 4) = dKappaDDeltaLambda;
        /**/
        answer.at(4, 1) = dFDInv.at(1);
        answer.at(4, 2) = dFDInv.at(2);
        answer.at(4, 3) = dFDKappa;
        answer.at(4, 4) = 0.;
        return answer;
    }







    double
    ConcreteDPM2RatePlastic::computeYieldValue(double sig,
                                               double rho,
                                               double theta,
                                               double tempKappa,
                                               double deltaTime,
                                               GaussPoint *gp) const
    {
        //compute yieldHard
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );


        double yieldHardOne = computeHardeningOne(tempKappa);
        double yieldHardTwo = computeHardeningTwo(tempKappa);

        //Xiaowei: Check this. Your way of calculating tempKappaRate might have caused problems. Now it uses giveKappaP() and the tempKappa that is input of the function. The same in dfDKappa.
        double tempKappaRate = ( tempKappa - status->giveKappaP() ) / deltaTime;

        //  printf("tempKappaRate = %e\n", tempKappaRate);

        //  compute elliptic function r
        double rFunction = ( 4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) +
                             pow( ( 2. * ecc - 1. ), 2. ) ) /
                           ( 2. * ( 1. - pow(ecc, 2.) ) * cos(theta) +
                             ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.)
                                                      + 5. * pow(ecc, 2.) - 4. * ecc) );

        double fcYield = this->fc * ( 1 + cCompression * log(1. + tempKappaRate / kappaRate0Compression) );
        double ftYield = this->ft * ( 1 + cTension * log(1. + tempKappaRate / kappaRate0Tension) );

        double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );

        //compute help function Al
        double Al = ( 1. - yieldHardOne ) * pow( ( sig / fcYield + rho / ( sqrt(6.) * fcYield ) ), 2. ) +
                    sqrt(3. / 2.) * rho / fcYield;


        //Compute yield equation
        return pow(Al, 2.) +
               pow(yieldHardOne, 2.) * yieldHardTwo * myield * ( sig / fcYield + rho * rFunction / ( sqrt(6.) * fcYield ) ) -
               pow(yieldHardOne, 2.) * pow(yieldHardTwo, 2.);
    }


    //    Xiaowei, you have extended computeDfDkappa to take more variables. Therefore, everywhere where computeDFDKappa is called you need to make the schange, too. For instance, computeJacobian, computeFullJacobian. In those funtions, you will need to access deltaTime.

    double
    ConcreteDPM2RatePlastic::computeDFDKappa(double sig,
                                             double rho,
                                             double theta,
                                             double tempKappa,
                                             double deltaTime,
                                             GaussPoint *gp) const
    {
        double dFDKappa;

        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );


        //compute yieldHard and yieldSoft
        double yieldHardOne = computeHardeningOne(tempKappa);
        double yieldHardTwo = computeHardeningTwo(tempKappa);
        // compute the derivative of the hardening and softening laws
        double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
        double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);


        double tempKappaRate = ( tempKappa - status->giveKappaP() ) / deltaTime;


        double fcYield = fc * ( 1 + cCompression * log(1 + tempKappaRate / kappaRate0Compression) );
        double ftYield = ft * ( 1 + cTension * log(1 + tempKappaRate / kappaRate0Tension) );
        double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );

        //compute elliptic function r
        double rFunction =
            ( 4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + pow( ( 2. * ecc - 1. ), 2 ) ) /
            ( 2 * ( 1. - pow(ecc, 2) ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + 5. * pow(ecc, 2) - 4. * ecc) );

        //compute help functions Al, Bl
        double Al = ( 1. - yieldHardOne ) * pow( ( sig / fcYield + rho / ( sqrt(6.) * fcYield ) ), 2.) + sqrt(3. / 2.) * rho / fcYield;


        double Bl = sig / fcYield + rho / ( fcYield * sqrt(6.) );
        double dFDYieldHardOne = -2. * Al * pow(Bl, 2.)
                                 + 2. * yieldHardOne * yieldHardTwo * myield * ( sig / fcYield + rho * rFunction / ( sqrt(6.) * fcYield ) ) - 2. * yieldHardOne * pow(yieldHardTwo, 2.);

        double dFDYieldHardTwo = pow(yieldHardOne, 2.) * myield * ( sig / fcYield + rho * rFunction / ( sqrt(6.) * fcYield ) ) - 2. * yieldHardTwo * pow(yieldHardOne, 2.);



        //Xiaowei, if possible can you give them names so that we know which these terms are? Something like "double dFcDkappaRate" or equivalent.
        //Also, follow the OOFEM coding rules. New word in names starts with captial letter. So, "fcYield, ..." For variables, small letters at the start. There is a document on the OOFEM webpage which explains everything.

        double dFtdKapparate = ft * cTension * 1 / ( tempKappaRate + kappaRate0Tension );
        double dFcdKapparate = fc * cCompression * 1 / ( tempKappaRate + kappaRate0Compression );

        double dm0dFt = ( ( -3 * ecc ) / ( ecc + 1 ) ) * ( ( ftYield * ftYield + fcYield * fcYield ) / ( fcYield * ftYield * ftYield ) );
        double dm0dFc = ( ( 3 * ecc ) / ( ecc + 1 ) ) * ( ( ftYield * ftYield + fcYield * fcYield ) / ( fcYield * fcYield * ftYield ) );

        double dFdFt = dm0dFt * pow(yieldHardOne, 2.) * yieldHardTwo * ( sig / fcYield + rho * rFunction / ( sqrt(6.) * fcYield ) );

        double dBldFc = -rho / ( sqrt(6.) * fcYield * fcYield ) - sig / ( fcYield * fcYield );

        double dFdFc = 2 * Al * ( 2 * ( 1 - yieldHardOne ) * Bl * dBldFc - sqrt(3. / 2.) * rho / ( fcYield * fcYield ) ) + dm0dFc * pow(yieldHardOne, 2.) * yieldHardTwo * ( sig / fcYield + rho * rFunction / ( sqrt(6.) * fcYield ) )
                       - pow(yieldHardOne, 2.) * yieldHardTwo * myield / ( fcYield * fcYield ) * ( sig  + rho * rFunction / sqrt(6.) );



        // compute dFDKappa
        dFDKappa = dFDYieldHardOne * dYieldHardOneDKappa + dFDYieldHardTwo * dYieldHardTwoDKappa
                   + ( dFdFt * dFtdKapparate + dFdFc * dFcdKapparate ) / deltaTime;

        //dFDKappa = dFDYieldHardOne * dYieldHardOneDKappa + dFDYieldHardTwo * dYieldHardTwoDKappa;

        /*
         * set dFDKappa to zero, if it becomes greater than zero.
         * dFDKappa can only be negative or zero in the converged state for
         * the case of hardenig and perfect plasticity. For trial stresses, however,
         * it might be psoitive, which may lead to convergency problems. Therefore,
         * it is set to zero in this cases.
         */
        if ( dFDKappa > 0. ) {
            dFDKappa = 0.;
        }

        return dFDKappa;
    }




    FloatArrayF < 2 >
    ConcreteDPM2RatePlastic::computeDFDInv(double sig, double rho, double theta, double tempKappa, GaussPoint * gp) const
    {

        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );

        double tempKappaRate = ( tempKappa - status->giveKappaP() ) / deltaTime;
        double fcYield = fc * ( 1 + cCompression * log(1 + tempKappaRate / kappaRate0Compression) );
        double ftYield = ft * ( 1 + cTension * log(1 + tempKappaRate / kappaRate0Tension) );
        double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );



        //compute yieldHard
        double yieldHardOne = computeHardeningOne(tempKappa);
        double yieldHardTwo = computeHardeningTwo(tempKappa);

        //compute elliptic function r
        double rFunction = ( 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + ( 2. * ecc - 1. ) * ( 2. * ecc - 1. ) ) /
            ( 2. * ( 1. - ecc * ecc ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta)
                                                           + 5. * ecc * ecc - 4. * ecc) );

        //compute help functions AL, BL
        double AL = ( 1. - yieldHardOne ) * pow( ( sig / fcYield + rho / ( sqrt(6.) * fcYield ) ), 2.) + sqrt(3. / 2.) * rho / fcYield;
        double BL = sig / fcYield + rho / ( fcYield * sqrt(6.) );

        //compute dfdsig
        double dfdsig = 4. * ( 1. - yieldHardOne ) / fcYield * AL * BL + yieldHardTwo * pow(yieldHardOne, 2.) * myield / fcYield;

        //compute dfdrho
        double dfdrho = AL / ( sqrt(6.) * fcYield ) * ( 4. * ( 1. - yieldHardOne ) * BL + 6. ) + rFunction * myield * yieldHardTwo * pow(yieldHardOne, 2.) / ( sqrt(6.) * fcYield );

        return {
            dfdsig, dfdrho
        };
    }

    FloatArrayF < 6 >
    ConcreteDPM2RatePlastic::computeDFDStress(const FloatArrayF < 6 > & stress, double tempKappa, GaussPoint * gp) const
    {
        double sig, rho, theta;
        computeCoordinates(stress, sig, rho, theta);
        auto dFDInv = computeDFDInv(sig, rho, theta, tempKappa, gp);

        double dRDCosTheta = computeDRDCosTheta(theta, this->ecc);

        double yieldHardOne = computeHardeningOne(tempKappa);
        double yieldHardTwo = computeHardeningTwo(tempKappa);

        //Compute dFDCosTheta. This was not needed for the stress return, but now for the tangent stiffness
        auto dFDCosTheta = dRDCosTheta * pow(yieldHardOne, 2.) * yieldHardTwo * m * rho / ( sqrt(6.) * fc );

        auto dSigDStress = computeDSigDStress();
        auto dRhoDStress = computeDRhoDStress(stress);
        auto dCosThetaDStress = computeDCosThetaDStress(stress);

        dSigDStress *= dFDInv.at(1);
        dRhoDStress *= dFDInv.at(2);
        dCosThetaDStress *= dFDCosTheta;

        auto dFDStress = dSigDStress + dRhoDStress + dCosThetaDStress;

        return dFDStress;
    }


    FloatMatrixF < 8, 8 >
    ConcreteDPM2RatePlastic::computeFullJacobian(const FloatArrayF < 6 > & stress,
                                                 const double deltaLambda,
                                                 GaussPoint * gp,
                                                 TimeStep * atTime,
                                                 const double tempKappa,
                                                 const double deltaTime) const
    {
        FloatMatrixF < 8, 8 > jacobian;
        auto dFDStress = computeDFDStress(stress, tempKappa, gp);
        auto dGDStress = computeDGDStress(stress, tempKappa);
        auto dDGDDStress = computeDDGDDStress(stress, tempKappa);
        auto dDGDStressDKappa = computeDDGDStressDKappa(stress, tempKappa);
        auto dDKappaDDeltaLambdaDStress = computeDDKappaDDeltaLambdaDStress(stress, tempKappa);

        double sig, rho, theta;
        computeCoordinates(stress, sig, rho, theta);

        double dFDKappa = computeDFDKappa(sig, rho, theta, tempKappa, deltaTime, gp);
        double dKappaDDeltaLambda  = computeDKappaDDeltaLambda(sig, rho, theta, tempKappa);

        double dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, theta, tempKappa);

        auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, atTime);

        auto helpA = dot(d, dDGDDStress);

        //Assign jacobian
        for ( int i = 0; i < 6; i++ ) {
            for ( int j = 0; j < 6; j++ ) {
                if ( i == j ) {
                    jacobian.at(i + 1, j + 1) = 1. + deltaLambda * helpA.at(i + 1, j + 1);
                } else {
                    jacobian.at(i + 1, j + 1) = deltaLambda * helpA.at(i + 1, j + 1);
                }
            }
        }

        FloatArrayF < 6 > helpB;
        helpB = dot(d, dDGDStressDKappa);

        for ( int i = 0; i < 6; i++ ) {
            jacobian.at(i + 1, 7) = deltaLambda * helpB.at(i + 1);
        }

        helpB = dot(d, dGDStress);

        for ( int i = 0; i < 6; i++ ) {
            jacobian.at(i + 1, 8) = helpB.at(i + 1);
        }

        for ( int i = 0; i < 6; i++ ) {
            jacobian.at(7, i + 1) = deltaLambda * dDKappaDDeltaLambdaDStress.at(i + 1);
        }


        jacobian.at(7, 7) = deltaLambda * dDKappaDDeltaLambdaDKappa - 1.;
        jacobian.at(7, 8) = dKappaDDeltaLambda;

        for ( int i = 0; i < 6; i++ ) {
            jacobian.at(8, i + 1) = dFDStress.at(i + 1);
        }

        jacobian.at(8, 7) = dFDKappa;
        jacobian.at(8, 8) = 0.;

        return jacobian;
    }


    FloatArrayF < 6 >
    ConcreteDPM2RatePlastic::giveRealStressVector_3d(const FloatArrayF < 6 > & fullStrainVector, GaussPoint * gp, TimeStep * tStep) const
    {
        auto status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );


        // Initialize temp variables for this gauss point
        status->initTempStatus();

        FloatArrayF < 6 > stress;

        if ( status->giveTempDeletionFlag() == 1 ) {
            //Deal with deleted element
            return stress;
        }

        //    Remove thermal/shrinkage strains
        auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);

        auto strainVector = fullStrainVector - thermalStrain;

        status->letTempReducedStrainBe(strainVector);

        //Calculate time increment
        double dt = deltaTime;
        if ( dt == -1 ) {
            if ( tStep->giveTimeIncrement() == 0 ) { //Problem with the first step. For some reason the time increment is zero
                dt = 1.;
            } else {
                dt = tStep->giveTimeIncrement();
            }
        }

        auto D = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, tStep);
        // perform plasticity return
        auto effectiveStress = performPlasticityReturn(gp, D, strainVector, dt);

        FloatArrayF < 6 > effectiveStressTension;
        FloatArrayF < 6 > effectiveStressCompression;


        double alpha = 0.;

        double tempKappa = status->giveTempKappaP();

        if ( this->damageFlag != 0 ) {//Apply damage
            alpha  = computeAlpha(effectiveStressTension, effectiveStressCompression, effectiveStress);
            auto damages = computeDamage(strainVector, D, dt, gp, tStep, alpha, effectiveStress);

            if ( this->damageFlag == 1 ) { //Default as described in IJSS CDPM2 article
                stress = effectiveStressTension * ( 1. - damages.at(1) ) + effectiveStressCompression * ( 1. - damages.at(2) );
            } else if ( this->damageFlag == 2 ) { //Simplified version without split of stress but two damage variables
                stress = effectiveStress * ( 1. - ( 1. - alpha ) * damages.at(1) ) * ( 1. - alpha * damages.at(2) );
            } else if ( this->damageFlag == 3 ) { //Consider only tensile damage. Reduction to a fully isotropic model. Similar to CDPM article.
                stress = effectiveStress * ( 1. - damages.at(1) );
            } else {
                OOFEM_ERROR("Unknown value of damage flag. Must be 0, 1, 2 or 3");
            }
        } else {
            stress = effectiveStress;
        }


        status->letTempStrainVectorBe(fullStrainVector);
        status->letTempAlphaBe(alpha);
        status->letTempStressVectorBe(stress);
        status->letTempEffectiveStressBe(effectiveStress);
#ifdef keep_track_of_dissipated_energy
        double gf = pow(ft, 2) / this->eM; //rough estimation only for this purpose
        status->computeWork(gp, gf);
#endif
        assignStateFlag(gp);
        return stress;
    }



    int
    ConcreteDPM2RatePlastic::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
    {
        ConcreteDPM2RatePlasticStatus *status = static_cast < ConcreteDPM2RatePlasticStatus * > ( this->giveStatus(gp) );
    }

    MaterialStatus *
    ConcreteDPM2RatePlastic::CreateStatus(GaussPoint *gp) const
    {
        return new ConcreteDPM2RatePlasticStatus(gp);
    }
}// end namespace oofem
