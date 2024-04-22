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

#include "concretedpm2Plasticrate1.h"
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
REGISTER_Material(ConcreteDPM2PlasticRate1);
  
  //****************************************
  //Status of ConcreteDPM2PlasticRate1
  //*************************************
  
ConcreteDPM2PlasticRate1Status::ConcreteDPM2PlasticRate1Status(GaussPoint *gp) : ConcreteDPM2Status(gp)
{
}


void
ConcreteDPM2PlasticRate1Status::initTempStatus()
{
    ConcreteDPM2Status::initTempStatus();
    this->tempRateTension = this->RateTension;
    this->tempRateCompression = this->RateCompression;
    this->tempBeta      = this->beta;
    this->tempRateFactorTension = this->rateFactorTension;
    this->tempRateFactorCompression = this->rateFactorCompression;
    //this->tempftYield = this->ftYield;
    //this->tempfcYield = this->fcYield;
}


void
ConcreteDPM2PlasticRate1Status::printOutputAt(FILE *file, TimeStep *tStep) const
{
    // Call corresponding function of the parent class to print
    ConcreteDPM2Status::printOutputAt(file, tStep);
    fprintf(file, " RateTension %.10e,", this->RateTension);
    fprintf(file, " rateFactorTension %.10e,", this->rateFactorTension);
    fprintf(file, " rateFactorCompression %.10e,", this->rateFactorCompression);
    fprintf(file, " RateCompression %.10e,", this->RateCompression);
}

void
ConcreteDPM2PlasticRate1Status::updateYourself(TimeStep *tStep)
{
    ConcreteDPM2Status::updateYourself(tStep);
    this->RateTension = this->tempRateTension;
    this->RateCompression = this->tempRateCompression;
    this->beta      = this->tempBeta;
    this->rateFactorTension = this->tempRateFactorTension;
    this->rateFactorCompression = this->tempRateFactorCompression;
    //this->tempftYield = this->ftYield;
    //this->tempfcYield = this->fcYield;
}


void
ConcreteDPM2PlasticRate1Status::saveContext(DataStream &stream, ContextMode mode)
{
    ConcreteDPM2Status::saveContext(stream, mode);

    //    contextIOResultType iores;
    contextIOResultType iores;

    if ( !stream.write(RateTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(RateCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactorTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactorCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }
/*
    if ( !stream.write(ftYield) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(fcYield) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    */
}

void
ConcreteDPM2PlasticRate1Status::restoreContext(DataStream &stream, ContextMode mode)
{
    ConcreteDPM2Status::restoreContext(stream, mode);
    contextIOResultType iores;

    if ( !stream.write(RateTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(RateCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateFactorTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateFactorCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }
/*
    if ( !stream.read(ftYield) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(fcYield) ) {
        THROW_CIOERR(CIO_IOERR);
    }*/
}
  
  //***************************
  //ConcreteDPM2Rate Class
  //***************************

#define IDM_ITERATION_LIMIT 1.e-8
ConcreteDPM2PlasticRate1::ConcreteDPM2PlasticRate1(int n, Domain *d) :
    ConcreteDPM2(n, d)
{}

  
void
ConcreteDPM2PlasticRate1::initializeFrom(InputRecord &ir)
{
/*
    ConcreteDPM2::initializeFrom(ir);
    this->cCompression = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->cCompression, _IFT_ConcreteDPM2PlasticRate1_cCompression);

    this->kappaRate0Compression = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, this->kappaRate0Compression, _IFT_ConcreteDPM2PlasticRate1_kappaRate0Compression);

    this->cTension = 0.01935;
    IR_GIVE_OPTIONAL_FIELD(ir, this->cTension, _IFT_ConcreteDPM2PlasticRate1_cTension);

    this->kappaRate0Tension = 0.00000434165;
    IR_GIVE_OPTIONAL_FIELD(ir, this->kappaRate0Tension, _IFT_ConcreteDPM2PlasticRate1_kappaRate0Tension);
    */
    ConcreteDPM2::initializeFrom(ir);
    this->atOne = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atOne, _IFT_ConcreteDPM2PlasticRate1_atOne);
    this->atTwo = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atTwo, _IFT_ConcreteDPM2PlasticRate1_atTwo);
    this->atThree = 0.018;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atThree, _IFT_ConcreteDPM2PlasticRate1_atThree);
    this->atFour = 0.0062;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atFour, _IFT_ConcreteDPM2PlasticRate1_atFour);
    this->atFive = 0.333;
    IR_GIVE_OPTIONAL_FIELD(ir, this->atFive, _IFT_ConcreteDPM2PlasticRate1_atFive);

    this->acOne = 30.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acOne, _IFT_ConcreteDPM2PlasticRate1_acOne);
    this->acTwo = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acTwo, _IFT_ConcreteDPM2PlasticRate1_acTwo);
    this->acThree = 0.014;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acThree, _IFT_ConcreteDPM2PlasticRate1_acThree);
    this->acFour = 0.012;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acFour, _IFT_ConcreteDPM2PlasticRate1_acFour);
    this->acFive = 0.333;
    IR_GIVE_OPTIONAL_FIELD(ir, this->acFive, _IFT_ConcreteDPM2PlasticRate1_acFive);

}

double
ConcreteDPM2PlasticRate1::computeYieldValue(double sig,
    double rho,
    double theta,
    double tempKappa,
    const double dt,
    GaussPoint *gp) const
{
    // compute yieldHard
    auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );


    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    //  compute elliptic function r
    double rFunction = ( 4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) + pow( ( 2. * ecc - 1. ), 2. ) ) / ( 2. * ( 1. - pow(ecc, 2.) ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) + 5. * pow(ecc, 2.) - 4. * ecc) );

    double ftYield    = computeFtYield( dt, gp);
    double fcYield    = computeFcYield( dt, gp);


    double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );

    // compute help function Al
    double Al = ( 1. - yieldHardOne ) * pow( ( sig / fcYield + rho / ( sqrt(6.) * fcYield ) ), 2. ) + sqrt(3. / 2.) * rho / fcYield;


    // Compute yield equation
    return pow(Al, 2.) + pow(yieldHardOne, 2.) * yieldHardTwo * myield * ( sig / fcYield + rho * rFunction / ( sqrt(6.) * fcYield ) ) - pow(yieldHardOne, 2.) * pow(yieldHardTwo, 2.);
}


double
ConcreteDPM2PlasticRate1::computeFcYield(const double dt, GaussPoint *gp) const
{

    auto status          = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );
    double tempRateFactorCompression = computeRateFactorCompression(gp, dt);

    status->setTempRateFactorCompression(tempRateFactorCompression);
    //double fcYield       = this->fc * ( 1 + this->cCompression * log(1. + tempRateCompression / this->kappaRate0Compression) );
    double fcYield       = this->fc * tempRateFactorCompression;
    return fcYield;
}

double
ConcreteDPM2PlasticRate1::computeRateFactorCompression(GaussPoint *gp,
    const double dt) const
{
    auto status          = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );
    double tempRateCompression = 0.;

    const auto &tempTotalstrain = status->giveTempReducedStrain();
    const auto &Totalstrain = status->giveReducedStrain();
    auto deltaTotalStrain = tempTotalstrain - Totalstrain;
    double deltaTotalStrainNorm = 0;
    deltaTotalStrainNorm= norm(deltaTotalStrain);

    if (deltaTotalStrainNorm == 0 ) {
        tempRateCompression = status->giveRateCompression();
    } else {
        tempRateCompression = deltaTotalStrainNorm / dt;
        //tempRateCompression = 50;
        printf("tempRateCompression = %e\n" , tempRateCompression);
    }
    tempRateCompression = status->giveRateCompression();
    status->setTempRateCompression(tempRateCompression);

    double strainRateValues= status->giveRateCompression();

    //For compression according to Model Code 2010
    double rateFactorCompression = 1.;
    double strainRateRatioCompression = strainRateValues / this->acOne ;
    if ( this->acOne < strainRateValues && strainRateValues < this->acTwo ) {
        rateFactorCompression = pow(strainRateRatioCompression, this->acThree);
    } else if ( this->acTwo <= strainRateValues ) {
        rateFactorCompression =  this->acFour * pow(strainRateRatioCompression, this->acFive);
    }

    return rateFactorCompression;
}






double
ConcreteDPM2PlasticRate1::computeFtYield( const double dt, GaussPoint *gp) const
{
    auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

    double tempRateFactorTension = computeRateFactorTension(gp, dt);

    status->setTempRateFactorTension(tempRateFactorTension);

    double ftYield = this->ft * tempRateFactorTension;
    return ftYield;
}

double
ConcreteDPM2PlasticRate1::computeRateFactorTension(GaussPoint *gp,
    const double dt) const
{
    auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

    // Rate factor in tension has to be made mesh independent once damage has started, because the model is based on the crack band approach.
    // It is assumed that damage starts once tempKappaP is greater than 1.
    double tempBeta                 = 0.;
    double tempRateTension     = 0.;
    double tempRateTensionTest = 0.;

    const auto &tempTotalstrain = status->giveTempReducedStrain();
    const auto &Totalstrain = status->giveReducedStrain();
    auto deltaTotalStrain = tempTotalstrain - Totalstrain;
    double deltaTotalStrainNorm = 0;
    deltaTotalStrainNorm= norm(deltaTotalStrain);

/*
    double le = status->giveLe();
    if ( le == 0 ) { // So that it works for only plasticity
        le = gp->giveElement()->computeMeanSize();
    }
*/
    if ( deltaTotalStrainNorm == 0) {
        tempRateTension = status->giveRateTension();
   } else {
        tempRateTension = deltaTotalStrainNorm / dt;
        printf("tempRateTension = %e\n" , tempRateTension);
        //tempRateTension = 50;

      /*
        if ( status->giveKappaP() <= 1.) {
            // Damage is zero
            tempRateTension = deltaTotalStrainNorm / dt;
        } else {
            // Damage in previous step is not zero
            tempBeta = status->giveBeta();
            if ( tempBeta == 0 ) {
                // Calculate tempBeta only once
                tempBeta = status->giveRateTension() / ( le * ( deltaTotalStrainNorm ) / dt );
            }
            tempRateTension = tempBeta * le * ( deltaTotalStrainNorm ) / dt;
        }
        */
   }
     status->setTempRateTension(tempRateTension);
    //status->setTempBeta(tempBeta);
    double strainRateValues= status->giveRateTension();

    //double ftYield = this->ft * ( 1 + this->cTension * log(1. + tempRateTension / this->kappaRate0Tension) );

    double rateFactorTension = 1.;
    double strainRateRatioTension = strainRateValues / this->atOne;

    if (strainRateValues >this->atOne &&  strainRateValues < this->atTwo ) {
        rateFactorTension = pow(strainRateRatioTension, this->atThree);
    }else if (strainRateValues >= this->atTwo){
        rateFactorTension =  this->atFour * pow(strainRateRatioTension, this->atFive);
    }
    //For tension according to Model Code 2010


    return rateFactorTension;
}



FloatArrayF < 6 >
ConcreteDPM2PlasticRate1::performPlasticityReturn(GaussPoint * gp, const FloatMatrixF < 6, 6 > & D, const FloatArrayF < 6 > & strain, const double dt) const
{
    auto status = static_cast < ConcreteDPM2Status * > ( this->giveStatus(gp) );

    ConcreteDPM2_ReturnResult returnResult = RR_Unknown;
    ConcreteDPM2_ReturnType returnType = RT_Unknown;

    //get plastic strain and kappa
    auto tempPlasticStrain = status->givePlasticStrain();
    double tempKappaP = status->giveKappaP();

    //this theta computed here should stay constant for the rest of procedure.`

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
        double yieldValue = computeYieldValue(sig, rho, theta, tempKappaP, dt, gp);

        apexStress = 0.;

        if ( yieldValue > 0. ) {
            checkForVertexCase(apexStress, returnType, sig, tempKappaP, gp);
            if ( returnType == RT_Tension || returnType == RT_Compression ) {
                tempKappaP = performVertexReturn(effectiveStress, returnResult, returnType, apexStress, tempKappaP, gp, dt);
                status->letTempKappaPBe(tempKappaP);
                if ( returnType == RT_Tension ) {
                    status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_VertexTension);
                } else if ( returnType == RT_Compression ) {
                    status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_VertexCompression);
                }
            }
            if ( returnType == RT_Regular ) {
                printf("in return \n ");
                tempKappaP = performRegularReturn(effectiveStress, returnResult, returnType, tempKappaP, gp, theta, dt);

                status->letTempKappaPBe(tempKappaP);
            }
        } else {
            printf("no yielding \n ");
            returnResult = RR_Converged;
            tempPlasticStrain = status->givePlasticStrain();
            status->letTempPlasticStrainBe(tempPlasticStrain);
            status->letTempKappaPBe(tempKappaP);
            break;
        }

        if ( returnResult == RR_NotConverged ) {
            subincrementcounter++;
            if ( subincrementcounter > 10 ) {
                OOFEM_LOG_INFO( "Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                OOFEM_LOG_INFO( "Old strain vector %g %g %g %g %g %g  \n", oldStrain.at(1), oldStrain.at(2), oldStrain.at(3), oldStrain.at(4), oldStrain.at(5), oldStrain.at(6) );

                const auto &help = status->giveTempPlasticStrain();
                OOFEM_LOG_INFO( "Old plastic strain vector %g %g %g %g %g %g  \n", help.at(1), help.at(2), help.at(3), help.at(4), help.at(5), help.at(6) );
                OOFEM_LOG_INFO( "New strain vector %g %g %g %g %g %g  \n", strain.at(1), strain.at(2), strain.at(3), strain.at(4), strain.at(5), strain.at(6) );

                computeCoordinates(effectiveStress, sig, rho, theta);
                double sig1, rho1, theta1;
                auto help1 = dot(D, oldStrain - help);
                computeCoordinates(help1, sig1, rho1, theta1);
                yieldValue = computeYieldValue(sig, rho, theta, tempKappaP, dt, gp);
                OOFEM_LOG_INFO("OLD Sig %g rho %g theta %g  \n", sig1, rho1, theta1);
                OOFEM_LOG_INFO("NEW Sig %g rho %g theta %g  \n", sig, rho, theta);
                if ( returnType == RT_Tension || returnType == RT_Compression ) {
                    OOFEM_LOG_INFO("Vertex case apexstress %g\n", apexStress);
                } else {
                    OOFEM_LOG_INFO("Regular case %g \n", 15.18);
                }
                OOFEM_LOG_INFO("KappaP old %g new %g yieldfun %g\n", status->giveTempKappaP(), tempKappaP, yieldValue);
                OOFEM_WARNING( "ConcreteDamagePlasticity2:: performPlasticityReturn: Could not reach convergence with small deltaStrain, giving up. Delete Element number %d", gp->giveElement()->giveNumber() );
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
ConcreteDPM2PlasticRate1::performVertexReturn(FloatArrayF < 6 > &effectiveStress,
    ConcreteDPM2_ReturnResult &returnResult,
    ConcreteDPM2_ReturnType &returnType,
    double apexStress, double tempKappaP,
    GaussPoint *gp, const double dt) const
{
    // auto [deviatoricStressTrial, sigTrial] = computeDeviatoricVolumetricSplit(effectiveStress); // c++17
    auto tmp                   = computeDeviatoricVolumetricSplit(effectiveStress);
    auto deviatoricStressTrial = tmp.first;
    auto sigTrial              = tmp.second;

    auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

    double rhoTrial = computeSecondCoordinate(deviatoricStressTrial);

    double kappaInitial = tempKappaP;

    double sig2 = apexStress;

    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial,gp,dt);

    double yieldValue = computeYieldValue(sigTrial, 0., 0., tempKappaP, dt, gp);

    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2,gp,dt);

    double yieldValueMid = computeYieldValue(sig2, 0., 0., tempKappaP, dt, gp);

    if ( yieldValue * yieldValueMid >= 0. ) {
        returnType   = RT_Regular;
        returnResult = RR_NotConverged;
        return kappaInitial;
    }

    double dSig, sigAnswer;
    if ( yieldValue < 0.0 ) {
        dSig      = sig2 - sigTrial;
        sigAnswer = sig2;
    } else {
        dSig      = sigTrial - sig2;
        sigAnswer = sig2;
    }

    for ( int j = 0; j < 250; j++ ) {
        dSig = 0.5 * dSig;

        double sigMid = sigAnswer + dSig;


        tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigMid,gp,dt);

        yieldValueMid = computeYieldValue(sigMid, 0., 0., tempKappaP, dt, gp);

        if ( yieldValueMid <= 0. ) {
            sigAnswer = sigMid;
        }

        if ( fabs(yieldValueMid) < yieldTol && yieldValueMid <= 0. ) {
            double ratioPotential = computeRatioPotential(sigAnswer, 0, tempKappaP,gp,dt);


            double ratioTrial = rhoTrial / ( sigTrial - sigAnswer );

            if ( ( ( ( ratioPotential >= ratioTrial ) && returnType == RT_Tension ) ) || ( ( ratioPotential <= ratioTrial ) && returnType == RT_Compression ) ) {
                for ( int i = 0; i < 3; i++ ) {
                    effectiveStress.at(i + 1) = sigAnswer;
                }

                for ( int i = 3; i < 6; i++ ) {
                    effectiveStress.at(i + 1) = 0.;
                }
                returnResult = RR_Converged;
                return tempKappaP;
            } else {
                returnType   = RT_Regular;
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
ConcreteDPM2PlasticRate1::performRegularReturn(FloatArrayF < 6 > &effectiveStress,
    ConcreteDPM2_ReturnResult &returnResult,
    ConcreteDPM2_ReturnType &returnType,
    double kappaP,
    GaussPoint *gp,
    double theta,
    const double dt) const
{
    auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

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
    double yieldValue = computeYieldValue(sig, rho, theta, tempKappaP, dt, gp);

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
            auto jacobian = computeJacobian(sig, rho, theta, tempKappaP, deltaLambda, dt,gp);

            try {
                auto deltaIncrement = solve(jacobian, FloatArrayF < 4 > ( residuals ) );
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
            auto dGDInv = computeDGDInv(sig, rho, tempKappaP,gp,dt);
            double dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, theta, tempKappaP,gp, dt);

            residuals.at(1) = sig - trialSig + this->kM * deltaLambda * dGDInv.at(1);
            residuals.at(2) = rho - trialRho + ( 2. * this->gM ) * deltaLambda * dGDInv.at(2);
            residuals.at(3) = -tempKappaP + kappaP + deltaLambda * dKappaDDeltaLambda;
            residuals.at(4) = computeYieldValue(sig, rho, theta, tempKappaP, dt, gp);
        }
    }


    //compute the principal directions of the stress
    //auto [helpStress, stressPrincipalDir] = StructuralMaterial :: computePrincipalValDir(from_voigt_stress(trialStress)); // c++17
    auto tmpEig = StructuralMaterial::computePrincipalValDir(from_voigt_stress(trialStress) );
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


 FloatArrayF < 6 >
 ConcreteDPM2PlasticRate1::giveRealStressVector_3d(const FloatArrayF < 6 > & fullStrainVector, GaussPoint * gp, TimeStep * tStep) const
 {
   auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );
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
                        double dt = this->deltaTime;
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

                           double ftYield = computeFtYield(dt, gp);
                   #ifdef keep_track_of_dissipated_energy
                           double gf = pow(ftYield, 2) / this->eM; //rough estimation only for this purpose
                           status->computeWork(gp, gf);
                   #endif
                           assignStateFlag(gp);
                           return stress;
}

FloatArrayF < 2 >
ConcreteDPM2PlasticRate1::computeDFDInv(double sig, double rho, double theta, double tempKappa, const double dt, GaussPoint * gp) const
{
                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

                           double fcYield = computeFcYield(dt, gp);
                           double ftYield = computeFtYield(dt, gp);
                           //double fcYield = fc;
                           //double ftYield = ft;
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

double
ConcreteDPM2PlasticRate1::computeEquivalentStrainP(double sig, double rho, double theta, const double dt, GaussPoint *gp) const
{
                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

                           double tempKappaP = status->giveTempKappaP();

                           double ftYield    = computeFtYield( dt, gp);
                           double fcYield    = computeFcYield( dt, gp);
                           //double ftYield    = ft;
                          // double fcYield    = fc;

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
ConcreteDPM2PlasticRate1::computeDFDKappa(double sig,
                                          double rho,
                                          double theta,
                                          double tempKappa,
                                          const double dt,
                                          GaussPoint *gp) const
{
                           double dFDKappa;

                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );


                           //compute yieldHard and yieldSoft
                           double yieldHardOne = computeHardeningOne(tempKappa);
                           double yieldHardTwo = computeHardeningTwo(tempKappa);
                           // compute the derivative of the hardening and softening laws
                           double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
                           double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);

                           double fcYield = computeFcYield( dt, gp);
                           double ftYield = computeFtYield(dt, gp);
                           //double fcYield = fc;
                           //double ftYield = ft;
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

                           // compute dFDKappa
                           dFDKappa = dFDYieldHardOne * dYieldHardOneDKappa + dFDYieldHardTwo * dYieldHardTwoDKappa;
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

int
ConcreteDPM2PlasticRate1::checkForUnAndReloadingP(double &tempEquivStrain,
                                                 double &minEquivStrain,
                                                 const FloatMatrixF < 6, 6 > &D,
                                                 const double dt,
                                                 GaussPoint *gp) const
{
                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

                           //Access old and new strains
                           const auto &oldStrain = status->giveReducedStrain();
                           const auto &strain = status->giveTempReducedStrain();

                           //Compute the temp equivalent strain
                           auto tempElasticStrain = strain - status->giveTempPlasticStrain();
                           auto tempEffectiveStress = dot(D, tempElasticStrain);

                           double sigEffective, rhoEffective, thetaEffective;
                           computeCoordinates(tempEffectiveStress, sigEffective, rhoEffective, thetaEffective);
                           tempEquivStrain = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective,dt,gp);
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
                           double equivStrainPlus = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective,dt,gp);

                           //For slightly less than temp effective stress
                           auto intermediateEffectiveStressMinus = effectiveStress + 0.99 * deltaEffectiveStress;
                           computeCoordinates(intermediateEffectiveStressMinus, sigEffective, rhoEffective, thetaEffective);
                           double tempEquivStrainMinus = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective,dt,gp);

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
                                    double midEquivStrain = computeEquivalentStrainP(sigEffective, rhoEffective, thetaEffective,dt,gp);

                                    if ( midEquivStrain <= minEquivStrain ) {
                minEquivStrain = midEquivStrain;
                                    } else {
                return unloadingFlag;
                                    }
                                }
                           }
                           return unloadingFlag;
}

FloatArrayF < 2 >
ConcreteDPM2PlasticRate1::computeDamage(const FloatArrayF < 6 > & strain,
    const FloatMatrixF < 6, 6 > & D,
    double dt,
    GaussPoint * gp,
    TimeStep * tStep,
    double tempAlpha,
    const FloatArrayF < 6 > & effectiveStress) const
{
                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

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

                           int unAndReloadingFlag = checkForUnAndReloadingP(tempEquivStrain, minEquivStrain, D, dt,gp);

                           double rateFactor;
                           if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) {
                                rateFactor = computeRateFactor(tempAlpha, deltaTime, gp, tStep);
                           } else {
                                rateFactor = status->giveRateFactor();
                           }


                           //Compute equivalent strains for  tension and compression
                           double tempEquivStrainTension = 0.;
                           double tempEquivStrainCompression = 0.;

                           tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;

                           if ( unAndReloadingFlag == 0 ) { //Standard way
                                tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
                           } else {
                                tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
                           }


                           //If damage threshold is exceeded determine the rate factor from the previous step
                           if ( ( tempEquivStrainTension > e0 || tempEquivStrainCompression > e0 ) &&
                               ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) && !tStep->isTheFirstStep() ) {
                                //Rate factor from last step
                                rateFactor = status->giveRateFactor();

                                tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;
                                if ( unAndReloadingFlag == 0 ) { //Standard way
                                    tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
                                } else {
                                    tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
                                }
                           }

                           status->letTempRateFactorBe(rateFactor);

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
                           } else if ( fTension >= -yieldTolDamage && fCompression < -yieldTolDamage ) { //Only tension is active
                                //Update tension history variables
                                tempKappaDTension = tempEquivStrainTension;
                                deltaPlasticStrainNorm = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
                                tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNorm / ductilityMeasure / rateFactor;
                                tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

                                //Nothing changes for compression history variables
                                tempKappaDCompression = status->giveKappaDCompression();
                                tempKappaDCompressionOne = status->giveKappaDCompressionOne();
                                tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

                                //Initialise damage with tensile history variable
                                this->initDamaged(tempKappaDTension, strain, gp);

                                tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), rateFactor);

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
                                tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
                                tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() + ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

                                //Determine damage parameters
                                tempDamageTension = status->giveDamageTension();
                                tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression(), rateFactor);
                           } else if ( fTension >= -yieldTolDamage && fCompression >= -yieldTolDamage ) {
                                //Both tension and compression is active

                                //Update tension history variables
                                tempKappaDTension = tempEquivStrainTension;
                                deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
                                tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ( ductilityMeasure * rateFactor );
                                tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

                                //Update the compression history variables
                                tempKappaDCompression = tempEquivStrainCompression;
                                deltaPlasticStrainNormCompression =
                                    computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
                                tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
                                tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() +
                                    ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

                                //Determine the damage parameters
                                this->initDamaged(tempKappaDTension, strain, gp);

                                tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), rateFactor);

                                tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression(), rateFactor);
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

FloatMatrixF < 4, 4 >
ConcreteDPM2PlasticRate1::computeJacobian(double sig,
    double rho,
    double theta,
    double kappa,
    double deltaLambda,
    const double dt,
    GaussPoint * gp) const
{
                           auto dFDInv = computeDFDInv(sig, rho, theta, kappa,dt,gp);
                           auto dGDInv = computeDGDInv(sig, rho, kappa,gp,dt);
                           auto dDGDDInv = computeDDGDDInv(sig, rho, kappa,gp,dt);

                           double dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, theta, kappa, gp,dt);
                           double dFDKappa = computeDFDKappa(sig, rho, theta, kappa,dt,gp);

                           auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, kappa,gp,dt);

                           double dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, theta, kappa, gp,dt);
                           auto dDKappaDDeltaLambdaDInv = computeDDKappaDDeltaLambdaDInv(sig, rho, theta, kappa,gp,dt);

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

FloatArrayF < 6 >
ConcreteDPM2PlasticRate1::computeDFDStress(const FloatArrayF < 6 > & stress,
    double tempKappa, const double dt, GaussPoint * gp) const
{
                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

                           double ftYield    = computeFtYield( dt, gp);
                           double fcYield    = computeFcYield( dt, gp);
                           //double ftYield    = ft;
                           //double fcYield    = fc;

                           double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );

                           double sig, rho, theta;
                           computeCoordinates(stress, sig, rho, theta);
                           auto dFDInv = computeDFDInv(sig, rho, theta, tempKappa,dt,gp);

                           double dRDCosTheta = computeDRDCosTheta(theta, this->ecc);

                           double yieldHardOne = computeHardeningOne(tempKappa);
                           double yieldHardTwo = computeHardeningTwo(tempKappa);

                           //Compute dFDCosTheta. This was not needed for the stress return, but now for the tangent stiffness
                           auto dFDCosTheta = dRDCosTheta * pow(yieldHardOne, 2.) * yieldHardTwo * myield * rho / ( sqrt(6.) * fcYield );

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
ConcreteDPM2PlasticRate1::computeFullJacobian(const FloatArrayF < 6 > & stress,
    const double deltaLambda,
    GaussPoint * gp,
    TimeStep * atTime,
    const double tempKappa, const double dt) const
{
                           FloatMatrixF < 8, 8 > jacobian;
                           auto dFDStress = computeDFDStress(stress, tempKappa,dt,gp);
                           auto dGDStress = computeDGDStress(stress, tempKappa,gp,dt);
                           auto dDGDDStress = computeDDGDDStress(stress, tempKappa,gp,dt);
                           auto dDGDStressDKappa = computeDDGDStressDKappa(stress, tempKappa,gp,dt);
                           auto dDKappaDDeltaLambdaDStress = computeDDKappaDDeltaLambdaDStress(stress, tempKappa,gp,dt);

                           double sig, rho, theta;
                           computeCoordinates(stress, sig, rho, theta);

                           double dFDKappa = computeDFDKappa(sig, rho, theta, tempKappa,dt,gp);
                           double dKappaDDeltaLambda  = computeDKappaDDeltaLambda(sig, rho, theta, tempKappa, gp,dt);

                           double dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, theta, tempKappa,gp ,dt);

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

FloatMatrixF < 6, 6 >
ConcreteDPM2PlasticRate1::compute3dTangentStiffness(GaussPoint * gp, TimeStep * tStep, const double dt) const
{
                           FloatMatrixF < 6, 6 > answer;

                           FloatArrayF < 6 > effectiveStress;
                           auto status = static_cast < ConcreteDPM2Status * > ( giveStatus(gp) );
                           effectiveStress = status->giveTempEffectiveStress();

                           double tempKappa = status->giveTempKappaP();
                           double deltaLambda = status->giveDeltaLambda();

                           //Computes only the plastic part of the tangent stiffness
                           auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, tStep);

                           // Debug
                           FloatMatrixF < 8, 8 > fullJacobian = computeFullJacobian(effectiveStress, deltaLambda, gp, tStep, tempKappa,dt);

                           FloatMatrixF < 8, 8 > invFullJacobian = inv(fullJacobian);

                           FloatMatrixF < 6, 6 > help;
                           for ( int i = 1; i <= 6; i++ ) {
                                for ( int j = 1; j <= 6; j++ ) {
                                    help.at(i, j) = invFullJacobian.at(i, j);
                                }
                           }

                           answer = dot(help, d);

                           //  Damage parameters
                           double omegaTension = min(status->giveTempDamageTension(), 0.999999);
                           double omegaCompression = min(status->giveTempDamageCompression(), 0.999999);
                           double alpha = status->giveTempAlpha();

                           if ( damageFlag == 2 ) {
                                answer *= ( 1. - ( 1. - alpha ) * omegaTension ) * ( 1. - alpha * omegaCompression );
                           } else if ( damageFlag == 3 ) {
                                answer *= ( 1. - omegaTension );
                           }

                           return answer;
}



double
ConcreteDPM2PlasticRate1::computeDuctilityMeasure(double sig,
    double rho,
    double theta,
    GaussPoint *gp,
    const double dt) const
{
                           auto status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );

                           double fcYield = computeFcYield( dt, gp);
                           //double fcYield = fc;

                           double thetaConst = pow(2. * cos(theta), 2.);
                           double ductilityMeasure;
                           double x = -( sig + fcYield / 3 ) / fcYield;
                           if ( x < 0. ) {
                                /*Introduce exponential help function which results in a smooth
             * transition. */
                                double EHard = BHard - DHard;
                                double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );
                                ductilityMeasure = ( EHard * exp(x / FHard) + DHard ) / thetaConst;
                           } else {
                                ductilityMeasure = ( AHard + ( BHard - AHard ) * exp( -x / ( CHard ) ) ) / thetaConst;
                           }

                           return ductilityMeasure;
}

double
ConcreteDPM2PlasticRate1::computeTempKappa(double kappaInitial,
    double sigTrial,
    double rhoTrial,
    double sig,
    GaussPoint *gp,
    const double dt) const
{
                           //This function is called, if stress state is in vertex case
                           double equivalentDeltaPlasticStrain = sqrt( 1. / 9. * pow( ( sigTrial - sig ) / ( kM ), 2.) +
                               pow(rhoTrial / ( 2. * gM ), 2.) );

                           double thetaVertex = M_PI / 3.;
                           double ductilityMeasure = computeDuctilityMeasure(sig, 0., thetaVertex,gp,dt);

                           return kappaInitial + equivalentDeltaPlasticStrain / ductilityMeasure;
}


double
ConcreteDPM2PlasticRate1::computeDKappaDDeltaLambda(double sig,
    double rho,
    double theta,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{
                           auto dGDInv = computeDGDInv(sig, rho, tempKappa,gp,dt);
                           double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );
                           double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta,gp,dt);
                           return equivalentDGDStress / ductilityMeasure; // dKappaDDeltaLambda
}

FloatArrayF < 2 >
ConcreteDPM2PlasticRate1::computeDDKappaDDeltaLambdaDInv(double sig,
    double rho,
    double theta,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{
                           //Compute first and second derivative of plastic potential
                           auto dGDInv = computeDGDInv(sig, rho, tempKappa,gp,dt);
                           auto dDGDDInv = computeDDGDDInv(sig, rho, tempKappa,gp,dt);

                           //Compute equivalentDGDStress
                           double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );

                           //computeDuctilityMeasure
                           double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta, gp,dt);

                           //Compute dEquivalentDGDStressDInv
                           FloatArrayF < 2 > dEquivalentDGDStressDInv;
                           dEquivalentDGDStressDInv [ 0 ] =
                               ( 2. / 3. * dGDInv [ 0 ] * dDGDDInv(0, 0) + 2. * dGDInv [ 1 ] * dDGDDInv(1, 0) ) / ( 2. * equivalentDGDStress );
                           dEquivalentDGDStressDInv [ 1 ] =
                               ( 2. / 3. * dGDInv [ 0 ] * dDGDDInv(0, 1) + 2. * dGDInv [ 1 ] * dDGDDInv(1, 1) ) / ( 2. * equivalentDGDStress );


                           // compute the derivative of
                           auto dDuctilityMeasureDInv = computeDDuctilityMeasureDInv(sig, rho, theta, tempKappa,gp,dt);

                           FloatArrayF < 2 > answer;
                           answer [ 0 ] = ( dEquivalentDGDStressDInv [ 0 ] * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv [ 0 ] ) / pow(ductilityMeasure, 2.);
                           answer [ 1 ] = ( dEquivalentDGDStressDInv [ 1 ] * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv [ 1 ] ) / pow(ductilityMeasure, 2.);
                           return answer;
}

FloatArrayF < 6 >
ConcreteDPM2PlasticRate1::computeDDKappaDDeltaLambdaDStress(const FloatArrayF < 6 > & stress, double tempKappa,GaussPoint *gp, const double dt) const
{
                           auto tmp = computeDeviatoricVolumetricSplit(stress);
                           auto deviatoricStress = tmp.first;
                           double sig = tmp.second;

                           double rho = computeSecondCoordinate(deviatoricStress);
                           double theta = computeThirdCoordinate(deviatoricStress);

                           auto dDKappaDDeltaLambdaDInv = computeDDKappaDDeltaLambdaDInv(sig, rho, theta, tempKappa,gp,dt);

                           // compute dDKappaDDeltaLambdaDCosTheta
                           auto dGDInv = computeDGDInv(sig, rho, tempKappa,gp,dt);

                           double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );

                           double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta,gp,dt);

                           //Reuse implementation to compute dKappaDDeltaLambdaDCosTheta
                           //Ductility measure has in denominator (2*cos(theta))^2
                           double dDKappaDDeltaLambdaDCosTheta = equivalentDGDStress / ductilityMeasure / cos(theta);
                           //  dDKappaDDeltaLambdaDCosTheta = 0.;

                           auto dSigDStress = computeDSigDStress();
                           auto dRhoDStress = computeDRhoDStress(stress);
                           auto dCosThetaDStress = computeDCosThetaDStress(stress);

                           dSigDStress *= dDKappaDDeltaLambdaDInv.at(1);

                           dRhoDStress *= dDKappaDDeltaLambdaDInv.at(2);

                           dCosThetaDStress *= dDKappaDDeltaLambdaDCosTheta;

                           auto dDKappaDDeltaLambdaDStress = dSigDStress + dRhoDStress + dCosThetaDStress;

                           return dDKappaDDeltaLambdaDStress;
}

double
ConcreteDPM2PlasticRate1::computeDDKappaDDeltaLambdaDKappa(double sig, double rho, double theta, double tempKappa,GaussPoint *gp, const double dt) const
{
                           //Compute first and second derivative of plastic potential
                           auto dGDInv = computeDGDInv(sig, rho, tempKappa,gp,dt);
                           auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, tempKappa,gp,dt);

                           double equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );

                           double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta,gp,dt);
                           //Compute dEquivalentDGDStressDKappa
                           double dEquivalentDGDStressDKappa =
                               ( 2. / 3. * dGDInv [ 0 ] * dDGDInvDKappa [ 0 ] + 2. * dGDInv [ 1 ] * dDGDInvDKappa [ 1 ] ) / ( 2. * equivalentDGDStress );

                           return dEquivalentDGDStressDKappa / ductilityMeasure;
}


double
ConcreteDPM2PlasticRate1::computeRatioPotential(double sig,
    double rho,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{
                           //compute yieldHard and yieldSoft
                           double yieldHardOne = computeHardeningOne(tempKappa);
                           double yieldHardTwo = computeHardeningTwo(tempKappa);

                           double ftYield    = computeFtYield( dt, gp);
                           double fcYield    = computeFcYield( dt, gp);

                           //double ftYield    = ft;
                           //double fcYield    = fc;

                           double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );

                           //Compute dilation parameter
                           double AGParam = this->ft * yieldHardTwo * 3 / this->fc + myield / 2;
                           double BGParam =
                               yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
                               ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

                           double R = ( sig - ftYield / 3. * yieldHardTwo ) / fcYield / BGParam;
                           double mQ = AGParam * exp(R);

                           double Bl = sig / fcYield + rho / ( fcYield * sqrt(6.) );

                           double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fcYield;

                           double dgdsig = 4. * ( 1. - yieldHardOne ) / fcYield * Al * Bl + pow(yieldHardOne, 2.) * mQ / fcYield;

                           double dgdrho = Al / ( sqrt(6.) * fcYield ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +

                               myield * pow(yieldHardOne, 2.) / ( sqrt(6.) * fcYield );

                           return dgdrho / dgdsig * 3. * ( 1. - 2. * nu ) / ( 1. + nu );
}

FloatArrayF < 2 >
ConcreteDPM2PlasticRate1::computeDDuctilityMeasureDInv(double sig,
    double rho,
    double theta,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{

                          double fcYield    = computeFcYield( dt, gp);

                           //double fcYield = fc;

                           double thetaConst = pow(2. * cos(theta), 2.);
                           double x = ( -( sig + fcYield  / 3. ) ) / fcYield ;

                           if ( x < 0. ) {
                                double dXDSig = -1. / fcYield ;
                                /* Introduce exponential help function which results in a
             * smooth transition. */
                                double EHard = BHard - DHard;
                                double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );

                                double dDuctilityMeasureDX = EHard / FHard * exp(x / FHard) / thetaConst;
                                return {
                                    dDuctilityMeasureDX *dXDSig, 0.
                                };
                           } else {
                                double dXDSig = -1. / fcYield ;
                                double dDuctilityMeasureDX = -( BHard - AHard ) / ( CHard ) / thetaConst * exp( -x / ( CHard ) );
                                return {
                                    dDuctilityMeasureDX *dXDSig, 0.
                                };
                           }
}


FloatArrayF < 2 >
ConcreteDPM2PlasticRate1::computeDGDInv(double sig,
    double rho,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{
                           //compute yieldHard and yieldSoft
                           double yieldHardOne = computeHardeningOne(tempKappa);
                           double yieldHardTwo = computeHardeningTwo(tempKappa);

                           double ftYield    = computeFtYield( dt, gp);
                           double fcYield    = computeFcYield( dt, gp);

                           //double ftYield    = ft;
                           //double fcYield    = fc;

                           double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );
                           //Compute dilation parameter
                           double AGParam = this->ft * yieldHardTwo * 3 / this->fc + myield / 2;
                           double BGParam =
                               yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
                               ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

                           double R = ( sig - ftYield / 3. * yieldHardTwo ) / fcYield / BGParam;
                           double mQ = AGParam * exp(R);

                           double Bl = sig / fcYield + rho / ( fcYield * sqrt(6.) );

                           double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fcYield;

                           double dgdsig = 4. * ( 1. - yieldHardOne ) / fcYield * Al * Bl + pow(yieldHardOne, 2.) * mQ / fcYield;

                           double dgdrho = Al / ( sqrt(6.) * fcYield ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                               myield * pow(yieldHardOne, 2.) / ( sqrt(6.) * fcYield );

                           return {
                               dgdsig, dgdrho
                           };
}

FloatArrayF < 6 >
ConcreteDPM2PlasticRate1::computeDGDStress(const FloatArrayF < 6 > & stress, const double tempKappa,GaussPoint *gp,const double dt) const
{
                           auto tmp = computeDeviatoricVolumetricSplit(stress);
                           auto deviatoricStress = tmp.first;
                           double sig = tmp.second;
                           double rho = computeSecondCoordinate(deviatoricStress);

                           //compute dGDSig*dSigDStress + dGDRho*dRhoDStress
                           auto dGDInv = computeDGDInv(sig, rho, tempKappa,gp,dt);
                           auto dSigDStress = computeDSigDStress();
                           auto dRhoDStress = computeDRhoDStress(stress);

                           dSigDStress *= dGDInv.at(1);
                           dRhoDStress *= dGDInv.at(2);
                           dSigDStress += dRhoDStress;

                           return dSigDStress;
}

FloatMatrixF < 6, 6 >
ConcreteDPM2PlasticRate1::computeDDGDDStress(const FloatArrayF < 6 > & stress, const double tempKappa,GaussPoint *gp,const double dt) const
{
                           FloatMatrixF < 6, 6 > answer;

                           auto tmp = computeDeviatoricVolumetricSplit(stress);
                           auto deviatoricStress = tmp.first;
                           double sig = tmp.second;
                           double rho = computeSecondCoordinate(deviatoricStress);

                           //compute deriviates with respect to invariants
                           auto dDGDDInv = computeDDGDDInv(sig, rho, tempKappa,gp,dt);
                           auto dGDInv = computeDGDInv(sig, rho, tempKappa,gp,dt);

                           //Compute derivatives of invariants with respect to stress
                           auto dRhoDStress =  computeDRhoDStress(stress);

                           auto dDRhoDDStress = computeDDRhoDDStress(stress);

                           auto dSigDStress =  computeDSigDStress();

                           //Assemble terms

                           //Compute (dDGDDSig*dSigDStress + dDGDSigDRho*dRhoDStress)*dSigDStress
                           FloatArrayF < 6 > temp1;
                           temp1 = dSigDStress;
                           temp1 *= dDGDDInv.at(1, 1);

                           FloatArrayF < 6 > temp2;
                           temp2 = dRhoDStress;
                           temp2 *= dDGDDInv.at(1, 2);

                           temp1 += temp2;

                           FloatMatrixF < 6, 6 > helpA;
                           helpA = dyad(temp1, dSigDStress);

                           //Compute (dDGDDRho*dRhoDStress + dDGDRhoDSig*dSigDStress)*dRhoDstress
                           temp1 = dRhoDStress;
                           temp1 *= dDGDDInv.at(2, 2);

                           temp2 = dSigDStress;
                           temp2 *= dDGDDInv.at(2, 1);

                           temp1 += temp2;

                           FloatMatrixF < 6, 6 > helpB;
                           helpB = dyad(temp1, dRhoDStress);

                           //compute dGDRho * dDRhoDDStress
                           FloatMatrixF < 6, 6 > helpC = dDRhoDDStress;
                           helpC *= dGDInv.at(2);

                           //The term dGDSig*dDSigDDStress is zero

                           //sum up all parts
                           answer = helpA;
                           answer += helpB;
                           answer += helpC;

                           return answer;
}

FloatArrayF < 2 >
ConcreteDPM2PlasticRate1::computeDDGDInvDKappa(double sig,
    double rho,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{
                           //Compute dilation parameter

                           //compute yieldHard and yieldSoft
                           double yieldHardOne = computeHardeningOne(tempKappa);
                           double yieldHardTwo = computeHardeningTwo(tempKappa);

                           double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
                           double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);

                           double ftYield    = computeFtYield( dt, gp);
                           double fcYield    = computeFcYield( dt, gp);

                           //double ftYield    = ft;
                           //double fcYield    = fc;

                           double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );


                           //Compute dilation parameter
                           double AGParam = this->ft * yieldHardTwo * 3 / this->fc + myield / 2;
                           double BGParam =
                               yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
                               ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

                           double R = ( sig - ftYield / 3. * yieldHardTwo ) / ( fcYield * BGParam );
                           double mQ = AGParam * exp(R);

                           //Compute the derivative of mQ with respect to kappa

                           //Derivative of AGParam
                           double dAGParamDKappa = dYieldHardTwoDKappa * 3. * this->ft / this->fc;

                           //Derivative of BGParam
                           double BGParamTop = yieldHardTwo / 3. * ( 1. + this->ft / this->fc );
                           double BGParamBottom = ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

                           double dBGParamTopDKappa = dYieldHardTwoDKappa / 3. * ( 1. + this->ft / this->fc );
                           double dBGParamBottomDKappa = 1. / AGParam * dAGParamDKappa - 3. * dYieldHardTwoDKappa / ( 3 * yieldHardTwo + myield / 2. );
                           double dBGParamDKappa = ( dBGParamTopDKappa * BGParamBottom - BGParamTop * dBGParamBottomDKappa ) / pow(BGParamBottom, 2.);

                           //Derivative of R
                           double RTop = ( sig - ftYield / 3. * yieldHardTwo );
                           double RBottom = fcYield * BGParam;
                           double dRTopDKappa = -this->ft / 3. * dYieldHardTwoDKappa;
                           double dRBottomDKappa = this->fc * dBGParamDKappa;
                           double dRDKappa = ( dRTopDKappa * RBottom - RTop * dRBottomDKappa ) / pow(RBottom, 2.);

                           double dMQDKappa = dAGParamDKappa * exp(R) + AGParam * dRDKappa * exp(R);

                           double Bl = sig / fcYield + rho / ( fcYield * sqrt(6.) );

                           double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fcYield;

                           double dAlDYieldHard = -pow(Bl, 2.);

                           const double dDGDSigDKappa =
                               ( -4. * Al * Bl / fcYield + 4. * ( 1 - yieldHardOne ) / fcYield * dAlDYieldHard * Bl ) * dYieldHardOneDKappa +
                               dYieldHardOneDKappa * 2 * yieldHardOne * mQ / fcYield + pow(yieldHardOne, 2.) * dMQDKappa / fcYield;

                           //dDGDSigDKappa = 0.;

                           const double dDGDRhoDKappa =
                               ( dAlDYieldHard / ( sqrt(6.) * fcYield) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) -
                                   4. * Al / ( sqrt(6.) * fcYield ) * Bl + myield / ( sqrt(6.) * fc ) * 2 * yieldHardOne ) * dYieldHardOneDKappa;

                           return {
                               dDGDSigDKappa, dDGDRhoDKappa
                           };
}

FloatArrayF < 6 >
ConcreteDPM2PlasticRate1::computeDDGDStressDKappa(const FloatArrayF < 6 > & stress, double tempKappa, GaussPoint *gp,const double dt) const
{
                           FloatArrayF < 6 > answer;

                           auto tmp = computeDeviatoricVolumetricSplit(stress);
                           auto deviatoricStress = tmp.first;
                           double sig = tmp.second;

                           double rho = computeSecondCoordinate(deviatoricStress);
                           //  double theta = computeThirdCoordinate(deviatoricStress);

                           auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, tempKappa,gp,dt);
                           auto dSigDStress = computeDSigDStress();
                           auto dRhoDStress = computeDRhoDStress(stress);

                           answer = dSigDStress;
                           answer *= dDGDInvDKappa.at(1);

                           FloatArrayF < 6 > temp1;
                           temp1 = dRhoDStress;
                           temp1 *= dDGDInvDKappa.at(2);

                           answer += temp1;

                           return answer;
}


FloatMatrixF < 2, 2 >
ConcreteDPM2PlasticRate1::computeDDGDDInv(double sig,
    double rho,
    double tempKappa,
    GaussPoint *gp,
    const double dt) const
{
                           //compute yieldHardOne and yieldSoft
                           double yieldHardOne = computeHardeningOne(tempKappa);
                           double yieldHardTwo = computeHardeningTwo(tempKappa);

                           double ftYield    = computeFtYield( dt, gp);
                           double fcYield    = computeFcYield( dt, gp);

                           //double ftYield    = ft;
                           //double fcYield    = fc;

                           double myield = 3. * ( pow(fcYield, 2.) - pow(ftYield, 2.) ) / ( fcYield * ftYield ) * this->ecc / ( this->ecc + 1. );

                           //CoQpute dilation parameter
                           double AGParam = this->ft * yieldHardTwo * 3 / this->fc + myield / 2;
                           double BGParam =
                               yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
                               ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

                           double R = ( sig - ftYield / 3. * yieldHardTwo ) / fcYield / BGParam;

                           double dMQDSig = AGParam / ( BGParam * fcYield ) * exp(R);

                           //compute help parameter Al and Bl and the corresponding derivatives
                           double Bl = sig / fcYield + rho / ( fcYield * sqrt(6.) );

                           double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) +
                               sqrt(3. / 2.) * rho / fcYield;

                           double dAlDSig = 2. * ( 1. - yieldHardOne ) * Bl / fcYield;
                           double dBlDSig = 1. / fcYield;

                           double dAlDRho = 2. * ( 1. - yieldHardOne ) * Bl / ( fcYield * sqrt(6.) ) + sqrt(3. / 2.) / fcYield;
                           double dBlDRho = 1. / ( fcYield * sqrt(6.) );

                           //compute second derivatives of g
                           double ddgddSig = 4. * ( 1. - yieldHardOne ) / fcYield * ( dAlDSig * Bl + Al * dBlDSig ) +
                               pow(yieldHardOne, 2.) * dMQDSig / fcYield;
                           // ddgddSig = 0.;

                           double ddgddRho = dAlDRho / ( sqrt(6.) * fcYield ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                               Al * dBlDRho * 4. * ( 1. - yieldHardOne ) / ( sqrt(6.) * fcYield );

                           double ddgdSigdRho = 4. * ( 1. - yieldHardOne ) / fcYield * ( dAlDRho * Bl + Al * dBlDRho );
                           // ddgdSigdRho = 0.;

                           double ddgdRhodSig = dAlDSig / ( sqrt(6.) * fcYield ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. )
                               + Al / ( sqrt(6.) * fcYield ) * ( 4. * ( 1. - yieldHardOne ) * dBlDSig );
                           // ddgdRhodSig = 0.;

                           FloatMatrixF < 2, 2 > answer;
                           answer.at(1, 1) = ddgddSig;
                           answer.at(1, 2) = ddgdSigdRho;
                           answer.at(2, 1) = ddgdRhodSig;
                           answer.at(2, 2) = ddgddRho;
                           return answer;
}






int
ConcreteDPM2PlasticRate1::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    ConcreteDPM2PlasticRate1Status *status = static_cast < ConcreteDPM2PlasticRate1Status * > ( this->giveStatus(gp) );
}

MaterialStatus *
ConcreteDPM2PlasticRate1::CreateStatus(GaussPoint *gp) const
{
    return new ConcreteDPM2PlasticRate1Status(gp);
}
}// end namespace oofem

