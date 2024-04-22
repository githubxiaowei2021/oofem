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

#ifndef concretedpm2Plasticrate1_h
#define concretedpm2Plasticrate1_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define _IFT_ConcreteDPM2PlasticRate1_Name "con2dpmPlasticrate1"

#define _IFT_ConcreteDPM2PlasticRate1_atOne "atone"
#define _IFT_ConcreteDPM2PlasticRate1_atTwo "attwo"
#define _IFT_ConcreteDPM2PlasticRate1_atThree "atthree"
#define _IFT_ConcreteDPM2PlasticRate1_atFour "atfour"
#define _IFT_ConcreteDPM2PlasticRate1_atFive "atfive"

#define _IFT_ConcreteDPM2PlasticRate1_acOne "acone"
#define _IFT_ConcreteDPM2PlasticRate1_acTwo "actwo"
#define _IFT_ConcreteDPM2PlasticRate1_acThree "acthree"
#define _IFT_ConcreteDPM2PlasticRate1_acFour "acfour"
#define _IFT_ConcreteDPM2PlasticRate1_acFive "acfive"


/*
#define _IFT_ConcreteDPM2PlasticRate1_cCompression "cc"
#define _IFT_ConcreteDPM2PlasticRate1_kappaRate0Compression "k0ratec"

#define _IFT_ConcreteDPM2PlasticRate1_cTension "ct"
#define _IFT_ConcreteDPM2PlasticRate1_kappaRate0Tension "k0ratet"
*/
namespace oofem {
/**
 * This class implements associated Material Status to ConcreteDPM2Rate.
 * Stores the characteristic length of the element.
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2PlasticRate1Status : public ConcreteDPM2Status
{
protected:

    double RateTension = 0.;
    double tempRateTension = 0.;
    double RateCompression = 0.;
    double tempRateCompression = 0.;

    double beta = 0.;
    double tempBeta = 0.;

    double rateFactorTension = 0.;
    double tempRateFactorTension = 0.;

    double rateFactorCompression = 0.;
    double tempRateFactorCompression = 0.;


public:
    /// Constructor
    ConcreteDPM2PlasticRate1Status(GaussPoint *g);

    void updateYourself(TimeStep *tStep) override;
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveRateTension() const { return RateTension; }
    double giveRateCompression() const { return RateCompression; }
    double giveBeta() const { return beta; }
    double giveRateFactorTension() const {return rateFactorTension;}
    double giveRateFactorCompression() const {return rateFactorCompression;}


    double giveTempRateTension() const { return tempRateTension; }
    double giveTempRateCompression() const { return tempRateCompression; }
    double giveTempBeta() const { return tempBeta; }
    double giveTempRateFactorTension() const { return tempRateFactorTension;}
    double giveTempRateFactorCompression() const { return tempRateFactorCompression;}


    void setTempRateTension(double newRate) { tempRateTension = newRate; }
    void setTempRateCompression(double newRate) { tempRateCompression = newRate; }
    void setTempBeta(double newBeta) { tempBeta = newBeta; }
    void setTempRateFactorTension(double newRateFactor) { tempRateFactorTension = newRateFactor; }
    void setTempRateFactorCompression(double newRateFactor) { tempRateFactorCompression = newRateFactor; }


    const char *giveClassName() const override { return "ConcreteDPM2PlasticRate1Status"; }

};



/**
 * This class implements rate dependence for ConcreteDPM2
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2PlasticRate1 : public ConcreteDPM2
{
 protected:
     //double cCompression;
     //double kappaRate0Compression;

     //double cTension;
     //double kappaRate0Tension;
     double atOne;
     double atTwo;
     double atThree;
     double atFour;
     double atFive;

     double acOne;
     double acTwo;
     double acThree;
     double acFour;
     double acFive;

 public:
    /// Constructor
    ConcreteDPM2PlasticRate1(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "ConcreteDPM2PlasticRate1"; }

    double computeYieldValue(double sig, double rho, double theta, double tempKappa, const double dt, GaussPoint *gp) const;

    double computeFcYield( const double dt, GaussPoint *gp) const;
    double computeFtYield( const double dt, GaussPoint *gp) const;
    double computeRateFactorCompression(GaussPoint *gp,const double dt) const;
    double computeRateFactorTension(GaussPoint *gp, const double dt) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatMatrixF< 6, 6 > &D, const FloatArrayF< 6 > &strain, const double dt) const;

    double performVertexReturn(FloatArrayF< 6 > &effectiveStress, ConcreteDPM2_ReturnResult &returnResult, ConcreteDPM2_ReturnType &returnType, double apexStress, double tempKappaP, GaussPoint *gp, const double dt) const;

    double performRegularReturn(FloatArrayF< 6 > &effectiveStress, ConcreteDPM2_ReturnResult &returnResult, ConcreteDPM2_ReturnType &returnType, double kappaP, GaussPoint *gp, double theta, const double dt) const;

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF< 2 >computeDFDInv(double sig, double rho, double theta, double tempKappa, const double dt, GaussPoint *gp) const;
    double computeEquivalentStrainP(double sig, double rho, double theta, const double dt, GaussPoint *gp) const;
    double computeDFDKappa(double sig, double rho, double theta, double tempKappa, const double dt, GaussPoint *gp) const;

    int checkForUnAndReloadingP(double &tempEquivStrain, double &minEquivStrain, const FloatMatrixF< 6, 6 > &D,const double dt, GaussPoint *gp) const;
    FloatArrayF< 2 >computeDamage(const FloatArrayF< 6 > &strain, const FloatMatrixF< 6, 6 > &D, const double dt, GaussPoint *gp, TimeStep *tStep, double tempAlpha, const FloatArrayF< 6 > &effectiveStress) const override;

    FloatMatrixF< 4, 4 >computeJacobian(double sig, double rho, double theta, double kappa, double deltaLambda, const double dt, GaussPoint *gp) const;
    FloatArrayF< 6 >computeDFDStress(const FloatArrayF< 6 > &stress, double tempKappa, const double dt, GaussPoint *gp) const;
    FloatMatrixF< 8, 8 >computeFullJacobian(const FloatArrayF< 6 > &stress, const double deltaLambda, GaussPoint *gp, TimeStep *atTime, const double tempKappa, const double dt) const;
    FloatMatrixF < 6, 6 >compute3dTangentStiffness(GaussPoint * gp, TimeStep * tStep, const double dt) const;

    double computeDuctilityMeasure(double sig,double rho,double theta,GaussPoint *gp,const double dt) const;
    double computeTempKappa(double kappaInitial,double sigTrial,double rhoTrial,double sig,GaussPoint *gp,const double dt) const;
    double computeDKappaDDeltaLambda(double sig,double rho,double theta,double tempKappa,GaussPoint *gp,const double dt) const;
    FloatArrayF < 2 >computeDDKappaDDeltaLambdaDInv(double sig,double rho,double theta,double tempKappa,GaussPoint *gp,const double dt) const;
    FloatArrayF < 6 >computeDDKappaDDeltaLambdaDStress(const FloatArrayF < 6 > & stress, double tempKappa,GaussPoint *gp,const double dt) const;
    double computeDDKappaDDeltaLambdaDKappa(double sig, double rho, double theta, double tempKappa,GaussPoint *gp, const double dt) const;

    double computeRatioPotential(double sig,double rho,double tempKappa,GaussPoint *gp,const double dt) const;

    FloatArrayF < 2 >computeDDuctilityMeasureDInv(double sig,double rho,double theta,double tempKappa,GaussPoint *gp,const double dt) const;

    FloatArrayF < 2 >computeDGDInv(double sig,double rho,double tempKappa,GaussPoint *gp,const double dt) const;
    FloatArrayF < 6 >computeDGDStress(const FloatArrayF < 6 > & stress, const double tempKappa,GaussPoint *gp,const double dt) const;
    FloatMatrixF < 6, 6 >computeDDGDDStress(const FloatArrayF < 6 > & stress, const double tempKappa,GaussPoint *gp,const double dt) const;

    FloatArrayF < 2 >computeDDGDInvDKappa(double sig,double rho,double tempKappa,GaussPoint *gp,const double dt) const;
    FloatArrayF < 6 >computeDDGDStressDKappa(const FloatArrayF < 6 > & stress, double tempKappa,  GaussPoint *gp,const double dt) const;

    FloatMatrixF < 2, 2 >computeDDGDDInv(double sig,double rho,double tempKappa,GaussPoint *gp,const double dt) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem
#endif // concretedpm2rate_h
