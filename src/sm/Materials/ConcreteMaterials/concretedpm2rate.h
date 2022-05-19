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

#ifndef concretedpm2rate_h
#define concretedpm2rate_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define _IFT_ConcreteDPM2Rate_Name "con2dpmrate"

#define _IFT_ConcreteDPM2Rate_atOne "atone"
#define _IFT_ConcreteDPM2Rate_atTwo "attwo"
#define _IFT_ConcreteDPM2Rate_atThree "atthree"
#define _IFT_ConcreteDPM2Rate_atFour "atfour"
#define _IFT_ConcreteDPM2Rate_atFive "atfive"

#define _IFT_ConcreteDPM2Rate_acOne "acone"
#define _IFT_ConcreteDPM2Rate_acTwo "actwo"
#define _IFT_ConcreteDPM2Rate_acThree "acthree"
#define _IFT_ConcreteDPM2Rate_acFour "acfour"
#define _IFT_ConcreteDPM2Rate_acFive "acfive"


namespace oofem {
/**
 * This class implements associated Material Status to ConcreteDPM2Rate.
 * Stores the characteristic length of the element.
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2RateStatus : public ConcreteDPM2Status
{
protected:

    double beta = 0.;
    double tempBeta = 0.;

    double strainRateTension = 0.;
    double tempStrainRateTension = 0.;

    double strainRateCompression = 0.;
    double tempStrainRateCompression = 0.;    
    
    double rateFactorTension = 0.;
    double tempRateFactorTension = 0.;

    double rateFactorCompression = 0.;
    double tempRateFactorCompression = 0.;
    
public:
    /// Constructor
    ConcreteDPM2RateStatus(GaussPoint *g);

    void updateYourself(TimeStep *tStep) override;
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveBeta() const { return beta;}
    double giveStrainRateTension() const {return strainRateTension;}    
    double giveRateFactorTension() const {return rateFactorTension;}
    double giveStrainRateCompression() const {return strainRateCompression;}    
    double giveRateFactorCompression() const {return rateFactorCompression;}
    
    double giveTempBeta() const { return tempBeta;}
    double giveTempStrainRateTension() const { return tempStrainRateTension;}
    double giveTempRateFactorTension() const { return tempRateFactorTension;}
    double giveTempStrainRateCompression() const { return tempStrainRateCompression;}
    double giveTempRateFactorCompression() const { return tempRateFactorCompression;}

    void setTempBeta(double newBeta) { tempBeta = newBeta; }
    void setTempStrainRateTension(double newStrainRate) { tempStrainRateTension = newStrainRate; }
    void setTempRateFactorTension(double newRateFactor) { tempRateFactorTension = newRateFactor; }
    void setTempStrainRateCompression(double newStrainRate) { tempStrainRateCompression = newStrainRate; }
    void setTempRateFactorCompression(double newRateFactor) { tempRateFactorCompression = newRateFactor; }

    const char *giveClassName() const override { return "ConcreteDPM2RateStatus"; }

};



/**
 * This class implements rate dependence for ConcreteDPM2
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2Rate : public ConcreteDPM2
{
 protected:
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
    ConcreteDPM2Rate(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "ConcreteDPM2Rate"; }

    
    FloatArrayF< 2 >computeDamage(const FloatArrayF< 6 > &strain, const FloatMatrixF< 6, 6 > &D, double timeFactor, GaussPoint *gp, TimeStep *tStep, double alpha, const FloatArrayF< 6 > &effectiveStress) const override;

    double computeDamageParamTension(double equivStrain, double KappaDTensionOne, double KappaDTensionTwo, double le, double omegaOld) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    
    double computeRateFactorTension(const double strainRate, GaussPoint *gp, TimeStep *deltaTime) const;

    double computeRateFactorCompression(const double strainRate, GaussPoint *gp, TimeStep *deltaTime) const;

    /// Compute damage parameter in compression.
    double computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld) const;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem
#endif // concretedpm2rate_h
