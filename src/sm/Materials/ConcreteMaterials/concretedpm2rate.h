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

#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define _IFT_ConcreteDPM2Rate_Name "idm1rate"
#define _IFT_ConcreteDPM2Rate_strengthratetype "sratetype"


namespace oofem {
/**
 * This class implements associated Material Status to ConcreteDPM2Rate.
 * Stores the characteristic length of the element.
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2Rate : public ConcreteDPM2Status
{
protected:

    FloatArray reducedStrain;
    FloatArray tempReducedStrain;

    double kappaOne = 0.;
    double kappaTwo = 0.;
    double tempKappaOne = 0.;
    double tempKappaTwo = 0.;

    double beta = 0.;
    double tempBeta = 0.;

    double strainRate = 0.;
    double tempStrainRate = 0.;
       
    double rateFactor = 0.;
    double tempRateFactor = 0.;

    double DamageTension = 0.0;
    double tempDamageTension = 0.0;

    double KappaDTension = 0.0;
    double tempKappaDTension = 0.0;

    double KappaDTensionOne = 0.0;
    double tempKappaDTensionOne = 0.0;

    double KappaDTensionTwo = 0.0;
    double tempKappaDTensionTwo = 0.0;

public:
    /// Constructor
    ConcreteDPM2RateStatus(GaussPoint *g);

    void updateYourself(TimeStep *tStep) override;
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveKappaDTensionOne() const { return KappaDTensionOne; }
    double giveKappaDTensionTwo() const { return KappaDTensionTwo; }
    double giveKappaOne() const { return kappaOne; }
    double giveKappaTwo() const { return kappaTwo; }
    double giveBeta() const { return beta;}
    double giveStrainRate() const {return strainRate;}    
    double giveRateFactor() const {return rateFactor;}
    
    double giveTempKappaDTensionOne() const { return tempKappaDTensionOne; }   
    double giveTempKappaDTensionTwo() const { return tempKappaDTensionTwo; }
    double giveTempKappaOne() const { return tempKappaOne; }   
    double giveTempKappaTwo() const { return tempKappaTwo; }
    double giveTempBeta() const { return tempBeta;}
    double giveTempStrainRate() const { return tempStrainRate;}
    double giveTempRateFactor() const { return tempRateFactor;}
    
    void setTempKappaDTensionOne(double newKappaDTensionOne) { tempKappaDTensionOne = newKappaDTensionOne; }
    void setTempKappaDTensionTwo(double newKappaDTensionTwo) { tempKappaTwo = newKappaDTensionTwo; }
    void setTempKappaOne(double newKappaOne) { tempKappaOne = newKappaOne; }
    void setTempKappaTwo(double newKappaTwo) { tempKappaTwo = newKappaTwo; }
    void setTempBeta(double newBeta) { tempBeta = newBeta; }
    void setTempStrainRate(double newStrainRate) { tempStrainRate = newStrainRate; }
    void setTempRateFactor(double newRateFactor) { tempRateFactor = newRateFactor; }

    const char *giveClassName() const override { return "ConcreteDPM2RateStatus"; }

    /**
     * Get the reduced strain vector from the material status.
     * @return Strain vector.
     */
    const FloatArray &giveReducedStrain() const { return reducedStrain; }

    void letTempReducedStrainBe(const FloatArray &v)
    { tempReducedStrain = v; }
};



/**
 * This class implements rate dependence for ConcreteDPM2
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2Rate : public ConcreteDPM2
{

public:
    /// Constructor
    ConcreteDPM2RateStatus(GaussPoint *gp);
    

    const char *giveClassName() const override { return "ConcreteDPM2Rate"; }

    void initTempStatus() override;

    void giveRealStressVector(const FloatArrayF< 6 > &strain,
                                       const FloatMatrixF< 6, 6 > &D,
                                       double deltaTime,
                                       GaussPoint *gp,
                                       TimeStep *tStep,
                                       double tempAlpha,
                                       const FloatArrayF< 6 > &effectiveStress) override;

    double computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    
    double computeRateFactor(const double strainRate, GaussPoint *gp, TimeStep *deltaTime) const;


    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    MaterialStatus *giveStatus(GaussPoint *gp) const override;


protected:
};
} // end namespace oofem
#endif // concretedpm2rate_h
