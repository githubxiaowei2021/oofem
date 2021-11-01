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

#ifndef idm1rate_h
#define idm1rate_h

#include "material.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "sm/Materials/ConcreteMaterials/idm1.h"
#include "sm/Materials/structuralms.h"


#define _IFT_IDM1Rate_Name "idm1rate"
#define _IFT_IDM1Rate_strengthratetype "sratetype"


namespace oofem {
/**
 * This class implements associated Material Status to IDM1Rate.
 * Stores the characteristic length of the element.
 * @author: Xiaowei Liu, Peter Grassl
 */
class IDM1RateStatus : public IsotropicDamageMaterial1Status
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


public:
    /// Constructor
    IDM1RateStatus(GaussPoint *g);

    void updateYourself(TimeStep *tStep) override;
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveKappaOne() const { return kappaOne; }
    double giveKappaTwo() const { return kappaTwo; }
    double giveBeta() const { return beta;}
    double giveStrainRate() const {return strainRate;}    
    double giveRateFactor() const {return rateFactor;}
    
    double giveTempKappaOne() const { return tempKappaOne; }   
    double giveTempKappaTwo() const { return tempKappaTwo; }
    double giveTempBeta() const { return tempBeta;}
    double giveTempStrainRate() const { return tempStrainRate;}
    double giveTempRateFactor() const { return tempRateFactor;}
    
    void setTempKappaOne(double newKappaOne) { tempKappaOne = newKappaOne; }
    void setTempKappaTwo(double newKappaTwo) { tempKappaTwo = newKappaTwo; }
    void setTempBeta(double newBeta) { tempBeta = newBeta; }
    void setTempStrainRate(double newStrainRate) { tempStrainRate = newStrainRate; }
    void setTempRateFactor(double newRateFactor) { tempRateFactor = newRateFactor; }

    const char *giveClassName() const override { return "IDM1RateStatus"; }

    /**
     * Get the reduced strain vector from the material status.
     * @return Strain vector.
     */
    const FloatArray &giveReducedStrain() const { return reducedStrain; }

    void letTempReducedStrainBe(const FloatArray &v)
    { tempReducedStrain = v; }
};



/**
 * This class implements rate dependence for idm1
 * @author: Xiaowei Liu, Peter Grassl
 */
class IDM1Rate : public IsotropicDamageMaterial1
{
protected:

    /** Type of strength strain rate dependence used.
     * 0 = no strain rate (default)
     * 1 = Model Code 2010 initial and second branch of strain rate effect for strength
     * 2 = Substitute crack opening for strain rate
     */
    int strengthRateType = 0;

public:
    /// Constructor
    IDM1Rate(int n, Domain *d);
    /// Destructor
    virtual ~IDM1Rate();

    const char *giveClassName() const override { return "IDM1Rate"; }

    void initializeFrom(InputRecord &ir) override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) override;

    double computeDamageParameter(double tempKappaOne, double tempKappaTwo, GaussPoint *gp) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    
    double computeRateFactor(const double strainRate, GaussPoint *gp, TimeStep *deltaTime) const;


    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    MaterialStatus *giveStatus(GaussPoint *gp) const override;


protected:
};
} // end namespace oofem
#endif // idm1rate_h
