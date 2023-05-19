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

#ifndef concretedpm2rateplastic_h
#define concretedpm2rateplastic_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define _IFT_ConcreteDPM2RatePlastic_Name "con2dpmrateplastic"

#define _IFT_ConcreteDPM2RatePlastic_cCompression "cc"
#define _IFT_ConcreteDPM2RatePlastic_kappaRate0Compression "k0ratec"

#define _IFT_ConcreteDPM2RatePlastic_cTension "ct"
#define _IFT_ConcreteDPM2RatePlastic_Kapparate0tension "k0ratet"



namespace oofem {
/**
 * This class implements associated Material Status to ConcreteDPM2Rate.
 * Stores the characteristic length of the element.
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2RatePlasticStatus : public ConcreteDPM2Status
{
protected:

public:
    /// Constructor
    ConcreteDPM2RatePlasticStatus(GaussPoint *g);

    void updateYourself(TimeStep *tStep) override;
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "ConcreteDPM2RatePlasticStatus"; }
};



/**
 * This class implements rate dependence in the plasticity part for ConcreteDPM2
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2RatePlastic : public ConcreteDPM2
{
protected:
    double cCompression;
    double kappaRate0Compression;

    double cTension;
    double kappaRate0Tension;


public:
    /// Constructor
    ConcreteDPM2RatePlastic(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "ConcreteDPM2RatePlastic"; }

    FloatArrayF< 2 >computeDamage(const FloatArrayF< 6 > &strain, const FloatMatrixF< 6, 6 > &D, double deltaTime, GaussPoint *gp, TimeStep *tStep, double tempAlpha, const FloatArrayF< 6 > &effectiveStress) const override;

   // double computeRateFactor(double alpha) const;

    double computeDeltaPlasticStrainNormTensionP(double tempKappaD, double kappaD, GaussPoint *gp) const;

    double computeDeltaPlasticStrainNormCompressionP(double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp, const double rho) const;

    int checkForUnAndReloadingP(double &tempEquivStrain, double &minEquivStrain, const FloatMatrixF< 6, 6 > &D, GaussPoint *gp) const;

    double computeEquivalentStrainP(double sig, double rho, double theta, GaussPoint *gp) const;

    double computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, GaussPoint *gp) const;
    double computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld, GaussPoint *gp) const;

    void initDamagedP(double kappaD, const FloatArrayF< 6 > &strain, GaussPoint *gp) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatMatrixF< 6, 6 > &D, const FloatArrayF< 6 > &strain, const double deltaTime) const;

    double performRegularReturn(FloatArrayF< 6 > &effectiveStress, ConcreteDPM2_ReturnResult &returnResult, ConcreteDPM2_ReturnType &returnType, double kappaP, GaussPoint *gp, double theta, const double deltaTime) const;

    double performVertexReturn(FloatArrayF< 6 > &effectiveStress, ConcreteDPM2_ReturnResult &returnResult, ConcreteDPM2_ReturnType &returnType, double apexStress, double tempKappaP, GaussPoint *gp, const double deltaTime) const;



    FloatMatrixF< 4, 4 >computeJacobian(double sig, double rho, double theta, double kappa, double deltaLambda, GaussPoint *gp, const double deltaTime) const;

    double computeYieldValue(double sig, double rho, double theta, double tempKappa, const double deltaTime, GaussPoint *gp) const;


    double computeDFDKappa(double sig, double rho, double theta, double tempKappa, double deltaTime, GaussPoint *gp) const;

    FloatArrayF < 2 >computeDFDInv(double sig, double rho, double theta, double tempKappa, GaussPoint * gp) const;

    FloatArrayF < 6 >computeDFDStress(const FloatArrayF < 6 > & stress, double tempKappa, GaussPoint * gp) const;

    FloatMatrixF< 8, 8 >computeFullJacobian(const FloatArrayF< 6 > &stress, const double deltaLambda, GaussPoint *gp, TimeStep *atTime, const double tempKappa, const double deltaTime) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;




protected:
};
} // end namespace oofem
#endif
