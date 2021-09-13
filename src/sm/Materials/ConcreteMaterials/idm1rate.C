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

#include "idm1rate.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "strainvector.h"
#include "stressvector.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"
#include "crosssection.h"


namespace oofem {
REGISTER_Material(IDM1Rate);


IDM1Rate :: IDM1Rate(int n, Domain *d) : IsotropicDamageMaterial1(n, d)
{

}


IDM1Rate :: ~IDM1Rate()
{

}

void
IDM1Rate :: initializeFrom(InputRecord &ir)
{
  IsotropicDamageMaterial1 :: initializeFrom(ir);

  this->strengthRateType = 0;
  IR_GIVE_OPTIONAL_FIELD(ir, this->strengthRateType, _IFT_IDM1Rate_strengthratetype);
  if ( this->strengthRateType < 0 || this->strengthRateType > 2 ) {
    OOFEM_ERROR("strengthRateType not implemented. Must be 0, 1 or 2\n");
  }

}



void
IDM1Rate :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    this->initTempStatus(gp);

    // subtract stress-independent part
    // note: eigenStrains (temperature) are present in strains stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);

    //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);

    // compute equivalent strain
    equivStrain = this->computeEquivalentStrain(reducedTotalStrainVector, gp, tStep);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = equivStrain - status->giveKappa();
double rateFactor = 1;
        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            omega     = status->giveDamage();
        } else {
            // damage grows
            //tempKappa = equivStrain;
            rateFactor = computeRateFactor(reducedTotalStrainVector, 0, gp, tStep);
            tempKappa = f/rateFactor + status->giveKappa();
            this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            omega = this->computeDamageParam(tempKappa, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = equivStrain;
        this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
        omega = this->computeDamageParam(tempKappa, reducedTotalStrainVector, gp);
        if ( omega < status->giveDamage() ) {
            // unloading takes place
            omega = status->giveDamage();
        }
    } else {
        OOFEM_ERROR("unsupported loading/unloading criterion");
    }


    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    //mj
    // permanent strain - so far implemented only in 1D
    if ( permStrain && reducedTotalStrainVector.giveSize() == 1 ) {
        double epsp = evaluatePermanentStrain(tempKappa, omega);
        reducedTotalStrainVector.at(1) -= epsp;
    }
    // damage deactivation in compression for 1D model
    if ( ( reducedTotalStrainVector.giveSize() > 1 ) || ( reducedTotalStrainVector.at(1) > 0. ) ) {
        //emj
        de.times(1.0 - omega);
    }

    answer.beProductOf(de, reducedTotalStrainVector);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
}


double
IDM1Rate ::computeRateFactor(FloatArray &strain,
			     double alpha,
			     GaussPoint *gp,
			     TimeStep *tStep) const
{

    // Peter: I don't think that you need to have multiple strain rate factors. Tension is sufficient, since the model uses only one damage variable. Also, I suggest to have the strain as an input into the function (as it is done above), because reduced strain is not stored in the status of imd1rate.
    if ( this->strengthRateType == 0 ) {
        return 1;
    }

    auto status = static_cast<IsotropicDamageMaterialStatus *>( this->giveStatus( gp ) );

    // Determine the principal values of the strain
    auto principalStrain = StructuralMaterial::computePrincipalValues( from_voigt_strain( strain ) ); ///@todo CHECK

    double strainRate;

    strainRate = ( principalStrain - status->giveprincipalStrain() ) / tStep;

    double rateFactor = 1.;

    // Peter: Add here the calculation of the rate factor
    double strainRateRatio = strainRate / 1.e-6;

    if ( this->strengthRateType == 1 ) {
        if ( strainRate < 1.e-6 ) {
            rateFactor = 1.;
        } else if ( 1.e-6 < strainRate ) {
            rateFactor = pow( strainRateRatio, 0.018 );
        }
    } else if ( this->strengthRateType == 2 ) {
        if ( strainRate < 1.e-6 ) {
            rateFactor = 1.;
        } else if ( 1.e-6 < strainRate && strainRate < 10 ) {
            rateFactor = pow( strainRateRatio, 0.018 );
        } else {
            rateFactor = 0.0062 * pow( strainRateRatio, 1. / 3. );
        }
    }

    return rateFactor;
}


}
// end namespace oofem
