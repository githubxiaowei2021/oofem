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
    IDM1RateStatus *status = static_cast< IDM1RateStatus * >( this->giveStatus(gp) );
    //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, oldEquivStrain, tempKappa = 0.0, omega = 0.0, deltaTime = 0.,strainRate = 0.;

    this->initTempStatus(gp);

    // subtract stress-independent part
    // note: eigenStrains (temperature) are present in strains stored in gp
    // therefore it is necessary to subtract always the total eigen strain value

    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);
   
    equivStrain = this->computeEquivalentStrain(reducedTotalStrainVector, gp, tStep);

    FloatArray oldReducedStrain = status->giveReducedStrain();
    
    //Peter: computeEquivalentStrain again. However, this time from the old strain.
    oldEquivStrain = this->computeEquivalentStrain(oldReducedStrain, gp, tStep);    
    
    f = equivStrain - status->giveKappa();
    double rateFactor = 1;
    if ( f <= 0.0 ) {
      // damage does not grow
      tempKappa = status->giveKappa();
      omega     = status->giveDamage();
    } else {
      // damage grows
      //Calculate strain rate.
      if ( tStep->giveTimeIncrement() == 0 ) { //Problem with the first step. For some reason the time increment is zero
	deltaTime = 1.;
        } else {
            deltaTime = tStep->giveTimeIncrement();
        }
      //This needs to be extended to cracking!
      strainRate = (equivStrain-oldEquivStrain)/deltaTime;
      rateFactor = computeRateFactor(strainRate, gp, tStep);     

      printf("strainRate = %e, rateFactor = %e\n", strainRate, rateFactor);
      
      tempKappa = f/rateFactor + status->giveKappa();
      this->initDamaged(tempKappa, reducedTotalStrainVector, gp);

      // evaluate damage parameter
      omega = this->computeDamageParam(tempKappa, reducedTotalStrainVector, gp);
    }
    
    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);


    if ( ( reducedTotalStrainVector.giveSize() > 1 ) || ( reducedTotalStrainVector.at(1) > 0. ) ) {
        //emj
        de.times(1.0 - omega);
    }

    answer.beProductOf(de, reducedTotalStrainVector);

    // update gp
    status->letTempReducedStrainBe(reducedTotalStrainVector);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);

#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
}


double
IDM1Rate ::computeRateFactor(double strainRate,
			     GaussPoint *gp,
			     TimeStep *tStep) const
{
  double rateFactor = 1.;
  
    if ( this->strengthRateType == 0 ) {
        return 1;
    }

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



MaterialStatus *
IDM1Rate :: CreateStatus(GaussPoint *gp) const
{
    return new IDM1RateStatus(gp);
}

MaterialStatus *
IDM1Rate :: giveStatus(GaussPoint *gp) const
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

  
  IDM1RateStatus :: IDM1RateStatus(GaussPoint *g) : IsotropicDamageMaterial1Status(g), reducedStrain(), tempReducedStrain()
{
  
    int rsize = IDM1Rate :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    reducedStrain.resize(rsize);

    // reset temp vars.
    tempReducedStrain = reducedStrain;
}

void
IDM1RateStatus :: initTempStatus()
{
    IsotropicDamageMaterial1Status :: initTempStatus();
    this->tempReducedStrain = this->reducedStrain;
}


void
IDM1RateStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    // Call corresponding function of the parent class to print
    IsotropicDamageMaterial1Status::printOutputAt(file, tStep);
}
  
void
IDM1RateStatus :: updateYourself(TimeStep *tStep)
{
    IsotropicDamageMaterial1Status :: updateYourself(tStep);
    this->reducedStrain = this->tempReducedStrain;
}


void
IDM1RateStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    IsotropicDamageMaterial1Status :: saveContext(stream, mode);

    contextIOResultType iores;
    
    if ( ( iores = reducedStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

}

void
IDM1RateStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    IsotropicDamageMaterial1Status :: restoreContext(stream, mode);

    contextIOResultType iores;
    
    if ( ( iores = reducedStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }    

}
  
}// end namespace oofem
