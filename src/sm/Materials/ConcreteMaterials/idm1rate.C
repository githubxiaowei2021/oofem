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
REGISTER_Material(IsotropicDamageMaterial1Rate);

IsotropicDamageMaterial1Rate :: IsotropicDamageMaterial1Rate(int n, Domain *d) :
     IsotropicDamageMaterial1(n, d),
    RandomMaterialExtensionInterface()
{
    // deleted by parent, where linearElasticMaterial instance declared
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}




double
IsotropicDamageMaterial1Rate :: computeDamageParam(double kappa, const FloatArray &strain, GaussPoint *gp) const
{
    if ( this->softType == ST_Disable_Damage ) { //dummy material with no damage
        return 0.;
    } else if ( isCrackBandApproachUsed() ) { // adjustment of softening law according to the element size, given crack opening or fracture energy
        return computeDamageParamForCohesiveCrack(kappa, gp);
    } else { // no adjustment according to element size, given fracturing strain
        return damageFunction(kappa, gp);
    }
}

double
IsotropicDamageMaterial1Rate :: computeDamageParamForCohesiveCrack(double kappa, GaussPoint *gp, TimeStep *tStep, double tempAlpha) const
{
    const double e0 = this->give(e0_ID, gp);  // e0 is the strain at the peak stress
    const double E = this->linearElasticMaterial->give('E', gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);     // wf is the crack opening
    double omega = 0.0;
    double rateFactor = computeRateFactor(tempAlpha, gp, tStep);
    

    if ( kappa > e0 ) {
        if ( this->gf != 0. ) { //cohesive crack model
            if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
                wf = this->gf / E / e0; // wf is the crack opening
            } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi)linear softening law
                wf = 2. * gf / E / e0; // wf is the crack opening
            } else {
                OOFEM_ERROR("Gf unsupported for softening type softType = %d", softType);
            }
        } else if ( softType == ST_BiLinear_Cohesive_Crack ) {
            wf = this->wk / ( e0 * E - this->sk ) * ( e0 * E );
        }


        auto status = static_cast< IsotropicDamageMaterial1RateStatus * >( this->giveStatus(gp) );
        double Le = status->giveLe();
        double ef = wf / Le;    //ef is the fracturing strain /// FIXME CHANGES BEHAVIOR!
        if ( ef < e0 ) { //check that no snapback occurs
            double minGf = 0.;
            OOFEM_WARNING("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f", ef, e0, gp->giveElement()->giveNumber(), Le);
            if ( gf != 0. ) { //cohesive crack
                if ( softType == ST_Exponential_Cohesive_Crack ) { //exponential softening
                    minGf = E * e0 * e0 * Le;
                } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { //(bi)linear softening law
                    minGf = E * e0 * e0 * Le / 2.;
                } else {
                    OOFEM_WARNING("Gf unsupported for softening type softType = %d", softType);
                }

                if ( checkSnapBack ) {
                    OOFEM_ERROR("Material number %d, decrease e0, or increase Gf from %f to Gf=%f", this->giveNumber(), gf, minGf);
                }
            }

            if ( checkSnapBack ) { //given fracturing strain
                OOFEM_ERROR("Material number %d, increase ef %f to minimum e0 %f", this->giveNumber(), ef, e0);
            }
        }

        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            if ( kappa < ef ) {
	      omega = ( ef / ( kappa * rateFactor) ) * ( kappa * rateFactor  - e0 ) / ( ef - e0 );
            } else {
                omega = 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
            }
        } else if (  this->softType == ST_BiLinear_Cohesive_Crack ) {
            double gft = this->give(gft_ID, gp);
            double ef, sigmak, epsf, ek;
            if ( gft > 0.0 ) {
                ek = this->give(ek_ID, gp);
                ef = 2 * gf / E / e0 / Le; //the first part corresponds to linear softening
                sigmak = E * e0 * ( ef - ek ) / ( ef - e0 );
                epsf = 2 * ( gft - gf ) / sigmak / Le + ef;

                if ( gft < gf ) {
                    OOFEM_ERROR("The total fracture energy gft %f must be greater than the initial fracture energy gf %f", gft, gf);
                }
            } else {
                ek     = this->wk / Le + ( this->sk ) / E;
                ef     = ( this->wk / ( e0 * E - this->sk ) * ( e0 * E ) ) / Le;
                sigmak = this->sk;
                epsf   = this->wf / Le;
            }
            if ( ( ek > ef ) || ( ek < e0 ) ) {
                OOFEM_WARNING("ek %f is not between e0 %f and ef %f", ek, e0, ef);
            }

            if ( kappa <= ek ) {
	      omega = 1.0 - ( ( e0 / ( kappa * rateFactor ) ) * ( ek - kappa * rateFactor ) / ( ek - e0 ) + ( ( sigmak / ( E * kappa * rateFactor ) ) * ( kappa * rateFactor - e0 ) / ( ek - e0 ) ) );
            } else if ( kappa > ek && kappa <= epsf ) {
                omega = 1.0 - ( ( sigmak / ( E * kappa * rateFactor ) ) * ( epsf - kappa * rateFactor ) / ( epsf - ek ) );
            } else if ( kappa <= e0 ) {
                omega = 0.0;
            } else {
                omega = maxOmega;
            }
        } else if (  this->softType == ST_Exponential_Cohesive_Crack ) {
            // exponential cohesive crack - iteration needed
            double R, Lhs, help;
            int nite = 0;
            // iteration to achieve objectivity
            // we are looking for a state in which the elastic stress is equal to
            // the stress from crack-opening relation
            // ef has now the meaning of strain
            do {
                nite++;
                help = omega * kappa / ef;
                R = ( 1. - omega ) * kappa - e0 *exp(-help); //residuum
                Lhs = kappa - e0 *exp(-help) * kappa / ef; //- dR / (d omega)
                omega += R / Lhs;
                if ( nite > 40 ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(R) >= e0 * IDM1_ITERATION_LIMIT );
        } else {
            OOFEM_ERROR("Unknown softening type for cohesive crack model.");
        }

        if ( omega > 1.0 ) {
            OOFEM_WARNING("damage parameter is %f, which is greater than 1, snap-back problems", omega);
            omega = maxOmega;
            if ( checkSnapBack ) {
                OOFEM_ERROR("");
            }
        }

        if ( omega < 0.0 ) {
            OOFEM_WARNING("damage parameter is %f, which is smaller than 0, snap-back problems", omega);
            omega = 0.0;
            if ( checkSnapBack ) {
                OOFEM_ERROR("");
            }
        }
    }
    return omega;
}

double
IsotropicDamageMaterial1Rate ::computeRateFactor(double alpha,
                                GaussPoint *gp,
                                TimeStep *tStep) const
{
    this->strengthRateType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->strengthRateType, _IFT_IsotropicDamageMaterial1Rate_strengthratetype);
    if ( this->strengthRateType < 0 || this->strengthRateType > 2 ) {
        OOFEM_ERROR("strengthRateType not implemented. Must be 0, 1 or 2\n");
    }

    if ( this->strengthRateType == 0 ) {
        return 1;
    }

    auto status = static_cast< IsotropicDamageMaterial1RateStatus * >( this->giveStatus(gp) );

    //Determine the principal values of the strain
    auto principalStrain = StructuralMaterial::computePrincipalValues(from_voigt_strain);   ///@todo CHECK

    //Determine max and min value;
    double maxStrain = -1.e20, minStrain = 1.e20;
    for ( int k = 1; k <= principalStrain.giveSize(); k++ ) {
        //maximum
        if ( principalStrain.at(k) > maxStrain ) {
            maxStrain = principalStrain.at(k);
        }

        //minimum
        if ( principalStrain.at(k) < minStrain ) {
            minStrain = principalStrain.at(k);
        }
    }


    //Tension
    //For tension according to Model Code 2010
    double rateFactorTension = 1.;
    double strainRateRatioTension = strainRate / 1.e-6;

    if ( this->strengthRateType == 1 ) {
        if ( strainRate < 1.e-6 ) {
            rateFactorTension = 1.;
        } else if ( 1.e-6 < strainRate ) {
            rateFactorTension = pow(strainRateRatioTension, 0.018);
        }
    } else if ( this->strengthRateType == 2 ) {
        if ( strainRate < 1.e-6 ) {
            rateFactorTension = 1.;
        } else if ( 1.e-6 < strainRate && strainRate < 10 ) {
            rateFactorTension = pow(strainRateRatioTension, 0.018);
        } else {
            rateFactorTension =  0.0062 * pow(strainRateRatioTension, 1. / 3.);
        }
    }

    //For compression according to Model Code 2010
    double rateFactorCompression = 1.;
    double strainRateRatioCompression = strainRate / ( -30.e-6 );
    if ( this->strengthRateType == 1 ) {
        if ( strainRate > -30.e-6 ) {
            rateFactorCompression = 1.;
        } else if ( -30.e-6 > strainRate ) {
            rateFactorCompression = pow(strainRateRatioCompression, 0.014);
        }
    } else if ( this->strengthRateType == 2 ) {
        if ( strainRate > -30.e-6 ) {
            rateFactorCompression = 1.;
        } else if ( -30.e-6 > strainRate && strainRate > -30 ) {
            rateFactorCompression = pow(strainRateRatioCompression, 0.014);
        } else if ( -30 > strainRate && strengthRateType == 2 ) {
            rateFactorCompression =  0.012 * pow(strainRateRatioCompression, 0.333);
        }
    }

    double rateFactor = ( 1. - alpha ) * rateFactorTension + alpha * rateFactorCompression;

    return rateFactor;
}

