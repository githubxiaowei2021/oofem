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
        //this->KappaRate = this->tempKappaRate;
    }


    void
    ConcreteDPM2RatePlasticStatus::saveContext(DataStream &stream, ContextMode mode)
    {
        ConcreteDPM2Status::saveContext(stream, mode);

        contextIOResultType iores;
    }

    void
    ConcreteDPM2RatePlasticStatus::restoreContext(DataStream &stream, ContextMode mode)
    {
        ConcreteDPM2Status::restoreContext(stream, mode);
        contextIOResultType iores;

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

        this->Ccompression = 1.e-6;
        IR_GIVE_OPTIONAL_FIELD(ir, this->cCompression, _IFT_ConcreteDPM2RatePlastic_Ccompression);
        this->Kapparate0compression = 10;
        IR_GIVE_OPTIONAL_FIELD(ir, this->kappaRate0Compression, _IFT_ConcreteDPM2RatePlastic_Kapparate0compression);

        //this->Ctension = 0.01935;
        this->Ctension = 0.01935;
        IR_GIVE_OPTIONAL_FIELD(ir, this->cTension, _IFT_ConcreteDPM2RatePlastic_Ctension);
        //this->Kapparate0tension = 0.00000434165;
        this->Kapparate0tension = 0.00000434165;
        IR_GIVE_OPTIONAL_FIELD(ir, this->kappaRate0Tension, _IFT_ConcreteDPM2RatePlastic_Kapparate0tension);
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

        double tempKappaP;
        double tempKappaRate = ( tempKappaP - status->giveTempKappaP() ) / deltaTime;


        //status->setTempKappaRate(tempKappaRate);

        //  compute elliptic function r
        double rFunction = ( 4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) +
                             pow( ( 2. * ecc - 1. ), 2.) ) /
                           ( 2. * ( 1. - pow(ecc, 2.) ) * cos(theta) +
                             ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.)
                                                      + 5. * pow(ecc, 2.) - 4. * ecc) );

        double fcyield = this->fc * ( 1 + Ccompression * log(1 + tempKappaRate / Kapparate0compression) );
        double ftyield = this->ft * ( 1 + Ctension * log(1 + tempKappaRate / Kapparate0tension) );
        double myield = 3. * ( pow(fcyield, 2.) - pow(ftyield, 2.) ) / ( fcyield * ftyield ) * this->ecc / ( this->ecc + 1. );

        //compute help function Al
        double Al = ( 1. - yieldHardOne ) * pow( ( sig / fcyield + rho / ( sqrt(6.) * fcyield ) ), 2.) +
                    sqrt(3. / 2.) * rho / fcyield;




        //Compute yield equation
        return pow(Al, 2.) +
               pow(yieldHardOne, 2.) * yieldHardTwo * myield * ( sig / fcyield + rho * rFunction / ( sqrt(6.) * fcyield ) ) -
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


        double tempKappaP;
        double tempKappaRate = ( tempKappaP - status->giveTempKappaP() ) / deltaTime;

        //status->setTempKappaRate(tempKappaRate);

        double fcyield = fc * ( 1 + Ccompression * log(1 + tempKappaRate / Kapparate0compression) );
        double ftyield = ft * ( 1 + Ctension * log(1 + tempKappaRate / Kapparate0tension) );
        double myield = 3. * ( pow(fcyield, 2.) - pow(ftyield, 2.) ) / ( fcyield * ftyield ) * this->ecc / ( this->ecc + 1. );

        //compute elliptic function r
        double rFunction =
            ( 4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + pow( ( 2. * ecc - 1. ), 2) ) /
            ( 2 * ( 1. - pow(ecc, 2) ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + 5. * pow(ecc, 2) - 4. * ecc) );

        //compute help functions Al, Bl
        double Al = ( 1. - yieldHardOne ) * pow( ( sig / fcyield + rho / ( sqrt(6.) * fcyield ) ), 2. ) + sqrt(3. / 2.) * rho / fcyield;


        double Bl = sig / fcyield + rho / ( fcyield * sqrt(6.) );
        double dFDYieldHardOne = -2. * Al * pow(Bl, 2.)
                                 + 2. * yieldHardOne * yieldHardTwo * myield * ( sig / fcyield + rho * rFunction / ( sqrt(6.) * fcyield ) ) - 2. * yieldHardOne * pow(yieldHardTwo, 2.);

        double dFDYieldHardTwo = pow(yieldHardOne, 2.) * myield * ( sig / fcyield + rho * rFunction / ( sqrt(6.) * fcyield ) ) - 2. * yieldHardTwo * pow(yieldHardOne, 2.);



	//Xiaowei, if possible can you give them names so that we know which these terms are? Something like "double dFcDkappaRate" or equivalent.
	//Also, follow the OOFEM coding rules. New word in names starts with captial letter. So, "fcYield, ..." For variables, small letters at the start. There is a document on the OOFEM webpage which explains everything.
	
        double T1 = ft * Ctension * 1 / ( tempKappaRate + Kapparate0tension );
        double C1 = fc * Ccompression * 1 / ( tempKappaRate + Kapparate0compression );

        double mT = ( ( -3 * ecc ) / ( ecc + 1 ) ) * ( ( ftyield * ftyield + fcyield * fcyield ) / ( fcyield * ftyield * ftyield ) );
        double mC = ( ( 3 * ecc ) / ( ecc + 1 ) ) * ( ( ftyield * ftyield + fcyield * fcyield ) / ( fcyield * fcyield * ftyield ) );

        double T2 = mT * pow(yieldHardOne, 2.) * yieldHardTwo * ( sig / fcyield + rho * rFunction / ( sqrt(6.) * fcyield ) );

        double C3 =  sig / fcyield + rho / ( sqrt(6.) * fcyield );
        double C4 = -rho / ( sqrt(6.) * fcyield * fcyield ) - sig / ( fcyield * fcyield );

        double C2 = 2 * Al * ( 2 * ( 1 - yieldHardOne ) * C3 * C4 - sqrt(3. / 2.) * rho / ( fcyield * fcyield ) ) + mC * pow(yieldHardOne, 2.) * yieldHardTwo * ( sig / fcyield + rho * rFunction / ( sqrt(6.) * fcyield ) )
                    - pow(yieldHardOne, 2.) * yieldHardTwo * myield / ( fcyield * fcyield ) * ( sig  + rho * rFunction / sqrt(6.) );



        // compute dFDKappa
        dFDKappa = dFDYieldHardOne * dYieldHardOneDKappa + dFDYieldHardTwo * dYieldHardTwoDKappa
                   + ( T2 * T1 + C2 * C1 ) / deltaTime;



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
