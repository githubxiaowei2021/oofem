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

#ifndef concretedpm2ratePlastic_h
#define concretedpm2ratePlastic_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define _IFT_ConcreteDPM2RatePlastic_Name "con2dpmratePlastic"

#define _IFT_ConcreteDPM2RatePlastic_Ccompression "Ccompression"
#define _IFT_ConcreteDPM2RatePlastic_Kapparate0compression "Kapparate0compression"

#define _IFT_ConcreteDPM2RatePlastic_Ctension "Ctension"
#define _IFT_ConcreteDPM2RatePlastic_Kapparate0tension "Kapparate0tension"



namespace oofem {
/**
 * This class implements associated Material Status to ConcreteDPM2Rate.
 * Stores the characteristic length of the element.
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2RatePlasticStatus : public ConcreteDPM2Status
{
protected:

    double KappaRate = 0.;
    double tempKappaRate = 0.;
    
public:
    /// Constructor
    ConcreteDPM2RatePlasticStatus(GaussPoint *g);

    void updateYourself(TimeStep *tStep) override;
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    //double giveKappaRate() const {return KappaRate;}



    //double giveTempKappaRate() const { return tempKappaRate;}
    //void setTempKappaRate(double newStrainRate) { tempKappaRate = newStrainRate; }

    const char *giveClassName() const override { return "ConcreteDPM2RatePlasticStatus"; }

};



/**
 * This class implements rate dependence for ConcreteDPM2
 * @author: Xiaowei Liu, Peter Grassl
 */
class ConcreteDPM2RatePlastic : public ConcreteDPM2
{
 protected:
    double Ccompression;
    double Kapparate0compression;

    double Ctension;
    double Kapparate0tension;

  
public:
    /// Constructor
    ConcreteDPM2RatePlastic(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "ConcreteDPM2RatePlastic"; }




    double computeYieldValue(double sig, double rho, double theta, double tempKappa, double deltaTime, GaussPoint *gp) const;


    double computeDFDKappa(double sig, double rho, double theta, double tempKappa, double deltaTime, GaussPoint *gp) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;




protected:
};
} // end namespace oofem
#endif
