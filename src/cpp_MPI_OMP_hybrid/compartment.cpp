#include <iostream>
#include <cmath>
#include <omp.h>
#include "parameters.h"
#include "utility.h"
#include "kernel.h"
#include "compartment.h"

using namespace std;
#define OPENMP_MASTER 0
#define INCLUDEBREAKAGE false
#define DUMP2D(varName) dump2DData(varName, #varName)
#define DUMP3D(varName) dump3DData(varName, #varName)
#define DUMP4D(varName) dump4DData(varName, #varName)

#define DUMP2DCSV(varName) dump2DCSV(varName, #varName)
#define DUMP3DCSV(varName) dump3DCSV(varName, #varName)
#define DUMP4DCSV(varName) dump4DCSV(varName, #varName)

CompartmentOut performCompartmentCalculations(PreviousCompartmentIn prevCompIn, CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double time, double timeStep)
{
    CompartmentOut compartmentOut;
    
    int s, ss, s1, ss1, s2, ss2, /*chunk,*/ tid = 0, /*nthreads,*/ a, b = 0;
    double totalSolidVolume = 0.0;
    double liquidAdditionRate = compartmentIn.liquidAdditionRate;

    #pragma omp parallel default(shared) private(s, ss, s1, ss1, s2, ss2, a, b, tid)
    {
    arrayOfInt2D sInd = vector2Array2D(compartmentIn.sInd);
    arrayOfInt2D ssInd = vector2Array2D(compartmentIn.ssInd);
    arrayOfInt2D sIndB = vector2Array2D(compartmentIn.sIndB);
    arrayOfInt2D ssIndB = vector2Array2D(compartmentIn.ssIndB);
    vector<double> vs = compartmentIn.vs;
    vector<double> vss = compartmentIn.vss;
    arrayOfDouble2D fAll = vector2Array2D(compartmentIn.fAll);
    arrayOfDouble2D fLiquid = vector2Array2D(compartmentIn.fLiquid);
    arrayOfDouble2D fGas = vector2Array2D(compartmentIn.fGas);
    arrayOfDouble2D sMeshXY = vector2Array2D(compartmentIn.sMeshXY);
    arrayOfDouble2D ssMeshXY = vector2Array2D(compartmentIn.ssMeshXY);
    arrayOfInt2D sAggregationCheck = vector2Array2D(compartmentIn.sAggregationCheck);
    arrayOfInt2D ssAggregationCheck = vector2Array2D(compartmentIn.ssAggregationCheck);
    arrayOfInt2D sCheckB = vector2Array2D(compartmentIn.sCheckB);
    arrayOfInt2D ssCheckB = vector2Array2D(compartmentIn.ssCheckB);
    arrayOfDouble2D sLow = vector2Array2D(compartmentIn.sLow);
    arrayOfDouble2D sHigh = vector2Array2D(compartmentIn.sHigh);
    arrayOfDouble2D ssLow = vector2Array2D(compartmentIn.ssLow);
    arrayOfDouble2D ssHigh = vector2Array2D(compartmentIn.ssHigh);
    //double liquidAdditionRate = compartmentIn.liquidAdditionRate;

    arrayOfDouble2D fAllComingIn = vector2Array2D(prevCompIn.fAllComingIn);
    arrayOfDouble2D fAllPreviousCompartment = vector2Array2D(prevCompIn.fAllPreviousCompartment);
    arrayOfDouble2D flPreviousCompartment = vector2Array2D(prevCompIn.flPreviousCompartment);
    arrayOfDouble2D fgComingIn = vector2Array2D(prevCompIn.fgComingIn);
    arrayOfDouble2D fgPreviousCompartment = vector2Array2D(prevCompIn.fgPreviousCompartment);

    //Declaration of vectors required initially for OMP implementation
    arrayOfDouble2D internalLiquid{{0.0}};
    arrayOfDouble2D externalLiquid{{0.0}};
    arrayOfDouble2D externalLiquidContent{{0.0}};
    arrayOfDouble2D externalVolumeBins{{0.0}};
    arrayOfDouble2D volumeBins{{0.0}};
    arrayOfDouble4D aggregationRate{{0.0}};
    arrayOfDouble2D depletionOfGasThroughAggregation{{0.0}};
    arrayOfDouble2D depletionOfLiquidThroughAggregation{{0.0}};
    arrayOfDouble2D birthThroughAggregation{{0.0}};
    arrayOfDouble2D firstSolidBirthThroughAggregation{{0.0}};
    arrayOfDouble2D secondSolidBirthThroughAggregation{{0.0}};
    arrayOfDouble2D liquidBirthThroughAggregation{{0.0}};
    arrayOfDouble2D gasBirthThroughAggregation{{0.0}};
    arrayOfDouble2D firstSolidVolumeThroughAggregation{{0.0}};
    arrayOfDouble2D secondSolidVolumeThroughAggregation{{0.0}};
    arrayOfDouble2D birthAggLowLow{{0.0}};
    arrayOfDouble2D birthAggHighHigh{{0.0}};
    arrayOfDouble2D birthAggLowHigh{{0.0}};
    arrayOfDouble2D birthAggHighLow{{0.0}};
    arrayOfDouble2D birthAggLowLowLiq{{0.0}};
    arrayOfDouble2D birthAggHighHighLiq{{0.0}};
    arrayOfDouble2D birthAggLowHighLiq{{0.0}};
    arrayOfDouble2D birthAggHighLowLiq{{0.0}};

    arrayOfDouble2D birthAggLowLowGas{{0.0}};
    arrayOfDouble2D birthAggHighHighGas{{0.0}};
    arrayOfDouble2D birthAggLowHighGas{{0.0}};
    arrayOfDouble2D birthAggHighLowGas{{0.0}};
    arrayOfDouble2D formationThroughAggregationCA{{0.0}};
    arrayOfDouble2D formationOfLiquidThroughAggregationCA{{0.0}};
    arrayOfDouble2D formationOfGasThroughAggregationCA{{0.0}};
    arrayOfDouble4D breakageRate{{0.0}};

    arrayOfDouble2D depletionThroughAggregation{{0.0}};
    arrayOfDouble2D depletionThroughBreakage{{0.0}};
    arrayOfDouble2D depletionOfGasThroughBreakage{{0.0}};
    arrayOfDouble2D depletionOfLiquidthroughBreakage{{0.0}};

    arrayOfDouble2D birthThroughBreakage1{{0.0}};
    arrayOfDouble2D birthThroughBreakage2{{0.0}};

    arrayOfDouble2D firstSolidBirthThroughBreakage{{0.0}};
    arrayOfDouble2D secondSolidBirthThroughBreakage{{0.0}};

    arrayOfDouble2D liquidBirthThroughBreakage1{{0.0}};
    arrayOfDouble2D gasBirthThroughBreakage1{{0.0}};

    arrayOfDouble2D liquidBirthThroughBreakage2{{0.0}};
    arrayOfDouble2D gasBirthThroughBreakage2{{0.0}};

    arrayOfDouble2D firstSolidVolumeThroughBreakage{{0.0}};
    arrayOfDouble2D secondSolidVolumeThroughBreakage{{0.0}};
    arrayOfDouble2D fractionBreakage00{{0.0}};
    arrayOfDouble2D fractionBreakage01{{0.0}};
    arrayOfDouble2D fractionBreakage10{{0.0}};
    arrayOfDouble2D fractionBreakage11{{0.0}};
    arrayOfDouble2D formationThroughBreakageCA{{0.0}};
    arrayOfDouble2D formationOfLiquidThroughBreakageCA{{0.0}};
    arrayOfDouble2D formationOfGasThroughBreakageCA{{0.0}};
    arrayOfDouble2D transferThroughLiquidAddition{{0.0}};
    arrayOfDouble2D transferThroughConsolidation{{0.0}};

    arrayOfDouble2D particleMovement{{0.0}};
    arrayOfDouble2D liquidMovement{{0.0}};
    arrayOfDouble2D gasMovement{{0.0}};

    //Calculation of liquid and gas bins
    //cout << "Begin liquidBins & gasBins" << endl;
    arrayOfDouble2D liquidBins{{0.0}};
    arrayOfDouble2D gasBins{{0.0}};
    arrayOfDouble2D internalVolumeBins{{0.0}};

    arrayOfDouble2D dfAlldt{{0.0}};
    arrayOfDouble2D dfLiquiddt{{0.0}};
    arrayOfDouble2D dfGasdt{{0.0}};
    
    arrayOfDouble2D meshXYSum{{0.0}};

    arrayOfDouble4D aggregationKernel{{0.0}};
    arrayOfDouble4D breakageKernel{{0.0}};// = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);;

    //double totalSolidVolume = 0.0;
    /*int s, ss, s1, ss1, s2, ss2, /*chunk, tid = 0, /*nthreads, a, b = 0;

    #pragma omp parallel default(shared) private(s, ss, s1, ss1, s2, ss2, a, b, tid)
    { */
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            //nthreads = omp_get_num_threads();
            //chunk = static_cast<int>(NUMBEROFFIRSTSOLIDBINS / nthreads);
        }
        //cout << "Before Schedule: I am omp_id = " << tid << endl;   
        
        #pragma omp for 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            {
                double temp = fAll[s][ss];
                if (fabs(temp) > EPSILON) //If there are particles in the bin...
                {
                    //New liquid bins are calculated as total amount liquid in that size class divided by
                    //the number of particle in that size class
                    liquidBins[s][ss] = fLiquid[s][ss] / temp;

                    // New gas bins are calculated as total amount liquid in that size class divided by
                    //the number of particle in that size class
                    gasBins[s][ss] = fGas[s][ss] / temp;
                }
                else
                {
                    liquidBins[s][ss] = 0.0; // No particles in the bin -> particles in that size class has no liquid
                    gasBins[s][ss] = 0.0;    // No particles in the bin -> particles in that size class has no gas
                }
            }
        //cout << "End liquidBins & gasBins" << endl;

        //Internal and external liquid demarcation
        //cout << "Begin Internal & External liquid" << endl;
    
        #pragma omp for 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
    	{
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            {
                double temp = liquidBins[s][ss];
                double temp1 = min(GRANULESATURATIONFACTOR * gasBins[s][ss], temp);
                double temp2 = max(0.0, temp - internalLiquid[s][ss]);
                //internalLiquid[s][ss] = min(GRANULESATURATIONFACTOR * gasBins[s][ss], liquidBins[s][ss]);
                //externalLiquid[s][ss] = max(0.0, liquidBins[s][ss] - internalLiquid[s][ss]);
                //#pragma omp critical
                {
                    externalLiquidContent[s][ss] = temp2 / temp;
                    internalLiquid[s][ss] = temp1;
                    externalLiquid[s][ss] = temp2;
                }
            }
        }
        //cout << "End Internal & External liquid" << endl;

        //cout << "Begin Internal, External & Volume Bins" << endl;

        //compartmentOut.internalVolumeBins = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        #pragma omp for 
        //if (tid == OPENMP_MASTER)
        //{
            for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
            {
                for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
                {
                    double t1 = sMeshXY[s][ss]; double t2 = ssMeshXY[s][ss]; double t3 = gasBins[s][ss];
                    double temp1 = t1 + t2 + internalLiquid[s][ss] + t3;
                    double temp2 = t1 + t2 + liquidBins[s][ss] + t3;
                    double temp3 = t1 + t2;

                    //#pragma omp critical
                    {
                        internalVolumeBins[s][ss] = temp1;
                        //internalVolumeBins[s][ss] += internalLiquid[s][ss] + gasBins[s][ss];
                        externalVolumeBins[s][ss] = temp2;
                        //externalVolumeBins[s][ss] += liquidBins[s][ss] + gasBins[s][ss];
                        volumeBins[s][ss] = temp3;
                    }
                }
            }
    
        //cout << "End Internal, External & Volume Bins" << endl;

        //cout << "Begin Aggregation Kernel" << endl;
        //    compartmentOut.aggregationKernel = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS,
        //                                       NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        //    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        //        for (int ss1  = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
        //            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
        //                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
        //                {
        //                    double expr1 = externalVolumeBins[s1][ss1]+externalVolumeBins[s2][ss2]-externalLiquid[s1][ss1]-externalLiquid[s2][ss2];
        //                    expr1 = pow(expr1, AGGREGATIONKERNELCONSTANTGAMMA);
        //
        //                    double expr2 = (externalLiquid[s1][ss1]/externalVolumeBins[s1][ss1] + externalLiquid[s1][ss1]/externalVolumeBins[s1][ss1])/2.0;
        //                    expr2 = pow(expr2, AGGREGATIONKERNELCONSTANTALPHA);
        //
        //                    double expr3 = 1-(externalLiquid[s1][ss1]/externalVolumeBins[s1][ss1] + externalLiquid[s2][ss2]/externalVolumeBins[s2][ss2])/2.0;
        //                    expr3 = pow(expr3, AGGREGATIONKERNELCONSTANTDELTA);
        //
        //                    compartmentOut.aggregationKernel[s1][ss1][s2][ss2] = AGGREGATIONKERNELCONSTANT*expr1*pow(expr2*expr3, AGGREGATIONKERNELCONSTANTALPHA);
        //                    // compartmentOut.aggregationKernel[s1][ss1][s2][ss2] = DEMAGGREGATIONKERNELCONST*DEMAGGREGATIONKERNELVALUE;
        //                }
        //compartmentOut.aggregationKernel = DEMDependentAggregationKernel(compartmentIn, compartmentDEMIn, externalLiquidContent, timeStep);
        aggregationKernel = DEMDependentAggregationKernel(compartmentIn, compartmentDEMIn, externalLiquidContent, timeStep);
        #pragma omp barrier
        //cout << "End Aggregation Kernel" << endl;

        //AGGREGATION
        //cout << "Begin aggregationRate" << endl;
        #pragma omp for 
        for (s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
            for (ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
                for (s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                {
                    //#pragma omp ordered simd safelen(8)
                    for (ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    {    
                        double temp1 = sAggregationCheck[s1][ss1] * ssAggregationCheck[s2][ss2] * aggregationKernel[s1][ss1][s2][ss2] * fAll[s1][ss1] * fAll[s2][ss2];
                        //#pragma omp critical
                        aggregationRate[s1][ss1][s2][ss2] = temp1;
                        depletionThroughAggregation[s1][ss1] += temp1;
                        depletionThroughAggregation[s2][ss2] += temp1;
                    
                    }
                }

        //cout << "End aggregationRate" << endl;

    
        //cout << "Begin depletionThroughAggregation" << endl;
        //#pragma omp for schedule(static, chunk)
        /*if (tid == 0) 
        {
            for (s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
                for (ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
                    for (s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                        for (ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                        {
                            depletionThroughAggregation[s1][ss1] += aggregationRate[s1][ss1][s2][ss2];
                            depletionThroughAggregation[s2][ss2] += aggregationRate[s1][ss1][s2][ss2];
                        }
        }*/
    
        //cout << "End depletionThroughAggregation" << endl;

        //cout << "Begin depletionOfGasThroughAggregation & depletionOfLiquidThroughAggregation" << endl;

        #pragma omp for 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
        {
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            {
                double temp3 = depletionThroughAggregation[s][ss];
                double temp1 = temp3 * gasBins[s][ss];
                double temp2 = temp3 * liquidBins[s][ss];
                //#pragma omp critical
                {
                    depletionOfGasThroughAggregation[s][ss] = temp1;
                    depletionOfLiquidThroughAggregation[s][ss] = temp2;
                }
            }

        }
        //cout << "End depletionOfGasThroughAggregation & depletionOfLiquidThroughAggregation" << endl << endl;

        //cout << "Begin birthThroughAggregation, firstSolidBirthThroughAggregation, secondSolidBirthThroughAggregation, ";
        //cout << "liquidBirthThroughAggregation & gasBirthThroughAggregation" << endl;
    
        //#pragma omp for
        if (tid == OPENMP_MASTER)
        {
            for (s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
                for (ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
                    for (s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                        for (ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                        {
                            for (a = 0; a < NUMBEROFFIRSTSOLIDBINS; a++)
                            {
                                for (b = 0; b < NUMBEROFSECONDSOLIDBINS; b++)
                                {
                                    if (sInd[s1][s2] == (a + 1) && ssInd[ss1][ss2] == (b + 1))
                                    {
                                        #pragma omp atomic update
                                        birthThroughAggregation[a][b] += aggregationRate[s1][ss1][s2][ss2];
                                        #pragma omp atomic update
                                        firstSolidBirthThroughAggregation[a][b] += (vs[s1] + vs[s2]) * aggregationRate[s1][ss1][s2][ss2];
                                        #pragma omp atomic update
                                        secondSolidBirthThroughAggregation[a][b] += (vss[ss1] + vss[ss2]) * aggregationRate[s1][ss1][s2][ss2];
                                        #pragma omp atomic update
                                        liquidBirthThroughAggregation[a][b] += (liquidBins[s1][ss1] + liquidBins[s2][ss2]) * aggregationRate[s1][ss1][s2][ss2];
                                        #pragma omp atomic update
                                        gasBirthThroughAggregation[a][b] += (gasBins[s1][ss1] + gasBins[s2][ss2]) * aggregationRate[s1][ss1][s2][ss2];
                                    }
                                }
                            }
                        }
        }
        #pragma omp barrier
        //cout << "End birthThroughAggregation, firstSolidBirthThroughAggregation, secondSolidBirthThroughAggregation, ";
        //cout << "liquidBirthThroughAggregation & gasBirthThroughAggregation" << endl << endl;

        //cout << "Begin firstSolidVolumeThroughAggregation & secondSolidVolumeThroughAggregation" << endl;
        //#pragma omp barrier
        #pragma omp for 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            {
                double temp1 = birthThroughAggregation[s][ss];
                if (fabs(temp1) > EPSILON)
                {
                    firstSolidVolumeThroughAggregation[s][ss] = firstSolidBirthThroughAggregation[s][ss] / temp1;
                    secondSolidVolumeThroughAggregation[s][ss] = secondSolidBirthThroughAggregation[s][ss] / temp1;
                }
                
            }
        //cout << "After VolBins: I am omp_id = " << tid << endl;   
    
    
        //cout << "End firstSolidVolumeThroughAggregation & secondSolidVolumeThroughAggregation" << endl << endl;

        //cout << "Begin birth_agg_low_low, birth_agg_high_high, birth_agg_low_high & birth_agg_high_low" << endl;

        #pragma omp for 
        //if (tid == 0)
        //{
            for (s = 0; s < NUMBEROFFIRSTSOLIDBINS - 1; s++)
            {    for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS - 1; ss++)
                {
                    double tvs1 = vs[s + 1]; double tvs = vs[s]; double tvss1 = vss[ss + 1]; double tvss = vss[ss];
                    double temp_first = firstSolidVolumeThroughAggregation[s][ss]; double temp_second = secondSolidVolumeThroughAggregation[s][ss]; 
                    double temp_birthagg = birthThroughAggregation[s][ss]; double temp_liqbirth = liquidBirthThroughAggregation[s][ss];
                    double temp_gasagg = gasBirthThroughAggregation[s][ss];
                	#pragma omp atomic write
                    birthAggLowLow[s][ss] = ((tvs1 - temp_first) / (tvs1 - tvs)) * \
                    ((tvss1 - temp_second) / (tvss1 - tvss)) * temp_birthagg;
                    //#pragma omp atomic update
                    //birthAggLowLow[s][ss] *= ((vss[ss + 1] - secondSolidVolumeThroughAggregation[s][ss]) / (vss[ss + 1] - vss[ss]));
                    //birthAggLowLow[s][ss] *= birthThroughAggregation[s][ss];
                	
                	#pragma omp atomic write
                    birthAggHighHigh[s + 1][ss + 1] = ((temp_first - tvs) / (tvs1 - tvs)) * \
                    ((temp_second - tvss) / (tvss1 - tvss)) * temp_birthagg;
                    //birthAggHighHigh[s + 1][ss + 1] *= (secondSolidVolumeThroughAggregation[s][ss] - vss[ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggHighHigh[s + 1][ss + 1] *= birthThroughAggregation[s][ss];
                	
                	#pragma omp atomic write
                    birthAggLowHigh[s][ss + 1] = ((tvs1 - temp_first) / (tvs1 - tvs)) * \
                    ((temp_second - tvss) / (tvss1 - tvss)) * temp_birthagg;
                    //birthAggLowHigh[s][ss + 1] *= (secondSolidVolumeThroughAggregation[s][ss] - vss[ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggLowHigh[s][ss + 1] *= birthThroughAggregation[s][ss];

                	#pragma omp atomic write
                    birthAggHighLow[s + 1][ss] = ((temp_first - tvs) / (tvs1 - tvs)) * \
                    ((tvss1 - temp_second) / (tvss1 - tvss)) * temp_birthagg;
                    //birthAggHighLow[s + 1][ss] *= (vss[ss + 1] - secondSolidVolumeThroughAggregation[s][ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggHighLow[s + 1][ss] *= birthThroughAggregation[s][ss];

            //cout << "End birth_agg_low_low, birth_agg_high_high, birth_agg_low_high & birth_agg_high_low" << endl << endl;

            //cout << "Begin birth_agg_low_low_liq, birth_agg_high_high_liq, birth_agg_low_high_liq & birth_agg_high_low_liq" << endl;

                	#pragma omp atomic write
                    birthAggLowLowLiq[s][ss] = ((tvs1 - temp_first) / (tvs1 - tvs)) * \
                    ((tvss1 - temp_second) / (tvss1 - tvss)) * temp_liqbirth;
                    //birthAggLowLowLiq[s][ss] *= (vss[ss + 1] - secondSolidVolumeThroughAggregation[s][ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggLowLowLiq[s][ss] *= liquidBirthThroughAggregation[s][ss];

                	#pragma omp atomic write
                    birthAggHighHighLiq[s + 1][ss + 1] = ((temp_first - tvs) / (tvs1 - tvs)) * \
                    ((temp_second - tvss) / (tvss1 - tvss)) * temp_liqbirth;
                    //birthAggHighHighLiq[s + 1][ss + 1] *= (secondSolidVolumeThroughAggregation[s][ss] - vss[ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggHighHighLiq[s + 1][ss + 1] *= liquidBirthThroughAggregation[s][ss];

                	#pragma omp atomic write
                    birthAggLowHighLiq[s][ss + 1] = ((tvs1 - temp_first) / (tvs1 - tvs)) * \
                    ((temp_second - tvss) / (tvss1 - tvss)) * temp_liqbirth;
                    //birthAggLowHighLiq[s][ss + 1] *= (secondSolidVolumeThroughAggregation[s][ss] - vss[ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggLowHighLiq[s][ss + 1] *= liquidBirthThroughAggregation[s][ss];

                	#pragma omp atomic write
                    birthAggHighLowLiq[s + 1][ss] = ((temp_first - tvs) / (tvs1 - tvs)) * \
                    ((tvss1 - temp_second) / (tvss1 - tvss)) * temp_liqbirth;
                    //birthAggHighLowLiq[s + 1][ss] *= (vss[ss + 1] - secondSolidVolumeThroughAggregation[s][ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggHighLowLiq[s + 1][ss] *= liquidBirthThroughAggregation[s][ss];
            //cout << "End birth_agg_low_low_liq, birth_agg_high_high_liq, birth_agg_low_high_liq & birth_agg_high_low_liq" << endl << endl;

            //cout << "Begin birth_agg_low_low_gas, birth_agg_high_high_gas, birth_agg_low_high_gas & birth_agg_high_low_gas" << endl;

                	#pragma omp atomic write
                    birthAggLowLowGas[s][ss] = ((tvs1 - temp_first) / (tvs1 - tvs)) * \
                    ((tvss1 - temp_second) / (tvss1 - tvss)) * temp_gasagg;
                    //birthAggLowLowGas[s][ss] *= (vss[ss + 1] - secondSolidVolumeThroughAggregation[s][ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggLowLowGas[s][ss] *= gasBirthThroughAggregation[s][ss];

                	#pragma omp atomic write
                    birthAggHighHighGas[s + 1][ss + 1] = ((temp_first - tvs) / (tvs1 - tvs)) * \
                    ((temp_second - tvss) / (tvss1 - tvss)) * temp_gasagg;
                    //birthAggHighHighGas[s + 1][ss + 1] *= (secondSolidVolumeThroughAggregation[s][ss] - vss[ss]) / (tvss1 - vss[ss]);
                    //birthAggHighHighGas[s + 1][ss + 1] *= gasBirthThroughAggregation[s][ss];

                	#pragma omp atomic write
                    birthAggLowHighGas[s][ss + 1] = ((tvs1 - temp_first) / (tvs1 - tvs)) * \
                    ((temp_second - tvss) / (tvss1 - tvss)) * temp_gasagg;
                    //birthAggLowHighGas[s][ss + 1] *= (secondSolidVolumeThroughAggregation[s][ss] - vss[ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggLowHighGas[s][ss + 1] *= gasBirthThroughAggregation[s][ss];

                    #pragma omp atomic write
                    birthAggHighLowGas[s + 1][ss] = ((temp_first - tvs) / (tvs1 - tvs)) * \
                    ((tvss1 - temp_second) / (tvss1 - tvss)) * temp_gasagg;
                    //birthAggHighLowGas[s + 1][ss] *= (vss[ss + 1] - secondSolidVolumeThroughAggregation[s][ss]) / (vss[ss + 1] - vss[ss]);
                    //birthAggHighLowGas[s + 1][ss] *= gasBirthThroughAggregation[s][ss];
                }
            }
    
        //}
        //#pragma omp barrier
        //cout << "End birth_agg_low_low_gas, birth_agg_high_high_gas, birth_agg_low_high_gas & birth_agg_high_low_gas" << endl << endl;

        //cout << "Begin formationThroughAggregationCA, formationOfLiquidThroughAggregationCA & formationOfGasThroughAggregationCA" << endl;
        #pragma omp for 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            {
                formationThroughAggregationCA[s][ss] = birthAggLowLow[s][ss] + birthAggHighHigh[s][ss] + birthAggLowHigh[s][ss] + birthAggHighLow[s][ss];
                formationOfLiquidThroughAggregationCA[s][ss] = birthAggLowLowLiq[s][ss] + birthAggHighHighLiq[s][ss] + birthAggLowHighLiq[s][ss] + birthAggHighLowLiq[s][ss];
                formationOfGasThroughAggregationCA[s][ss] = birthAggLowLowGas[s][ss] + birthAggHighHighGas[s][ss] + birthAggLowHighGas[s][ss] + birthAggHighLowGas[s][ss];
            }
    
        //cout << "End formationThroughAggregationCA, formationOfLiquidThroughAggregationCA & formationOfGasThroughAggregationCA" << endl << endl;

        //cout << "************End of Aggregation**************" << endl << endl;

        //BREAKAGE
        //Breakage kernel Calculation inside time loop becasue gas and liquid bins change with time
        //cout << "Begin Breakage Kernel" << endl;
        //compartmentOut.breakageKernel = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);//, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        //compartmentOut.breakageKernel = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        //    for (int s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
        //        for (int ss  = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
        //        {
        //            double expr1 = sqrt(4/(15*M_PI));
        //            double expr2 = pow(SHEARRATE, 2)* cbrt(3*externalVolumeBins[s][ss]/(4*M_PI));
        //            compartmentOut.breakageKernel[s][ss] = expr1*SHEARRATE*exp(-BREAKAGEKERNELCONSTANT/expr2);
        //            //compartmentOut.breakageKernel[s][ss] = DEMBREAKAGEKERNELCONST * DEMBREAKAGEKERNELVALUE;
        //        }
        //    //reset to zero values
        //    compartmentOut.breakageKernel = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        //compartmentOut.breakageKernel = DEMDependentBreakageKernel(compartmentIn, compartmentDEMIn, timeStep);
        if (tid == OPENMP_MASTER)
        {
            if (INCLUDEBREAKAGE)
                breakageKernel = DEMDependentBreakageKernel(compartmentIn, compartmentDEMIn, timeStep);
            //else
                //breakageKernel = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

            //cout << "End Breakage Kernel" << endl;

            //cout << "Begin breakageRate, depletionThroughBreakage, depletionOfLiquidthroughBreakage & depletionOfGasThroughBreakage" << endl;
            if (INCLUDEBREAKAGE)
            {
                for (s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
                    for (ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
                        for (s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                            for (ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                            {
                                breakageRate[s1][ss1][s2][ss2] = sCheckB[s1][s2] * ssCheckB[ss1][ss2] * breakageKernel[s1][ss1][s2][ss2] * fAll[s1][ss1];
                                //breakageRate[s1][ss1][s2][ss2] *= breakageKernel[s1][ss1][s2][ss2] * fAll[s1][ss1];
                                depletionThroughBreakage[s1][ss1] += breakageRate[s1][ss1][s2][ss2];
                                depletionOfLiquidthroughBreakage[s1][ss1] = depletionThroughBreakage[s1][ss1] * liquidBins[s1][ss1];
                                depletionOfGasThroughBreakage[s1][ss1] = depletionThroughBreakage[s1][ss1] * gasBins[s1][ss1];

                                birthThroughBreakage1[s1][ss1] += breakageRate[s1][ss1][s2][ss2];
                            }
            }
            //cout << "End breakageRate, depletionThroughBreakage, depletionOfLiquidthroughBreakage & depletionOfGasThroughBreakage" << endl << endl;
            //cout << "Begin birthThroughBreakage2, firstSolidBirthThroughBreakage, secondSolidBirthThroughBreakage, ";
            //cout << "liquidBirthThroughBreakage1, gasBirthThroughBreakage1, liquidBirthThroughBreakage2, ";
            //cout << "gasBirthThroughBreakage2, firstSolidVolumeThroughBreakage & secondSolidVolumeThroughBreakage" << endl;
            if (INCLUDEBREAKAGE)
            {
                for (s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
                    for (ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
                        for (s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                            for (ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                            {
                                for (a = 0; a < NUMBEROFFIRSTSOLIDBINS; a++)
                                    for (b = 0; b < NUMBEROFSECONDSOLIDBINS; b++)
                                    {
                                        if (sIndB[s1][s2] == (a + 1) && ssIndB[ss1][ss2] == (b + 1))
                                        {
                                            birthThroughBreakage2[a][b] += breakageRate[s1][ss1][s2][ss2];

                                            firstSolidBirthThroughBreakage[a][b] += (vs[s1] - vs[s2]) * breakageRate[s1][ss1][s2][ss2];
                                            secondSolidBirthThroughBreakage[a][b] += (vss[ss1] - vss[ss2]) * breakageRate[s1][ss1][s2][ss2];

                                            liquidBirthThroughBreakage2[a][b] += (liquidBins[s1][ss1] * (1 - (volumeBins[s2][ss2] / volumeBins[s1][ss1]))) * breakageRate[s1][ss1][s2][ss2];
                                            gasBirthThroughBreakage2[a][b] += (gasBins[s1][ss1] * (1 - (volumeBins[s2][ss2] / volumeBins[s1][ss1]))) * breakageRate[s1][ss1][s2][ss2];
                                            if (fabs(birthThroughBreakage2[a][b]) > 1e-16)
                                            {
                                                firstSolidVolumeThroughBreakage[a][b] = firstSolidBirthThroughBreakage[a][b] / birthThroughBreakage2[a][b];
                                                secondSolidVolumeThroughBreakage[a][b] = secondSolidBirthThroughBreakage[a][b] / birthThroughBreakage2[a][b];
                                            }
                                        }
                                    }
                                liquidBirthThroughBreakage1[s2][ss2] += (liquidBins[s1][ss1] * (volumeBins[s2][ss2] / volumeBins[s1][ss1])) * breakageRate[s1][ss1][s2][ss2];
                                gasBirthThroughBreakage1[s2][ss2] += (gasBins[s1][ss1] * (volumeBins[s2][ss2] / volumeBins[s1][ss1])) * breakageRate[s1][ss1][s2][ss2];
                            }
            }
            //cout << "End birthThroughBreakage2, firstSolidBirthThroughBreakage, secondSolidBirthThroughBreakage, ";
            //cout << "liquidBirthThroughBreakage1, gasBirthThroughBreakage1, liquidBirthThroughBreakage2, ";
            //cout << "gasBirthThroughBreakage2, firstSolidVolumeThroughBreakage & secondSolidVolumeThroughBreakage" << endl << endl;

            //cout << "Begin fractionBreakage00, fractionBreakage01, fractionBreakage10 & fractionBreakage11" << endl;

            if (INCLUDEBREAKAGE)
            {
                for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
                {
                    double value1 = 0.0;
                    double value2 = 0.0;
                    for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
                    {
                        value1 = fabs(sLow[s][ss] - firstSolidVolumeThroughBreakage[s][ss]);
                        value1 = sHigh[s][ss] - sLow[s][ss] - value1;
                        value1 /= (sHigh[s][ss] - sLow[s][ss]);
                        value2 = fabs(ssLow[s][ss] - secondSolidVolumeThroughBreakage[s][ss]);
                        value2 = ssHigh[s][ss] - ssLow[s][ss] - value2;
                        value2 /= (ssHigh[s][ss] - ssLow[s][ss]);
                        fractionBreakage00[s][ss] = value1 / value2;

                        value2 = fabs(ssHigh[s][ss] - secondSolidVolumeThroughBreakage[s][ss]);
                        value2 = ssHigh[s][ss] - ssLow[s][ss] - value2;
                        value2 /= (ssHigh[s][ss] - ssLow[s][ss]);
                        fractionBreakage01[s][ss] = value1 / value2;

                        value1 = fabs(sHigh[s][ss] - firstSolidVolumeThroughBreakage[s][ss]);
                        value1 = sHigh[s][ss] - sLow[s][ss] - value1;
                        value1 /= (sHigh[s][ss] - sLow[s][ss]);
                        fractionBreakage11[s][ss] = value1 / value2;

                        value2 = fabs(ssLow[s][ss] - secondSolidVolumeThroughBreakage[s][ss]);
                        value2 = ssHigh[s][ss] - ssLow[s][ss] - value2;
                        value2 /= (ssHigh[s][ss] - ssLow[s][ss]);
                        fractionBreakage10[s][ss] = value1 / value2;
                    }
                }
            }

            //cout << "End fractionBreakage00, fractionBreakage01, fractionBreakage10 & fractionBreakage11" << endl << endl;

            //cout << "Begin formationThroughBreakageCA, formationOfLiquidThroughBreakageCA & formationOfGasThroughBreakageCA" << endl;
            if (INCLUDEBREAKAGE)
            {
                for (s = 0; s < NUMBEROFFIRSTSOLIDBINS - 1; s++)
                    for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS - 1; ss++)
                    {
                        formationThroughBreakageCA[s][ss] += birthThroughBreakage2[s][ss] * fractionBreakage00[s][ss];
                        formationThroughBreakageCA[s][ss + 1] += birthThroughBreakage2[s][ss] * fractionBreakage01[s][ss];
                        formationThroughBreakageCA[s + 1][ss] += birthThroughBreakage2[s][ss] * fractionBreakage10[s][ss];
                        formationThroughBreakageCA[s + 1][ss + 1] += birthThroughBreakage2[s][ss] * fractionBreakage11[s][ss];

                        formationOfLiquidThroughBreakageCA[s][ss] += liquidBirthThroughBreakage2[s][ss] * fractionBreakage00[s][ss];
                        formationOfLiquidThroughBreakageCA[s][ss + 1] += liquidBirthThroughBreakage2[s][ss] * fractionBreakage01[s][ss];
                        formationOfLiquidThroughBreakageCA[s + 1][ss] += liquidBirthThroughBreakage2[s][ss] * fractionBreakage10[s][ss];
                        formationOfLiquidThroughBreakageCA[s + 1][ss + 1] += liquidBirthThroughBreakage2[s][ss] * fractionBreakage11[s][ss];

                        formationOfGasThroughBreakageCA[s][ss] += gasBirthThroughBreakage2[s][ss] * fractionBreakage00[s][ss];
                        formationOfGasThroughBreakageCA[s][ss + 1] += gasBirthThroughBreakage2[s][ss] * fractionBreakage01[s][ss];
                        formationOfGasThroughBreakageCA[s + 1][ss] += gasBirthThroughBreakage2[s][ss] * fractionBreakage10[s][ss];
                        formationOfGasThroughBreakageCA[s + 1][ss + 1] += gasBirthThroughBreakage2[s][ss] * fractionBreakage11[s][ss];
                    }
            }
        }

        //cout << "End formationThroughBreakageCA, formationOfLiquidThroughBreakageCA & formationOfGasThroughBreakageCA" << endl << endl;

        //cout << "**************End of Breakage****************" << endl << endl;

        //cout << "Begin Particle Transfer" << endl;
        //cout << "Finding Max of s_meshxy+ss_meshxy " << endl;

        
        #pragma omp for 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
                meshXYSum[s][ss] = sMeshXY[s][ss] + ssMeshXY[s][ss];
        
        double maxMeshXY = getMaximumOf2DArray(meshXYSum);
        double value = PARTICLEAVERAGEVELOCITY * timeStep / DISTANCEBETWEENCOMPARTMENTS;
        
        #pragma omp for 
        //if (tid == 0)
        //{    
            for (s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
                for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
                {
                    double valueMeshXY = 1 - (sMeshXY[s][ss] + ssMeshXY[s][ss]) / maxMeshXY;

                    particleMovement[s][ss] = fAllComingIn[s][ss] + fAllPreviousCompartment[s][ss] * value * valueMeshXY - fAll[s][ss] * value;
                    //particleMovement[s][ss] += fAllPreviousCompartment[s][ss] * value * valueMeshXY;
                    //particleMovement[s][ss] -= fAll[s][ss] * value;

                    liquidMovement[s][ss] = flPreviousCompartment[s][ss] * value * valueMeshXY - fLiquid[s][ss] * value;
                    //liquidMovement[s][ss] -= fLiquid[s][ss] * value;

                    gasMovement[s][ss] = fgComingIn[s][ss] + fgPreviousCompartment[s][ss] * value * valueMeshXY - fGas[s][ss] * value;
                    //gasMovement[s][ss] += fgPreviousCompartment[s][ss] * value * valueMeshXY;
                    //gasMovement[s][ss] -= fGas[s][ss] * value;
                }
       
            //cout << "End Particle Transfer" << endl;
        //}
        //cout << "Begin rate calculations" << endl << endl;
        if (time >= PREMIXINGTIME && time <= PREMIXINGTIME + LIQUIDADDITIONTIME)
            liquidAdditionRate *= timeStep;
        else
            liquidAdditionRate = 0.0;

        #pragma omp for reduction(+:totalSolidVolume) 
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS - 1; s++)
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS - 1; ss++)
                totalSolidVolume = totalSolidVolume + fAll[s][ss] * (vs[s] + vss[ss]);
    
        //compartmentOut.dfAlldt = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        //compartmentOut.dfLiquiddt = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
        //compartmentOut.dfGasdt = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

        #pragma omp for 
        //if (tid == 0)
        //{
        for (s = 0; s < NUMBEROFFIRSTSOLIDBINS - 1; s++)
            for (ss = 0; ss < NUMBEROFSECONDSOLIDBINS - 1; ss++)
            {
               	#pragma omp atomic write
                dfAlldt[s][ss] = particleMovement[s][ss] + formationThroughAggregationCA[s][ss] - depletionThroughAggregation[s][ss] + birthThroughBreakage1[s][ss] \
                + formationThroughBreakageCA[s][ss] - depletionThroughBreakage[s][ss];
                //dfAlldt[s][ss] += formationThroughAggregationCA[s][ss] - depletionThroughAggregation[s][ss];
                //dfAlldt[s][ss] += birthThroughBreakage1[s][ss] + formationThroughBreakageCA[s][ss] - depletionThroughBreakage[s][ss];

                if (totalSolidVolume > EPSILON)
                {
                    double value = (vs[s] + vss[ss]) / totalSolidVolume;
                    transferThroughLiquidAddition[s][ss] = liquidAdditionRate * value;
                }

                #pragma omp atomic write
                dfLiquiddt[s][ss] = liquidMovement[s][ss] + fAll[s][ss] * transferThroughLiquidAddition[s][ss] + formationOfLiquidThroughAggregationCA[s][ss] - depletionOfLiquidThroughAggregation[s][ss] \
                + liquidBirthThroughBreakage1[s][ss] + formationOfLiquidThroughBreakageCA[s][ss] - depletionOfLiquidthroughBreakage[s][ss];
                //dfLiquiddt[s][ss] += fAll[s][ss] * transferThroughLiquidAddition[s][ss];
                //dfLiquiddt[s][ss] += formationOfLiquidThroughAggregationCA[s][ss] - depletionOfLiquidThroughAggregation[s][ss];
                //dfLiquiddt[s][ss] += liquidBirthThroughBreakage1[s][ss] + formationOfLiquidThroughBreakageCA[s][ss];
                //dfLiquiddt[s][ss] -= depletionOfLiquidthroughBreakage[s][ss];

               	#pragma omp atomic write
                dfGasdt[s][ss] = gasMovement[s][ss] + fAll[s][ss] * transferThroughConsolidation[s][ss] + formationOfGasThroughAggregationCA[s][ss] - depletionOfGasThroughAggregation[s][ss] \
                + gasBirthThroughBreakage1[s][ss] + formationOfGasThroughBreakageCA[s][ss] - depletionOfGasThroughBreakage[s][ss];
                //dfGasdt[s][ss] += fAll[s][ss] * transferThroughConsolidation[s][ss];
                //dfGasdt[s][ss] += formationOfGasThroughAggregationCA[s][ss] - depletionOfGasThroughAggregation[s][ss];
                //dfGasdt[s][ss] += gasBirthThroughBreakage1[s][ss] + formationOfGasThroughBreakageCA[s][ss];
                //dfGasdt[s][ss] -= depletionOfGasThroughBreakage[s][ss];
            }
    
        //}
        //cout << "**************End of Rate Calculations****************" << endl << endl;

        //cout << "Return to Model code" << endl << endl;

    //#pragma omp master 
    if (tid == OPENMP_MASTER)
    {
    
	    compartmentOut.liquidBins = array2Vector2D(liquidBins);
	    compartmentOut.gasBins = array2Vector2D(gasBins);
	    compartmentOut.internalVolumeBins = array2Vector2D(internalVolumeBins);
	    //compartmentOut.aggregationKernel = aggregationKernel;
	    //compartmentOut.breakageKernel = breakageKernel;
	    compartmentOut.dfAlldt = array2Vector2D(dfAlldt);
	    compartmentOut.dfLiquiddt = array2Vector2D(dfLiquiddt);
	    compartmentOut.dfGasdt = array2Vector2D(dfGasdt);
	}
    //#pragma omp barrier
    } // omp parallel closed 
    return compartmentOut;
}
