/*
 * OutputWriter.hpp
 *
 *  Created on: Sep 2, 2023
 *      Author: matej
 */

#ifndef OUTPUTWRITER_HPP_
#define OUTPUTWRITER_HPP_

#include "ElementsAndInput_ginac.hpp"
#include "TimeOrderingsAndTimeCuts_ginac.hpp"
#include "FeynmanAndAlphaParam.hpp"
#include "SectorDecomposition.hpp"

//===========================================================================
// code related to rewriting ginac expressions into wolfram mathematica syntax

std::string rewriteLogsIntoMathematicaFormat(GiNaC::ex _ex);

std::string rewriteGammaFunctionMathematicaFormat(GiNaC::ex _ex);
//===========================================================================

//===========================================================================
// code related to writing WLS output file - different NIntegrate methods
// coefficients from different FinalIntegrals are first summed and only then
// integrated - to be used with findDivergentPartPropToExtFreqWLS_together()
//
// NOTE: options of numerical method to be fixed inside functions

void writeWLSBasic(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeWLSLocalAdaptive(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeWLSMC(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeWLSCombinedMCAndEx(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeWLSForExactHighestPoleCoef(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

//===========================================================================

//===========================================================================
// code related to writing WLS output file - employing Vegas algorithm from
// cuba package through MathLink
// coefficients from different FinalIntegrals are first summed and only then
// integrated - to be used with findDivergentPartPropToExtFreqWLS_together()
//
// NOTE: options of numerical method to be fixed inside functions

void writeWLSVegas(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeWLSCombinedVegasAndEx(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeWLSCombinedVegasAndEx_separately(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);
//===========================================================================

//===========================================================================
// code related to writing WLS output file - different NIntegrate methods
// coefficients from different FinalIntegrals are integrated and only after that
// summed - to be used with findDivergentPartPropToExtFreqWLS()
//
// NOTE: options of numerical method to be fixed inside functions

void writeBegginingWLS(Diagram _diag, GiNaC::ex _overallFactor,
		std::string _path);
void writePartialWLSBasic(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		GiNaC::ex _poleCoef, int _poleCoefsSize, int _indexInPoleCoefs,
		std::string _path);
void writePartialWLSLocalAdaptive(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars, GiNaC::ex _poleCoef,
		int _poleCoefsSize, int _indexInPoleCoefs, std::string _path);
void writePartialWLSMC(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		GiNaC::ex _poleCoef, int _poleCoefsSize, int _indexInPoleCoefs,
		std::string _path);
void writeEndingWLS(Diagram _diag, std::string _path);

//===========================================================================

//===========================================================================
// Dividing the calculation into more parts and paralelization
void writeWLSCombinedExVegas_Nparts(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);
int countParts(std::string _poleCoef);
void divideTermIntoNParts(std::string _poleCoef, int _N,
		std::vector<std::string> &_vector);
std::vector<std::string> dividePoleCoefIntoNParts(std::string _poleCoef,
		int _N);

//===========================================================================

//===========================================================================
// Dividing the calculation into more parts C vegas
void writeCVegas_Eps0_Nparts(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path, int _parts);

int countParts_CVegas(std::string _poleCoef);

void divideTermIntoNParts_CVegas(std::string _poleCoef, int _N,
		std::vector<std::string> &_vector);

std::vector<std::string> dividePoleCoefIntoNParts_CVegas(std::string _poleCoef,
		int _N);

//===========================================================================

//===========================================================================
// code related to writing output file for Vegas implementation in C

void writeCombinedCVegasAndEx_3loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeCombinedCVegasAndEx_2loop(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeCombinedCVegasAndEx_1loop(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);


void writeCVegas_eps1_3loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeCVegas_eps0_3loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeCVegas_eps0_2loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path);

//===========================================================================

//===========================================================================
// code related to writing output file for Vegas implementation in C for finite parts
void writeCVegas_Finite_2loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _finitePartCoefs, GiNaC::ex _overallFactor,
		std::string _path, int _order);

void writeFiniteCVegas_eps0(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _finitePartCoefs, GiNaC::ex _overallFactor,
		std::string _path);

void writeFiniteCVegas_eps1(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _finitePartCoefs, GiNaC::ex _overallFactor,
		std::string _path);
//===========================================================================


#endif /* OUTPUTWRITER_HPP_ */
