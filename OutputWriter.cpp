/*
 * OutputWriter.cpp
 *
 *  Created on: Sep 2, 2023
 *      Author: matej
 */

#include <sstream>
#include <iostream>
#include "OutputWriter.hpp"

//===========================================================================
// code related to rewriting ginac expressions into wolfram mathematica syntax

std::string rewriteLogsIntoMathematicaFormat(GiNaC::ex _ex) {
	// the point is ginac uses log() while mathematica requires Log[]
	// this function rewrites log() -> Log[] and returns string that
	// can be then given to wolfram mathematica

	std::string result;
	std::stringstream a;
	a << _ex;
	std::string input = a.str();

	char current;
	bool startOfLog = false;
	bool parenthesesStarted = false;
	int innerParenthesisCount = 0;

	// create scanner on input expression (in string format)
	std::istringstream scanner(input);
	for (int i = 0; i < input.size(); i++) {

		current = scanner.get();

		if (current == 'l') {
			startOfLog = true;
			result.push_back('L');
			continue;
		}

		if (startOfLog) {
			// if startOfLog = true and you find '(' or ')'

			if (current == '(' && !parenthesesStarted) {
				result.push_back('[');
				parenthesesStarted = true;
				continue;
			}

			if (current == '(' && parenthesesStarted) {
				result.push_back('(');
				innerParenthesisCount++;
				continue;
			}

			if (current == ')' && innerParenthesisCount != 0) {
				result.push_back(')');
				innerParenthesisCount--;
				//startOfLog = false;
				continue;
			}

			if (current == ')' && innerParenthesisCount == 0) {
				result.push_back(']');
				startOfLog = false;
				parenthesesStarted = false;
				innerParenthesisCount = 0;
				continue;
			}

			// if it happens to be 'o', 'g' or somehting from argument
			// of the log, just push it back
			result.push_back(current);

		} else {
			// startOfLog = flase
			result.push_back(current);
		}

	}

	return result;
}

std::string rewriteGammaFunctionMathematicaFormat(GiNaC::ex _ex) {
	// the point is ginac uses tgamma() while mathematica requires EulerGamma[]
	// this function rewrites tgamma() -> EulerGamma[] and returns string that
	// can be then given to wolfram mathematica

	// Note works if there is no possibility of combination tg appearing in the term
	// at any other place than in tgamma -> our case

	std::string result;
	std::stringstream a;
	a << _ex;
	std::string input = a.str();

	char current;
	int startOfGamma = -10;
	bool tgEncountered = false;
	int innerParenthesisCount = -1;

	// create scanner on input expression (in string format)
	std::istringstream scanner(input);
	for (int i = 0; i < input.size(); i++) {

		current = scanner.get();

		if (current == 't') {
			startOfGamma = i;
			//result.push_back('L');
			continue;
		}

		if (current != 'g' && startOfGamma == i - 1) {
			result.push_back('t');
			result.push_back('^');
			continue;
		}

		if (current == 'g' && startOfGamma == i - 1) {
			tgEncountered = true;

			result = result + "EulerG";
			continue;
		}

		if (tgEncountered) {
			// if startOfLog = true and you find '(' or ')'

			if (current == '(' && innerParenthesisCount == -1) {
				result.push_back('[');
				innerParenthesisCount++;
				continue;
			}

			if (current == '(' && innerParenthesisCount != -1) {
				result.push_back(current);
				innerParenthesisCount++;
				continue;
			}

			if (current == ')' && innerParenthesisCount != 0) {
				result.push_back(current);
				innerParenthesisCount--;
				continue;
			}

			if (current == ')' && innerParenthesisCount == 0) {
				result.push_back(']');
				tgEncountered = false;
				innerParenthesisCount = -1;
				continue;
			}

			// if it happens to be amma or somehting from argument
			// of the tgamma, just push it back
			result.push_back(current);

		} else {
			// tgEncountered = false
			result.push_back(current);
		}

	}

	return result;
}

//===========================================================================

//===========================================================================
// code related to writing WLS output file - different NIntegrate methods
// coefficients from different FinalIntegrals are first summed and only then
// integrated - to be used with findDivergentPartPropToExtFreqWLS_together()
//
// NOTE: options of numerical method to be fixed inside functions

void writeWLSBasic(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::vector<std::string> strPoleCoefs = { };
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	for (int i = 0; i < _integVars.size(); i++) {
		std::cout << _integVars.at(i) << ", ";
	}

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::app);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "polePart = 0" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		strPoleCoefs.push_back(
				rewriteLogsIntoMathematicaFormat(_poleCoefs.at(i)));
	}

	for (int i = 0; i < strPoleCoefs.size(); i++) {

		currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

		writer << "polePart += 1*(eps^" << currentPowerOfEps << ") * "
				<< "NIntegrate[";
		writer << strPoleCoefs.at(i) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}
		writer << "]" << '\n';

	}

	writer << "Write[OpenWrite[\"" << _path << _diag.getName()
			<< "_basic.out\"], (overallFactor*polePart) ]" << '\n';
	writer << "Print[overallFactor*polePart]" << '\n';

	writer.close();

	return;
}

void writeWLSLocalAdaptive(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::vector<std::string> strPoleCoefs = { };
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);
	//writer.open(s.str(), std::fstream::app);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "polePart = 0" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		strPoleCoefs.push_back(
				rewriteLogsIntoMathematicaFormat(_poleCoefs.at(i)));
	}

	for (int i = 0; i < strPoleCoefs.size(); i++) {

		currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

		writer << "polePart += 1*(eps^" << currentPowerOfEps << ") * "
				<< "NIntegrate[";
		writer << strPoleCoefs.at(i) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}

		writer << ", Method->\"LocalAdaptive\" ]" << '\n';

	}

	writer << "Write[OpenWrite[\"" << _path << _diag.getName()
			<< "_LocAdapt.out\"], (overallFactor*polePart) ]" << '\n';
	writer << "Print[overallFactor*polePart]" << '\n';

	writer.close();

	return;
}

void writeWLSMC(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {
	//, int PRECISIONGOAL, int WORKINGPRECISION, int MAXRECURSION, int MAXPOINTS) {

	std::vector<std::string> strPoleCoefs = { };
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);
	//writer.open(s.str(), std::fstream::app);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "polePart" << -currentPowerOfEps << "= 0" << '\n';
		writer << "err" << -currentPowerOfEps << " = 0" << '\n';
	}

	for (int i = 0; i < _poleCoefs.size(); i++) {
		strPoleCoefs.push_back(
				rewriteLogsIntoMathematicaFormat(_poleCoefs.at(i)));
	}

	for (int i = 0; i < strPoleCoefs.size(); i++) {

		currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

		writer << "polePart" << -currentPowerOfEps << "= 1*(eps^"
				<< currentPowerOfEps << ") * " << "NIntegrate[";
		writer << strPoleCoefs.at(i) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}

		writer
				<< ", Method->\"MonteCarlo\", PrecisionGoal -> 4, MaxPoints->1000000, IntegrationMonitor :> ((errors = Through[#1@\"Error\"]) &)]"
				<< '\n';
		//	<< ", Method->\"AdaptiveMonteCarlo\", PrecisionGoal->4, WorkingPrecision ->7,"
		//			" MaxRecursion->50, MaxPoints -> 100000, IntegrationMonitor :> ((errors = Through[#1@\"Error\"]) &)]"
		//	<< '\n';

		writer << "err" << -currentPowerOfEps << "=  Total@errors" << "\n";

		writer << "Print[polePart" << -currentPowerOfEps << "]" << "\n";
		writer << "Print[err" << -currentPowerOfEps << "]" << "\n";
	}

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_MC.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "Write[op, polePart" << -currentPowerOfEps << "]" << '\n';
		//writer << "WriteString[op, \t ]"
		// << '\n';
	}

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "WriteString[op, \"err" << -currentPowerOfEps << "= \"]"
				<< '\n';
		writer << "Write[op, err" << -currentPowerOfEps << "]" << '\n';
	}

	writer << "Print[overallFactor*polePart]" << '\n';
	writer << "Close[op]" << '\n';

	writer.close();

	return;
}

void writeWLSCombinedMCAndEx(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	// writes WLS output file where highest pole is calculated exactly using
	// Integrate and other poles are calculated numerically using NIntegrate and
	// in particular using MC method

	std::vector<std::string> strPoleCoefs = { };
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);
	//writer.open(s.str(), std::fstream::app);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "start = Now" << "\n";

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "polePart" << -currentPowerOfEps << "= 0" << '\n';
		writer << "err" << -currentPowerOfEps << " = 0" << '\n';
	}

	for (int i = 0; i < _poleCoefs.size(); i++) {
		strPoleCoefs.push_back(
				rewriteLogsIntoMathematicaFormat(_poleCoefs.at(i)));
	}

	for (int i = 0; i < strPoleCoefs.size(); i++) {

		if (i != 0) {
			// if i != 0 -> subleading poles, Numerical calculation
			currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

			writer << "polePart" << -currentPowerOfEps << "= 1*(eps^"
					<< currentPowerOfEps << ") * " << "NIntegrate[";
			writer << strPoleCoefs.at(i) << ", ";

			for (int j = 0; j < _integVars.size(); j++) {
				writer << "{" << _integVars.at(j) << ",0,1}";
				if (j != _integVars.size() - 1) {
					writer << ", ";
				}
			}

			writer
					<< ", Method->\"MonteCarlo\", PrecisionGoal -> 4, MaxPoints->1000000, IntegrationMonitor :> ((errors = Through[#1@\"Error\"]) &)]"
					<< '\n';
			//	<< ", Method->\"AdaptiveMonteCarlo\", PrecisionGoal->4, WorkingPrecision ->7,"
			//			" MaxRecursion->50, MaxPoints -> 100000, IntegrationMonitor :> ((errors = Through[#1@\"Error\"]) &)]"
			//	<< '\n';

			writer << "err" << -currentPowerOfEps << "=  Total@errors" << "\n";
		} else {
			// if i ==0 -> highest pole
			// exact calculation using Integrate
			currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

			writer << "polePart" << -currentPowerOfEps << "= 1*(eps^"
					<< currentPowerOfEps << ") * " << "Integrate[";
			writer << strPoleCoefs.at(i) << ", ";

			for (int j = 0; j < _integVars.size(); j++) {
				writer << "{" << _integVars.at(j) << ",0,1}";
				if (j != _integVars.size() - 1) {
					writer << ", ";
				}
			}

			writer << "]" << '\n';

			writer << "err" << -currentPowerOfEps << "=  0" << "\n";
		}

		writer << "Print[polePart" << -currentPowerOfEps << "]" << "\n";
		writer << "Print[err" << -currentPowerOfEps << "]" << "\n";
	}

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_MC.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "Write[op, polePart" << -currentPowerOfEps << "]" << '\n';
		//writer << "WriteString[op, \t ]"
		// << '\n';
	}

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "WriteString[op, \"err" << -currentPowerOfEps << "= \"]"
				<< '\n';
		writer << "Write[op, err" << -currentPowerOfEps << "]" << '\n';
	}

	//writer << "Print[overallFactor*polePart]" << '\n';
	writer << "Close[op]" << '\n';
	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	return;
}

void writeWLSForExactHighestPoleCoef(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	std::stringstream s;
	int currentPowerOfEps;

	s << _path << "Diagram_" << _diag.getName() << "_eps2" << ".wls";

	std::fstream writer;
	writer.open(s.str(), std::fstream::out);
	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';
	currentPowerOfEps = -(_poleCoefs.size() - (0 + 1));
	writer << "coef" << -currentPowerOfEps << " = " << "Integrate[";
	writer << rewriteLogsIntoMathematicaFormat(_poleCoefs.at(0)) << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}
	}
	writer << "] \n";
	writer << "Print[coef" << -currentPowerOfEps << "]" << "\n";

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_eps2.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';
	//writer << "WriteString[op, \"exact coef. of pole eps^(" << currentPowerOfEps << "): \"]"
	//				<< '\n';
	writer << "Write[op, coef" << -currentPowerOfEps << "]" << '\n';
	writer << "Close[op]" << '\n';

	writer.close();

	return;
}

void writeWLSForExactHighestPoleCoef_NParts(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path, int _parts) {

	int N0 = _parts; // divide eps0 part into N0 parts
	std::vector<std::string> poleCoefParts = { };

	std::string poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(0));
	poleCoefParts = dividePoleCoefIntoNParts(poleCoef, N0);

	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	std::stringstream s;
	int currentPowerOfEps;

	for (int m = 0; m < N0; m++) {

		s.clear();
		s.str("");

		s << _path << "Diagram_" << _diag.getName() << "_eps2_part_" << m << ".wls";

		std::fstream writer;
		writer.open(s.str(), std::fstream::out);
		writer << "#!/usr/bin/env wolframscript" << '\n';

		writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';
		currentPowerOfEps = -(_poleCoefs.size() - (0 + 1));
		writer << "coef" << -currentPowerOfEps << " = " << "Integrate[";
		writer << poleCoefParts.at(m) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}
		writer << "] \n";
		writer << "Print[coef" << -currentPowerOfEps << "]" << "\n";

		writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
				<< "_eps2_part_"<< m <<".out\"]" << '\n';
		writer << "Write[op, overallFactor]" << '\n';
		//writer << "WriteString[op, \"exact coef. of pole eps^(" << currentPowerOfEps << "): \"]"
		//				<< '\n';
		writer << "Write[op, coef" << -currentPowerOfEps << "]" << '\n';
		writer << "Close[op]" << '\n';

		writer.close();
	}
	return;
}

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
		std::string _path) {

	std::vector<std::string> strPoleCoefs = { };
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);
	//writer.open(s.str(), std::fstream::app);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	// install Vegas
	writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]" << '\n';

	writer << "start = Now" << "\n";
	// parallelization
	//writer<<"MapSample = ParallelMap"<<'\n';

	// now actual things inside wls
	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "polePart" << -currentPowerOfEps << "= 0" << '\n';
	}

	for (int i = 0; i < _poleCoefs.size(); i++) {
		strPoleCoefs.push_back(
				rewriteLogsIntoMathematicaFormat(_poleCoefs.at(i)));
	}

	for (int i = 0; i < strPoleCoefs.size(); i++) {

		currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

		//writer << "polePart" << -currentPowerOfEps << "= 1*(eps^"
		//		<< currentPowerOfEps << ") * " << "Vegas[";

		writer << "polePart" << -currentPowerOfEps << "=" << "Vegas[";
		writer << strPoleCoefs.at(i) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}

		writer << ", PrecisionGoal -> 4, MaxPoints->100000]" << '\n';

		writer << "Print[polePart" << -currentPowerOfEps << "]" << "\n";
	}

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_Vegas.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));

		writer << "Write[op, 1*(eps^" << currentPowerOfEps << ") *polePart"
				<< -currentPowerOfEps << "]" << '\n';

		//writer << "Write[op, polePart" << -currentPowerOfEps << "]" << '\n';
	}

	writer << "Print[overallFactor*polePart]" << '\n';
	writer << "Close[op]" << '\n';

	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	return;
}

void writeWLSCombinedVegasAndEx(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	// writes WLS output file where highest pole is calculated exactly using
	// Integrate and other poles are calculated numerically using Vegas.
	// One common WLS file is created where all the poles are calculated
	// in series

	std::vector<std::string> strPoleCoefs = { };
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);
	//writer.open(s.str(), std::fstream::app);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]" << '\n';

	writer << "start = Now" << "\n";

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "polePart" << -currentPowerOfEps << "= 0" << '\n';
		//writer << "err" << -currentPowerOfEps << " = 0" << '\n';
	}

	for (int i = 0; i < _poleCoefs.size(); i++) {
		strPoleCoefs.push_back(
				rewriteLogsIntoMathematicaFormat(_poleCoefs.at(i)));
	}

	for (int i = 0; i < strPoleCoefs.size(); i++) {

		if (i != 0) {
			// if i != 0 -> subleading poles, Numerical calculation
			currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

			writer << "polePart" << -currentPowerOfEps << "= 1*(eps^"
					<< currentPowerOfEps << ") * " << "Vegas[";
			writer << strPoleCoefs.at(i) << ", ";

			for (int j = 0; j < _integVars.size(); j++) {
				writer << "{" << _integVars.at(j) << ",0,1}";
				if (j != _integVars.size() - 1) {
					writer << ", ";
				}
			}

			writer << ", PrecisionGoal -> 8, MaxPoints->10000]" << '\n';
			//	<< ", Method->\"AdaptiveMonteCarlo\", PrecisionGoal->4, WorkingPrecision ->7,"
			//			" MaxRecursion->50, MaxPoints -> 100000, IntegrationMonitor :> ((errors = Through[#1@\"Error\"]) &)]"
			//	<< '\n';

			//writer << "err" << -currentPowerOfEps << "=  Total@errors" << "\n";
		} else {
			// if i ==0 -> highest pole
			// exact calculation using Integrate
			currentPowerOfEps = -(strPoleCoefs.size() - (i + 1));

			writer << "polePart" << -currentPowerOfEps << "= 1*(eps^"
					<< currentPowerOfEps << ") * " << "Integrate[";
			writer << strPoleCoefs.at(i) << ", ";

			for (int j = 0; j < _integVars.size(); j++) {
				writer << "{" << _integVars.at(j) << ",0,1}";
				if (j != _integVars.size() - 1) {
					writer << ", ";
				}
			}

			writer << "]" << '\n';

			//writer << "err" << -currentPowerOfEps << "=  0"<<"\n";
		}

		writer << "Print[polePart" << -currentPowerOfEps << "]" << "\n";
		//writer << "Print[err" << -currentPowerOfEps << "]" << "\n";
	}

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_Vegas.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	for (int i = 0; i < _poleCoefs.size(); i++) {
		currentPowerOfEps = -(_poleCoefs.size() - (i + 1));
		writer << "Write[op, polePart" << -currentPowerOfEps << "]" << '\n';
		//writer << "WriteString[op, \t ]"
		// << '\n';
	}

	//writer << "Print[overallFactor*polePart]" << '\n';
	writer << "Close[op]" << '\n';
	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	return;
}

void writeWLSCombinedVegasAndEx_separately(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	// writes WLS output file where highest pole is calculated exactly using
	// Integrate and other poles are calculated numerically using Vegas.
	// Every pole is wrote into separate wls file.

	int mcSteps = 100000;
	std::string poleCoef = "";
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//-------------------------------------------------------------------
	// eps^{-2}
	//-------------------------------------------------------------------

	s << _path << "Diagram_" << _diag.getName() << "_eps2" << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]" << '\n';

	writer << "start = Now" << "\n";

	currentPowerOfEps = 2;
	writer << "polePart" << currentPowerOfEps << "= 0" << '\n';

	poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(0));
	writer << "polePart" << currentPowerOfEps << "= 1*(eps^"
			<< -currentPowerOfEps << ") * " << "Integrate[";
	writer << poleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}

	}

	writer << "]" << '\n';
	writer << "Print[polePart" << currentPowerOfEps << "]" << "\n";

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_Vegas_eps2.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	writer << "Write[op, polePart" << 2 << "]" << '\n';

	writer << "Close[op]" << '\n';
	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	s.clear();
	s.str("");
	s << _path << "Diagram_" << _diag.getName() << "_eps1" << ".wls";

	writer.open(s.str(), std::fstream::out);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]" << '\n';

	writer << "start = Now" << "\n";

	currentPowerOfEps = 1;
	poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(1));
	writer << "polePart" << currentPowerOfEps << "= 1*(eps^"
			<< -currentPowerOfEps << ") * " << "Vegas[";
	writer << poleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}
	}

	writer << ", PrecisionGoal -> 8, MaxPoints->" << mcSteps << "]" << '\n';

	writer << "Print[polePart" << currentPowerOfEps << "]" << "\n";

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_Vegas_eps1.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	writer << "Write[op, polePart" << 1 << "]" << '\n';

	writer << "Close[op]" << '\n';
	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// eps^{0}
	//-------------------------------------------------------------------

	s.clear();
	s.str("");

	s << _path << "Diagram_" << _diag.getName() << "_eps0" << ".wls";

	writer.open(s.str(), std::fstream::out);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]" << '\n';

	writer << "start = Now" << "\n";

	currentPowerOfEps = 0;
	poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(2));
	writer << "polePart" << currentPowerOfEps << "= 1*(eps^"
			<< currentPowerOfEps << ") * " << "Vegas[";
	writer << poleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}
	}

	writer << ", PrecisionGoal -> 8, MaxPoints->" << mcSteps << "]" << '\n';

	writer << "Print[polePart" << currentPowerOfEps << "]" << "\n";

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_Vegas_eps0.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	writer << "Write[op, polePart" << 0 << "]" << '\n';

	writer << "Close[op]" << '\n';
	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	return;
}

//===========================================================================

//===========================================================================
// code related to writing WLS output file - different NIntegrate methods
// coefficients from different FinalIntegrals are integrated and only after that
// summed - to be used with findDivergentPartPropToExtFreqWLS()
//
// NOTE: options of numerical method to be fixed inside functions

void writeBegginingWLS(Diagram _diag, GiNaC::ex _overallFactor,
		std::string _path) {

	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::app);
	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';
	writer << "polePart = 0" << '\n';
	writer.close();
	return;
}

void writePartialWLSBasic(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		GiNaC::ex _poleCoef, int _poleCoefsSize, int _indexInPoleCoefs,
		std::string _path) {

	std::string strPoleCoef;

	int currentPowerOfEps;

	for (int i = 0; i < _integVars.size(); i++) {
		//std::cout << _integVars.at(i) << ", ";
	}
	//std::cout<<std::endl;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::app);

	strPoleCoef = rewriteLogsIntoMathematicaFormat(_poleCoef);

	currentPowerOfEps = -(_poleCoefsSize - (_indexInPoleCoefs + 1));

	writer << "polePart += 1*(eps^" << currentPowerOfEps << ") * "
			<< "NIntegrate[";
	writer << strPoleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}
	}
	writer << "]" << '\n';

	writer.close();

	return;
}

void writePartialWLSLocalAdaptive(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars, GiNaC::ex _poleCoef,
		int _poleCoefsSize, int _indexInPoleCoefs, std::string _path) {

	std::string strPoleCoef;

	int currentPowerOfEps;

	for (int i = 0; i < _integVars.size(); i++) {
		//std::cout << _integVars.at(i) << ", ";
	}
	//std::cout<<std::endl;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::app);

	strPoleCoef = rewriteLogsIntoMathematicaFormat(_poleCoef);

	currentPowerOfEps = -(_poleCoefsSize - (_indexInPoleCoefs + 1));

	writer << "polePart += 1*(eps^" << currentPowerOfEps << ") * "
			<< "NIntegrate[";
	writer << strPoleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}
	}

	writer << ", Method->\"LocalAdaptive\" ]" << '\n';

	writer.close();

	return;
}

void writePartialWLSMC(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		GiNaC::ex _poleCoef, int _poleCoefsSize, int _indexInPoleCoefs,
		std::string _path) {

	std::string strPoleCoef;

	int currentPowerOfEps;

	for (int i = 0; i < _integVars.size(); i++) {
		//std::cout << _integVars.at(i) << ", ";
	}
	//std::cout<<std::endl;

	// write záhlavok
	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::app);

	strPoleCoef = rewriteLogsIntoMathematicaFormat(_poleCoef);

	currentPowerOfEps = -(_poleCoefsSize - (_indexInPoleCoefs + 1));

	writer << "polePart += 1*(eps^" << currentPowerOfEps << ") * "
			<< "NIntegrate[";
	writer << strPoleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}
	}

	writer << ", Method->\"MonteCarlo\"]" << '\n';

	writer.close();

	return;
}

void writeEndingWLS(Diagram _diag, std::string _path) {

	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << _path << "Diagram_" << _diag.getName() << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::app);

	writer << "Write[OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_basic.out\"], (overallFactor*polePart) ]" << '\n';
	writer << "Print[overallFactor*polePart]" << '\n';

	writer.close();
	return;
}

//===========================================================================

//===========================================================================
// Dividing the calculation into more parts Mathematica

void writeWLSCombinedExVegas_Nparts(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	int N1 = 3; //divide eps1 into N1 parts
	int N0 = 10; // divide eps0 part into N0 parts
	int mcSteps = 10000;

	std::string poleCoef = "";
	std::string overallFactor = rewriteGammaFunctionMathematicaFormat(
			_overallFactor);
	int currentPowerOfEps;

	// write záhlavok
	std::stringstream s;

	//-------------------------------------------------------------------
	// eps^{-2}
	//-------------------------------------------------------------------

	s << _path << "Diagram_" << _diag.getName() << "_eps2" << ".wls";

	std::fstream writer;

	writer.open(s.str(), std::fstream::out);

	writer << "#!/usr/bin/env wolframscript" << '\n';

	writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

	writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]" << '\n';

	writer << "start = Now" << "\n";

	currentPowerOfEps = 2;
	writer << "polePart" << currentPowerOfEps << "= 0" << '\n';

	poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(0));
	writer << "polePart" << currentPowerOfEps << "= 1*(eps^"
			<< -currentPowerOfEps << ") * " << "Integrate[";
	writer << poleCoef << ", ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "{" << _integVars.at(j) << ",0,1}";
		if (j != _integVars.size() - 1) {
			writer << ", ";
		}

	}

	writer << "]" << '\n';
	writer << "Print[polePart" << currentPowerOfEps << "]" << "\n";

	writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
			<< "_Vegas_eps2.out\"]" << '\n';
	writer << "Write[op, overallFactor]" << '\n';

	writer << "Write[op, polePart" << 2 << "]" << '\n';

	writer << "Close[op]" << '\n';
	writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
	writer.close();

	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	std::vector<std::string> poleCoefParts;

	currentPowerOfEps = 1;
	poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(1));
	poleCoefParts = dividePoleCoefIntoNParts(poleCoef, N1);

	for (int i = 0; i < N1; i++) {
		s.clear();
		s.str("");
		s << _path << "Diagram_" << _diag.getName() << "_eps1_part" << i
				<< ".wls";

		writer.open(s.str(), std::fstream::out);

		writer << "#!/usr/bin/env wolframscript" << '\n';

		writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

		writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]"
				<< '\n';

		writer << "start = Now" << "\n";

		writer << "polePart" << currentPowerOfEps << "= 1*(eps^"
				<< -currentPowerOfEps << ") * " << "Vegas[";
		writer << poleCoefParts.at(i) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}

		writer << ", PrecisionGoal -> 8, MaxPoints->" << mcSteps << "]" << '\n';

		writer << "Print[polePart" << currentPowerOfEps << "]" << "\n";

		writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
				<< "_Vegas_eps1_" << i << ".out\"]" << '\n';
		writer << "Write[op, overallFactor]" << '\n';

		writer << "Write[op, polePart" << 1 << "]" << '\n';

		writer << "Close[op]" << '\n';
		writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
		writer.close();
	}
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// eps^{0}
	//-------------------------------------------------------------------

	poleCoefParts = { }; // reset the vector
	currentPowerOfEps = 0;
	poleCoef = rewriteLogsIntoMathematicaFormat(_poleCoefs.at(2));
	poleCoefParts = dividePoleCoefIntoNParts(poleCoef, N0);

	for (int i = 0; i < N0; i++) {

		s.clear();
		s.str("");

		s << _path << "Diagram_" << _diag.getName() << "_eps0_part" << i
				<< ".wls";

		writer.open(s.str(), std::fstream::out);

		writer << "#!/usr/bin/env wolframscript" << '\n';

		writer << "overallFactor = HoldForm[" << overallFactor << "]" << '\n';

		writer << "Install[ \"/home/matej/Downloads/Cuba-4.2.2/Vegas\"]"
				<< '\n';

		writer << "start = Now" << "\n";

		writer << "polePart" << currentPowerOfEps << "= 1*(eps^"
				<< currentPowerOfEps << ") * " << "Vegas[";
		writer << poleCoefParts.at(i) << ", ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "{" << _integVars.at(j) << ",0,1}";
			if (j != _integVars.size() - 1) {
				writer << ", ";
			}
		}

		writer << ", PrecisionGoal -> 8, MaxPoints->" << mcSteps << "]" << '\n';

		writer << "Print[polePart" << currentPowerOfEps << "]" << "\n";

		writer << "op = OpenWrite[\"" << _path << "Diagram_" << _diag.getName()
				<< "_Vegas_eps0_part" << i << ".out\"]" << '\n';
		writer << "Write[op, overallFactor]" << '\n';

		writer << "Write[op, polePart" << 0 << "]" << '\n';

		writer << "Close[op]" << '\n';
		writer << "Print[DateDifference[start, Now, \"Second\"] ]" << "\n";
		writer.close();
	}
}

int countParts(std::string _poleCoef) {
	// Function counts how many terms make up _pole coef
	// (e.g. term1 + term2 + term3 + ... + term N) -> function would return N

	int numOfParts = 0;

	char current;
	bool startOfLog = false;
	bool parenthesesStarted = false;
	int innerParenthesisCount = 0;
	int aux = 0;

	std::istringstream scanner(_poleCoef);

	for (int i = 0; i < _poleCoef.size(); i++) {
		current = scanner.get();

		// ------
		// find out if you are in a log
		if (current == 'L') {
			aux = i;
			continue;
		}

		if (current == 'o' && aux == i - 1) {
			continue;
		}

		if (current == 'g' && aux == i - 2) {
			startOfLog = true;
			continue;
		}

		if (startOfLog && current == ']') {
			startOfLog = false;
			continue;
		}
		// ------

		// ------
		// find if you are inside a parenthesis
		if (!startOfLog && current == '(') {
			innerParenthesisCount++;
			continue;
		}

		if (!startOfLog && current == ')' && innerParenthesisCount != 0) {
			innerParenthesisCount--;
			continue;
		}

		// ------

		// now count pluses
		if (current == '+' && !startOfLog && innerParenthesisCount == 0) {
			numOfParts++;
			continue;
		}

	}
	return numOfParts;
}

void divideTermIntoNParts(std::string _poleCoef, int _N,
		std::vector<std::string> &_vector) {

	int currentPlus = 0;
	std::string currentTerm = "";
	std::string auxString = "";
	int numOfSubterms = countParts(_poleCoef);
	int numOfPushedBackTerms = 0;

	char current;
	bool startOfLog = false;
	bool parenthesesStarted = false;
	int innerParenthesisCount = 0;
	int aux = 0;

	std::istringstream scanner(_poleCoef);

	for (int i = 0; i < _poleCoef.size(); i++) {
		current = scanner.get();
		currentTerm.push_back(current);
		// ------
		// find out if you are in a log
		if (current == 'L') {
			aux = i;
			continue;
		}

		if (current == 'o' && aux == i - 1) {
			continue;
		}

		if (current == 'g' && aux == i - 2) {
			startOfLog = true;
			continue;
		}

		if (startOfLog && current == ']') {
			startOfLog = false;
			continue;
		}
		// ------

		// ------
		// find if you are inside a parenthesis
		if (!startOfLog && current == '(') {
			innerParenthesisCount++;
			continue;
		}

		if (!startOfLog && current == ')' && innerParenthesisCount != 0) {
			innerParenthesisCount--;
			continue;
		}

		// ------

		// now count pluses
		if (current == '+' && !startOfLog && innerParenthesisCount == 0) {
			currentPlus++;
			if (currentPlus == numOfSubterms / _N
					&& numOfPushedBackTerms < (_N - 1)) {
				for (int m = 0; m < currentTerm.size() - 1; m++) {
					auxString.push_back(currentTerm.at(m));
				}
				_vector.push_back(auxString);
				numOfPushedBackTerms++;
				currentTerm = "";
				auxString = "";
				currentPlus = 0;
			}
			continue;
		}

	}

	for (int m = 0; m < currentTerm.size(); m++) {
		auxString.push_back(currentTerm.at(m));
	}
	_vector.push_back(auxString);

	if (_vector.size() != _N) {
		std::cout << _vector.size() << std::endl;
		std::cout
				<< "Something went wrong in divide into N parts, OutputWriter.cpp"
				<< std::endl;
	}

	return;
}

std::vector<std::string> dividePoleCoefIntoNParts(std::string _poleCoef,
		int _N) {

	std::vector<std::string> result;

	//int numOfParts = countParts(_poleCoef);

	divideTermIntoNParts(_poleCoef, _N, result);

	return result;
}

//===========================================================================

//===========================================================================
//===========================================================================
// Dividing the calculation into more parts C vegas

void writeCVegas_Eps0_Nparts(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path, int _parts) {

	int N0 = _parts; // divide eps0 part into N0 parts

	std::stringstream aux;
	std::string poleCoef = "";

	std::stringstream a;
	a << "\"" << _overallFactor << "\"";
	std::string overallFactor = a.str();

	int currentPowerOfEps;

	std::fstream writer;

	// write záhlavok
	std::stringstream s;

	//-------------------------------------------------------------------
	// eps^{0}
	//-------------------------------------------------------------------

	std::vector<std::string> poleCoefParts = { }; // reset the vector
	currentPowerOfEps = 0;
	_poleCoefs.at(2).print(GiNaC::print_csrc(aux));
	poleCoef = aux.str();
	poleCoefParts = dividePoleCoefIntoNParts_CVegas(poleCoef, N0);

	std::string testPoleCoefs = "";
	for (int m = 0; m < N0; m++) {
		testPoleCoefs += poleCoefParts.at(m);
		testPoleCoefs += "+";
	}
	testPoleCoefs.pop_back();
	if (testPoleCoefs.compare(poleCoef) != 0) {
		std::cout << "Something went wrong with division into N parts."
				<< std::endl;
	}

	for (int m = 0; m < N0; m++) {

		s.clear();
		s.str("");

		s << _path << "Diagram_" << _diag.getName() << "_eps0_part_" << m
				<< ".cpp";

		writer.open(s.str(), std::fstream::out);

		writer << "#include <stdio.h> \n"
				"#include <stdlib.h> \n"
				"#include <math.h> \n"
				"#include <iostream> \n"
				"#include <fstream> \n"
				"#include <sstream> \n"
				"#include\"/home/matej/Downloads/Cuba-4.2.2/cuba.h\" \n \n";

		writer << "#define NDIM " << _integVars.size() << "\n"
				"#define NCOMP 1 \n"
				"#define USERDATA NULL \n"
				"#define NVEC 1 \n"
				"#define EPSREL 1e-7 \n"
				"#define EPSABS 1e-12 \n"
				"#define VERBOSE 2 \n"
				"#define LAST 4 \n"
				"#define SEED 0 \n"
				"#define MINEVAL 0 \n"
				"#define MAXEVAL 10000000 \n"
				"\n"
				"#define NSTART 1000 \n"
				"#define NINCREASE 500 \n"
				"#define NBATCH 1000 \n"
				"#define GRIDNO 0 \n"
				"#define STATEFILE NULL \n"
				"#define SPIN NULL \n"
				"\n"
				"#define KEY 0 \n";

		writer
				<< "static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) { \n \n ";

		for (int j = 0; j < _integVars.size(); j++) {
			writer << "#define " << _integVars.at(j) << " xx[" << j << "] \n";
		}

		writer << "#define f ff[0] \n \n";

		writer << "f = " << poleCoefParts.at(m) << "; \n \n";
		//writer << "f = 1./4*t0*pow((1+2*t5*t0*pow(t1,2)+6*t4*pow(t2,2)*t0+8*t2*t5*t0*t1+t0*t1+2*t4+8*t4*t2*t0*t1+2*t2+2*t2*t0+6*pow(t2,2)*t5*t0+t0*pow(t1,2)+2*t5*t0*t1+2*t5+4*t4*t2+2*t4*t0*pow(t1,2)+4*t2*t0*t1+2*t5*t1+3*pow(t2,2)*t0+t1+4*t4*t2*t0+2*t4*t1+4*t2*t5*t0+2*t4*t0*t1+4*t2*t5),(-3))" << "; \n \n";

		writer << " return 0; \n"
				"}\n";

		writer
				<< "int main() { \n"
						"clock_t begin = clock(); \n"
						"int comp, nregions, neval, fail; cubareal integral[NCOMP], error[NCOMP], prob[NCOMP]; \n"
						"\n"
						"Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, \n"
						"   EPSREL, EPSABS, VERBOSE, SEED, \n"
						"	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, \n"
						"    GRIDNO, STATEFILE, SPIN, \n"
						"   &neval, &fail, integral, error, prob); "
						"\n"
						"printf(\"VEGAS RESULT:\\tneval %d\\tfail %d\\n\", neval, fail); \n"
						"for( comp = 0; comp < NCOMP; ++comp ) \n"
						"			printf(\"VEGAS RESULT:\\t%.8f +- %.8f\\tp = %.3f\\n\", "
						"(double)integral[comp], (double)error[comp], (double)prob[comp]); \n \n"
						"\n"
						"clock_t end = clock(); \n"
						"std::cout << \"Elapsed time: \" << double(end - begin) / CLOCKS_PER_SEC"
						"<<\"s.\"<< std::endl; \n"
						"std::stringstream s; \n"
						"std::fstream writer;\n";
		writer << "s<<\"" << _path << "Diagram_" << _diag.getName()
				<< "_eps0_part_" << m << ".out\"; \n";
		writer << "writer.open(s.str(), std::fstream::out); \n"
				"writer<<std::fixed; \n"
				"writer.precision(10); \n";
		writer << "writer<<" << overallFactor << "<< \" \\n \"; \n";
		writer
				<< "for (comp = 0; comp < NCOMP; ++comp){ \n"
						"	writer<<(double) integral[comp]<<\"\t\"<<(double) error[comp]<<\"\t\"<<(double) prob[comp]; \n"
						"}; \n";
		writer << "return 0; \n} \n";
		writer.close();

	}

	return;
}

int countParts_CVegas(std::string _poleCoef) {
	// Function counts how many terms make up _pole coef
	// (e.g. term1 + term2 + term3 + ... + term N) -> function would return N

	int numOfParts = 0;

	char current;
	bool startOfLog = false;
	bool parenthesesStarted = false;
	int innerParenthesisCount = 0;
	int aux = 0;

	std::istringstream scanner(_poleCoef);

	for (int i = 0; i < _poleCoef.size(); i++) {
		current = scanner.get();

		// ------
		// find if you are inside a parenthesis
		if (current == '(') {
			innerParenthesisCount++;
			continue;
		}

		if (current == ')' && innerParenthesisCount != 0) {
			innerParenthesisCount--;
			continue;
		}

		// ------

		// now count pluses
		//if ((current == '+' || current == '-') && innerParenthesisCount == 0) {
		if (current == '+' && innerParenthesisCount == 0) {

			numOfParts++;
			continue;
		}

	}
	return numOfParts;
}

void divideTermIntoNParts_CVegas(std::string _poleCoef, int _N,
		std::vector<std::string> &_vector) {

	int currentPlus = 0;
	std::string currentTerm = "";
	std::string auxString = "";
	int numOfSubterms = countParts_CVegas(_poleCoef);
	int numOfPushedBackTerms = 0;

	char current;
	bool startOfLog = false;
	bool parenthesesStarted = false;
	int innerParenthesisCount = 0;
	int aux = 0;

	//std::istringstream scanner(_poleCoef);

	for (int i = 0; i < _poleCoef.size(); i++) {
		current = _poleCoef.at(i);
		currentTerm.push_back(current);
		// ------

		// ------
		// find if you are inside a parenthesis
		if (current == '(') {
			innerParenthesisCount++;
			continue;
		}

		if (current == ')' && innerParenthesisCount != 0) {
			innerParenthesisCount--;
			continue;
		}

		// ------

		// now count pluses
		//if ((current == '+' || current == '-') && innerParenthesisCount == 0) {
		if (current == '+' && innerParenthesisCount == 0) {
			currentPlus++;
			if (currentPlus == numOfSubterms / _N
					&& numOfPushedBackTerms < (_N - 1)) {
				for (int m = 0; m < currentTerm.size() - 1; m++) {
					auxString.push_back(currentTerm.at(m));
				}
				_vector.push_back(auxString);
				numOfPushedBackTerms++;
				currentTerm = "";
				auxString = "";
				currentPlus = 0;
			}
			continue;
		}

	}

	for (int m = 0; m < currentTerm.size(); m++) {
		auxString.push_back(currentTerm.at(m));
	}
	_vector.push_back(auxString);

	if (_vector.size() != _N) {
		std::cout << _vector.size() << std::endl;
		std::cout
				<< "Something went wrong in divide into N parts, OutputWriter.cpp"
				<< std::endl;
	}

	return;
}

std::vector<std::string> dividePoleCoefIntoNParts_CVegas(std::string _poleCoef,
		int _N) {

	std::vector<std::string> result;

	//int numOfParts = countParts(_poleCoef);

	divideTermIntoNParts_CVegas(_poleCoef, _N, result);

	return result;
}

//===========================================================================

//===========================================================================

//===========================================================================
// code related to writing output file for Vegas implementation in C

void writeCombinedCVegasAndEx_3loop(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	//writeWLSForExactHighestPoleCoef_NParts(_diag, _integVars, _poleCoefs, _overallFactor, _path, 4);
	writeWLSForExactHighestPoleCoef(_diag, _integVars, _poleCoefs,
			_overallFactor, _path);
	writeCVegas_eps1_3loop(_diag, _integVars, _poleCoefs, _overallFactor, _path);
	//writeCVegas_eps0_3loop(_diag, _integVars, _poleCoefs, _overallFactor, _path);
	writeCVegas_Eps0_Nparts(_diag, _integVars, _poleCoefs, _overallFactor,
			_path, 5);
}

void writeCombinedCVegasAndEx_2loop(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	writeWLSForExactHighestPoleCoef(_diag, _integVars, _poleCoefs,
			_overallFactor, _path);
	writeCVegas_eps0_2loop(_diag, _integVars, _poleCoefs, _overallFactor, _path);
}
void writeCombinedCVegasAndEx_1loop(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	writeWLSForExactHighestPoleCoef(_diag, _integVars, _poleCoefs,
			_overallFactor, _path);
	//writeCVegas_eps0_2loop(_diag, _integVars, _poleCoefs, _overallFactor, _path);
}

void writeCVegas_eps1_3loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::stringstream aux;
	std::string poleCoef = "";

	std::stringstream a;
	a << "\"" << _overallFactor << "\"";
	std::string overallFactor = a.str();

	int currentPowerOfEps;

	std::fstream writer;

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	currentPowerOfEps = 1;
	_poleCoefs.at(1).print(GiNaC::print_csrc(aux));
	poleCoef = aux.str();
	// write záhlavok
	// TODO

	std::stringstream s;
	s << _path << "Diagram_" << _diag.getName() << "_eps1" << ".cpp";

	writer.open(s.str(), std::fstream::out);

	writer << "#include <stdio.h> \n"
			"#include <stdlib.h> \n"
			"#include <math.h> \n"
			"#include <iostream> \n"
			"#include <fstream> \n"
			"#include <sstream> \n"
			"#include\"/home/matej/Downloads/Cuba-4.2.2/cuba.h\" \n \n";

	writer << "#define NDIM " << _integVars.size() << "\n"
			"#define NCOMP 1 \n"
			"#define USERDATA NULL \n"
			"#define NVEC 1 \n"
			"#define EPSREL 1e-7 \n"
			"#define EPSABS 1e-12 \n"
			"#define VERBOSE 2 \n"
			"#define LAST 4 \n"
			"#define SEED 0 \n"
			"#define MINEVAL 0 \n"
			"#define MAXEVAL 10000000 \n"
			"\n"
			"#define NSTART 1000 \n"
			"#define NINCREASE 500 \n"
			"#define NBATCH 1000 \n"
			"#define GRIDNO 0 \n"
			"#define STATEFILE NULL \n"
			"#define SPIN NULL \n"
			"\n"
			"#define KEY 0 \n";

	writer
			<< "static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) { \n \n ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "#define " << _integVars.at(j) << " xx[" << j << "] \n";
	}

	writer << "#define f ff[0] \n \n";

	writer << "f = " << poleCoef << "; \n \n";
	//writer << "f = 1./4*t0*pow((1+2*t5*t0*pow(t1,2)+6*t4*pow(t2,2)*t0+8*t2*t5*t0*t1+t0*t1+2*t4+8*t4*t2*t0*t1+2*t2+2*t2*t0+6*pow(t2,2)*t5*t0+t0*pow(t1,2)+2*t5*t0*t1+2*t5+4*t4*t2+2*t4*t0*pow(t1,2)+4*t2*t0*t1+2*t5*t1+3*pow(t2,2)*t0+t1+4*t4*t2*t0+2*t4*t1+4*t2*t5*t0+2*t4*t0*t1+4*t2*t5),(-3))" << "; \n \n";

	writer << " return 0; \n"
			"} \n";

	writer
			<< "int main() { \n"
					"clock_t begin = clock(); \n"
					"int comp, nregions, neval, fail; cubareal integral[NCOMP], error[NCOMP], prob[NCOMP]; \n"
					"\n"
					"Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, \n"
					"   EPSREL, EPSABS, VERBOSE, SEED, \n"
					"	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, \n"
					"    GRIDNO, STATEFILE, SPIN, \n"
					"   &neval, &fail, integral, error, prob); "
					"\n"
					"printf(\"VEGAS RESULT:\\tneval %d\\tfail %d\\n\", neval, fail); \n"
					"for( comp = 0; comp < NCOMP; ++comp ) \n"
					"			printf(\"VEGAS RESULT:\\t%.8f +- %.8f\\tp = %.3f\\n\", "
					"(double)integral[comp], (double)error[comp], (double)prob[comp]); \n \n"
					"\n"
					"clock_t end = clock(); \n"
					"std::cout << \"Elapsed time: \" << double(end - begin) / CLOCKS_PER_SEC"
					"<<\"s.\"<< std::endl; \n"
					"std::stringstream s; \n"
					"std::fstream writer;\n";
	writer << "s<<\"" << _path << "Diagram_" << _diag.getName() << "_eps1"
			<< ".out\"; \n";
	writer << "writer.open(s.str(), std::fstream::out); \n"
			"writer<<std::fixed; \n"
			"writer.precision(10); \n";
	writer << "writer<<" << overallFactor << "<< \" \\n \"; \n";
	writer
			<< "for (comp = 0; comp < NCOMP; ++comp){ \n"
					"	writer<<(double) integral[comp]<<\"\t\"<<(double) error[comp]<<\"\t\"<<(double) prob[comp]; \n"
					"}; \n"
					"return 0; \n} "
					"\n";
	//-------------------------------------------------------------------

	writer.close();
	return;
}

void writeCVegas_eps0_3loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::stringstream aux;
	std::string poleCoef = "";

	std::stringstream a;
	a << "\"" << _overallFactor << "\"";
	std::string overallFactor = a.str();

	int currentPowerOfEps;

	std::fstream writer;

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	currentPowerOfEps = 0;
	_poleCoefs.at(2).print(GiNaC::print_csrc(aux));
	poleCoef = aux.str();
	// write záhlavok
	// TODO

	std::stringstream s;
	s << _path << "Diagram_" << _diag.getName() << "_eps0" << ".cpp";

	writer.open(s.str(), std::fstream::out);

	writer << "#include <stdio.h> \n"
			"#include <stdlib.h> \n"
			"#include <math.h> \n"
			"#include <iostream> \n"
			"#include <fstream> \n"
			"#include <sstream> \n"
			"#include\"/home/matej/Downloads/Cuba-4.2.2/cuba.h\" \n \n";

	writer << "#define NDIM " << _integVars.size() << "\n"
			"#define NCOMP 1 \n"
			"#define USERDATA NULL \n"
			"#define NVEC 1 \n"
			"#define EPSREL 1e-7 \n"
			"#define EPSABS 1e-12 \n"
			"#define VERBOSE 2 \n"
			"#define LAST 4 \n"
			"#define SEED 0 \n"
			"#define MINEVAL 0 \n"
			"#define MAXEVAL 10000000 \n"
			"\n"
			"#define NSTART 1000 \n"
			"#define NINCREASE 500 \n"
			"#define NBATCH 1000 \n"
			"#define GRIDNO 0 \n"
			"#define STATEFILE NULL \n"
			"#define SPIN NULL \n"
			"\n"
			"#define KEY 0 \n";

	writer
			<< "static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) { \n \n ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "#define " << _integVars.at(j) << " xx[" << j << "] \n";
	}

	writer << "#define f ff[0] \n \n";

	writer << "f = " << poleCoef << "; \n \n";
	//writer << "f = 1./4*t0*pow((1+2*t5*t0*pow(t1,2)+6*t4*pow(t2,2)*t0+8*t2*t5*t0*t1+t0*t1+2*t4+8*t4*t2*t0*t1+2*t2+2*t2*t0+6*pow(t2,2)*t5*t0+t0*pow(t1,2)+2*t5*t0*t1+2*t5+4*t4*t2+2*t4*t0*pow(t1,2)+4*t2*t0*t1+2*t5*t1+3*pow(t2,2)*t0+t1+4*t4*t2*t0+2*t4*t1+4*t2*t5*t0+2*t4*t0*t1+4*t2*t5),(-3))" << "; \n \n";

	writer << " return 0; \n"
			"} \n";

	writer
			<< "int main() { \n"
					"clock_t begin = clock(); \n"
					"int comp, nregions, neval, fail; cubareal integral[NCOMP], error[NCOMP], prob[NCOMP]; \n"
					"\n"
					"Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, \n"
					"   EPSREL, EPSABS, VERBOSE, SEED, \n"
					"	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, \n"
					"    GRIDNO, STATEFILE, SPIN, \n"
					"   &neval, &fail, integral, error, prob); "
					"\n"
					"printf(\"VEGAS RESULT:\\tneval %d\\tfail %d\\n\", neval, fail); \n"
					"for( comp = 0; comp < NCOMP; ++comp ) \n"
					"			printf(\"VEGAS RESULT:\\t%.8f +- %.8f\\tp = %.3f\\n\", "
					"(double)integral[comp], (double)error[comp], (double)prob[comp]); \n \n"
					"\n"
					"clock_t end = clock(); \n"
					"std::cout << \"Elapsed time: \" << double(end - begin) / CLOCKS_PER_SEC"
					"<<\"s.\"<< std::endl; \n"
					"std::stringstream s; \n"
					"std::fstream writer;\n";
	writer << "s<<\"" << _path << "Diagram_" << _diag.getName() << "_eps0"
			<< ".out\"; \n";
	writer << "writer.open(s.str(), std::fstream::out); \n"
			"writer<<std::fixed; \n"
			"writer.precision(10); \n";
	writer << "writer<<" << overallFactor << "<< \" \\n \"; \n";
	writer
			<< "for (comp = 0; comp < NCOMP; ++comp){ \n"
					"	writer<<(double) integral[comp]<<\"\t\"<<(double) error[comp]<<\"\t\"<<(double) prob[comp]; \n"
					"}; \n"
					"return 0; \n} "
					"\n";
	//-------------------------------------------------------------------

	writer.close();
	return;
}

void writeCVegas_eps0_2loop(Diagram _diag, std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _poleCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::stringstream aux;
	std::string poleCoef = "";

	std::stringstream a;
	a << "\"" << _overallFactor << "\"";
	std::string overallFactor = a.str();

	int currentPowerOfEps;

	std::fstream writer;

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	currentPowerOfEps = 1;
	_poleCoefs.at(1).print(GiNaC::print_csrc(aux));
	poleCoef = aux.str();
	// write záhlavok
	// TODO

	std::stringstream s;
	s << _path << "Diagram_" << _diag.getName() << "_eps1" << ".cpp";

	writer.open(s.str(), std::fstream::out);

	writer << "#include <stdio.h> \n"
			"#include <stdlib.h> \n"
			"#include <math.h> \n"
			"#include <iostream> \n"
			"#include <fstream> \n"
			"#include <sstream> \n"
			"#include\"/home/matej/Downloads/Cuba-4.2.2/cuba.h\" \n \n";

	writer << "#define NDIM " << _integVars.size() << "\n"
			"#define NCOMP 1 \n"
			"#define USERDATA NULL \n"
			"#define NVEC 1 \n"
			"#define EPSREL 1e-7 \n"
			"#define EPSABS 1e-12 \n"
			"#define VERBOSE 2 \n"
			"#define LAST 4 \n"
			"#define SEED 0 \n"
			"#define MINEVAL 0 \n"
			"#define MAXEVAL 10000000 \n"
			"\n"
			"#define NSTART 1000 \n"
			"#define NINCREASE 500 \n"
			"#define NBATCH 1000 \n"
			"#define GRIDNO 0 \n"
			"#define STATEFILE NULL \n"
			"#define SPIN NULL \n"
			"\n"
			"#define KEY 0 \n";

	writer
			<< "static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) { \n \n ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "#define " << _integVars.at(j) << " xx[" << j << "] \n";
	}

	writer << "#define f ff[0] \n \n";

	writer << "f = " << poleCoef << "; \n \n";
	//writer << "f = 1./4*t0*pow((1+2*t5*t0*pow(t1,2)+6*t4*pow(t2,2)*t0+8*t2*t5*t0*t1+t0*t1+2*t4+8*t4*t2*t0*t1+2*t2+2*t2*t0+6*pow(t2,2)*t5*t0+t0*pow(t1,2)+2*t5*t0*t1+2*t5+4*t4*t2+2*t4*t0*pow(t1,2)+4*t2*t0*t1+2*t5*t1+3*pow(t2,2)*t0+t1+4*t4*t2*t0+2*t4*t1+4*t2*t5*t0+2*t4*t0*t1+4*t2*t5),(-3))" << "; \n \n";

	writer << " return 0; \n"
			"} \n";

	writer
			<< "int main() { \n"
					"clock_t begin = clock(); \n"
					"int comp, nregions, neval, fail; cubareal integral[NCOMP], error[NCOMP], prob[NCOMP]; \n"
					"\n"
					"Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, \n"
					"   EPSREL, EPSABS, VERBOSE, SEED, \n"
					"	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, \n"
					"    GRIDNO, STATEFILE, SPIN, \n"
					"   &neval, &fail, integral, error, prob); "
					"\n"
					"printf(\"VEGAS RESULT:\\tneval %d\\tfail %d\\n\", neval, fail); \n"
					"for( comp = 0; comp < NCOMP; ++comp ) \n"
					"			printf(\"VEGAS RESULT:\\t%.8f +- %.8f\\tp = %.3f\\n\", "
					"(double)integral[comp], (double)error[comp], (double)prob[comp]); \n \n"
					"\n"
					"clock_t end = clock(); \n"
					"std::cout << \"Elapsed time: \" << double(end - begin) / CLOCKS_PER_SEC"
					"<<\"s.\"<< std::endl; \n"
					"std::stringstream s; \n"
					"std::fstream writer;\n";
	writer << "s<<\"" << _path << "Diagram_" << _diag.getName() << "_eps1"
			<< ".out\"; \n";
	writer << "writer.open(s.str(), std::fstream::out); \n"
			"writer<<std::fixed; \n"
			"writer.precision(10); \n";
	writer << "writer<<" << overallFactor << "<< \" \\n \"; \n";
	writer
			<< "for (comp = 0; comp < NCOMP; ++comp){ \n"
					"	writer<<(double) integral[comp]<<\"\t\"<<(double) error[comp]<<\"\t\"<<(double) prob[comp]; \n"
					"}; \n"
					"return 0; \n} "
					"\n";
	//-------------------------------------------------------------------

	writer.close();
	return;
}

//===========================================================================
// code related to writing output file for Vegas implementation in C
void writeCVegas_Finite_2loop(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _finitePartCoefs, GiNaC::ex _overallFactor,
		std::string _path, int _order) {

	if (_order > 1) {
		std::cout
				<< "You need to write implementation for it. Up to this moment it only works"
						"for Finite parts overall eps^0 up to overall eps^1. \n -outputWriter.cpp"
				<< std::endl;
		return;
	}

	if (_order == 0) {
		writeFiniteCVegas_eps0(_diag, _integVars, _finitePartCoefs,
				_overallFactor, _path);
	}

	if (_order == 1) {
		writeFiniteCVegas_eps0(_diag, _integVars, _finitePartCoefs,
				_overallFactor, _path);
		writeFiniteCVegas_eps1(_diag, _integVars, _finitePartCoefs,
				_overallFactor, _path);
	}

	return;
}

void writeFiniteCVegas_eps0(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _finitePartCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::stringstream aux;
	std::string coef = "";

	std::stringstream a;
	// the *eps^1 to make sure the result including overall factor is really order eps^0
	a << "\"" << _overallFactor << "*eps^1 \"";
	std::string overallFactor = a.str();

	//int currentPowerOfEps;

	std::fstream writer;

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	//currentPowerOfEps = 1;
	_finitePartCoefs.at(0).print(GiNaC::print_csrc(aux));
	coef = aux.str();
	// write záhlavok
	// TODO

	std::stringstream s;
	s << _path << "Diagram_" << _diag.getName() << "_finite_eps0" << ".cpp";

	writer.open(s.str(), std::fstream::out);

	writer << "#include <stdio.h> \n"
			"#include <stdlib.h> \n"
			"#include <math.h> \n"
			"#include <iostream> \n"
			"#include <fstream> \n"
			"#include <sstream> \n"
			"#include\"/home/matej/Downloads/Cuba-4.2.2/cuba.h\" \n \n";

	writer << "#define NDIM " << _integVars.size() << "\n"
			"#define NCOMP 1 \n"
			"#define USERDATA NULL \n"
			"#define NVEC 1 \n"
			"#define EPSREL 1e-7 \n"
			"#define EPSABS 1e-12 \n"
			"#define VERBOSE 2 \n"
			"#define LAST 4 \n"
			"#define SEED 0 \n"
			"#define MINEVAL 0 \n"
			"#define MAXEVAL 10000000 \n"
			"\n"
			"#define NSTART 1000 \n"
			"#define NINCREASE 500 \n"
			"#define NBATCH 1000 \n"
			"#define GRIDNO 0 \n"
			"#define STATEFILE NULL \n"
			"#define SPIN NULL \n"
			"\n"
			"#define KEY 0 \n";

	writer
			<< "static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) { \n \n ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "#define " << _integVars.at(j) << " xx[" << j << "] \n";
	}

	writer << "#define f ff[0] \n \n";

	writer << "f = " << coef << "; \n \n";
	//writer << "f = 1./4*t0*pow((1+2*t5*t0*pow(t1,2)+6*t4*pow(t2,2)*t0+8*t2*t5*t0*t1+t0*t1+2*t4+8*t4*t2*t0*t1+2*t2+2*t2*t0+6*pow(t2,2)*t5*t0+t0*pow(t1,2)+2*t5*t0*t1+2*t5+4*t4*t2+2*t4*t0*pow(t1,2)+4*t2*t0*t1+2*t5*t1+3*pow(t2,2)*t0+t1+4*t4*t2*t0+2*t4*t1+4*t2*t5*t0+2*t4*t0*t1+4*t2*t5),(-3))" << "; \n \n";

	writer << " return 0; \n"
			"} \n";

	writer
			<< "int main() { \n"
					"clock_t begin = clock(); \n"
					"int comp, nregions, neval, fail; cubareal integral[NCOMP], error[NCOMP], prob[NCOMP]; \n"
					"\n"
					"Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, \n"
					"   EPSREL, EPSABS, VERBOSE, SEED, \n"
					"	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, \n"
					"    GRIDNO, STATEFILE, SPIN, \n"
					"   &neval, &fail, integral, error, prob); "
					"\n"
					"printf(\"VEGAS RESULT:\\tneval %d\\tfail %d\\n\", neval, fail); \n"
					"for( comp = 0; comp < NCOMP; ++comp ) \n"
					"			printf(\"VEGAS RESULT:\\t%.8f +- %.8f\\tp = %.3f\\n\", "
					"(double)integral[comp], (double)error[comp], (double)prob[comp]); \n \n"
					"\n"
					"clock_t end = clock(); \n"
					"std::cout << \"Elapsed time: \" << double(end - begin) / CLOCKS_PER_SEC"
					"<<\"s.\"<< std::endl; \n"
					"std::stringstream s; \n"
					"std::fstream writer;\n";
	writer << "s<<\"" << _path << "Diagram_" << _diag.getName()
			<< "_finite_eps0" << ".out\"; \n";
	writer << "writer.open(s.str(), std::fstream::out); \n"
			"writer<<std::fixed; \n"
			"writer.precision(10); \n";

	writer << "writer<<" << overallFactor << "<< \" \\n \"; \n";
	writer
			<< "for (comp = 0; comp < NCOMP; ++comp){ \n"
					"	writer<<(double) integral[comp]<<\"\t\"<<(double) error[comp]<<\"\t\"<<(double) prob[comp]; \n"
					"}; \n"
					"return 0; \n} "
					"\n";
	//-------------------------------------------------------------------

	writer.close();
	return;
}

void writeFiniteCVegas_eps1(Diagram _diag,
		std::vector<GiNaC::symbol> _integVars,
		std::vector<GiNaC::ex> _finitePartCoefs, GiNaC::ex _overallFactor,
		std::string _path) {

	std::stringstream aux;
	std::string coef = "";

	std::stringstream a;
	// the *eps^2 is to make sure the result including overall factor is really order eps^1
	a << "\"" << _overallFactor << "*eps^2\"";
	std::string overallFactor = a.str();

	//int currentPowerOfEps;

	std::fstream writer;

	//-------------------------------------------------------------------
	// eps^{-1}
	//-------------------------------------------------------------------

	//currentPowerOfEps = 1;
	_finitePartCoefs.at(1).print(GiNaC::print_csrc(aux));
	coef = aux.str();
	// write záhlavok
	// TODO

	std::stringstream s;
	s << _path << "Diagram_" << _diag.getName() << "_finite_eps1" << ".cpp";

	writer.open(s.str(), std::fstream::out);

	writer << "#include <stdio.h> \n"
			"#include <stdlib.h> \n"
			"#include <math.h> \n"
			"#include <iostream> \n"
			"#include <fstream> \n"
			"#include <sstream> \n"
			"#include\"/home/matej/Downloads/Cuba-4.2.2/cuba.h\" \n \n";

	writer << "#define NDIM " << _integVars.size() << "\n"
			"#define NCOMP 1 \n"
			"#define USERDATA NULL \n"
			"#define NVEC 1 \n"
			"#define EPSREL 1e-7 \n"
			"#define EPSABS 1e-12 \n"
			"#define VERBOSE 2 \n"
			"#define LAST 4 \n"
			"#define SEED 0 \n"
			"#define MINEVAL 0 \n"
			"#define MAXEVAL 10000000 \n"
			"\n"
			"#define NSTART 1000 \n"
			"#define NINCREASE 500 \n"
			"#define NBATCH 1000 \n"
			"#define GRIDNO 0 \n"
			"#define STATEFILE NULL \n"
			"#define SPIN NULL \n"
			"\n"
			"#define KEY 0 \n";

	writer
			<< "static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) { \n \n ";

	for (int j = 0; j < _integVars.size(); j++) {
		writer << "#define " << _integVars.at(j) << " xx[" << j << "] \n";
	}

	writer << "#define f ff[0] \n \n";

	writer << "f = " << coef << "; \n \n";
	//writer << "f = 1./4*t0*pow((1+2*t5*t0*pow(t1,2)+6*t4*pow(t2,2)*t0+8*t2*t5*t0*t1+t0*t1+2*t4+8*t4*t2*t0*t1+2*t2+2*t2*t0+6*pow(t2,2)*t5*t0+t0*pow(t1,2)+2*t5*t0*t1+2*t5+4*t4*t2+2*t4*t0*pow(t1,2)+4*t2*t0*t1+2*t5*t1+3*pow(t2,2)*t0+t1+4*t4*t2*t0+2*t4*t1+4*t2*t5*t0+2*t4*t0*t1+4*t2*t5),(-3))" << "; \n \n";

	writer << " return 0; \n"
			"} \n";

	writer
			<< "int main() { \n"
					"clock_t begin = clock(); \n"
					"int comp, nregions, neval, fail; cubareal integral[NCOMP], error[NCOMP], prob[NCOMP]; \n"
					"\n"
					"Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, \n"
					"   EPSREL, EPSABS, VERBOSE, SEED, \n"
					"	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, \n"
					"    GRIDNO, STATEFILE, SPIN, \n"
					"   &neval, &fail, integral, error, prob); "
					"\n"
					"printf(\"VEGAS RESULT:\\tneval %d\\tfail %d\\n\", neval, fail); \n"
					"for( comp = 0; comp < NCOMP; ++comp ) \n"
					"			printf(\"VEGAS RESULT:\\t%.8f +- %.8f\\tp = %.3f\\n\", "
					"(double)integral[comp], (double)error[comp], (double)prob[comp]); \n \n"
					"\n"
					"clock_t end = clock(); \n"
					"std::cout << \"Elapsed time: \" << double(end - begin) / CLOCKS_PER_SEC"
					"<<\"s.\"<< std::endl; \n"
					"std::stringstream s; \n"
					"std::fstream writer;\n";
	writer << "s<<\"" << _path << "Diagram_" << _diag.getName()
			<< "_finite_eps1" << ".out\"; \n";
	writer << "writer.open(s.str(), std::fstream::out); \n"
			"writer<<std::fixed; \n"
			"writer.precision(10); \n";

	writer << "writer<<" << overallFactor << "<< \" \\n \"; \n";
	writer
			<< "for (comp = 0; comp < NCOMP; ++comp){ \n"
					"	writer<<(double) integral[comp]<<\"\t\"<<(double) error[comp]<<\"\t\"<<(double) prob[comp]; \n"
					"}; \n"
					"return 0; \n} "
					"\n";
	//-------------------------------------------------------------------

	writer.close();
	return;
}

//===========================================================================

