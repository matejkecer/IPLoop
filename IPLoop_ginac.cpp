//============================================================================
// Name        : IPLoop.cpp
// Author      : M. Kecer
// Version     : 1.2
// Description : symmetry factor of diagram now works
//============================================================================

#include "ElementsAndInput_ginac.hpp"
#include "TimeOrderingsAndTimeCuts_ginac.hpp"
#include "FeynmanAndAlphaParam.hpp"
#include "SectorDecomposition.hpp"
#include "OutputWriter.hpp"
#include<iostream>
#include<sstream>
#include<vector>
#include <ginac/ginac.h>

/* Limitations:
 * 1.) the code is designed for IP process -> phi3 topology, type of propagator expressions, cuts, causal propags. only,
 * this code cannot deal with non-causal props. as of yet (closed loops and tadpoles)
 * 2.) the code breaks down if there is vertex indexed with two-digit number in nickel indices (say e...|....|e10)
 * where the last thing says e and 10 as vertices -> input from nickel breaks down
 * 3.) In spanning trees function  makeAllPossibleNTuplesOfProps(...) I have only given implementation
 * for at most triples -> only three loop graphs, if more are looked at 4-tuple needs to be added and so on
 * 4.) In multiple parts of code I only consider working with 2 point functions (main goal of this project
 * is to calculate exponent "z" which can be found by considering 2 point functions, namely their parts that
 * are proportional to external frequency.) If problem is encountered - console should std::cout a message about this.
 * Implementation of the functionality for 3 point function can be straightforwardly created following general procedures
 * for 2 point functions already working. This is to be done if such situation should arise.
 * 5.) In some parts of code related to feynman parametrization we specifically assume, that we are calculating
 * part of diagrams proportional to external frequency in IP (ext. momenta set to zero). If this is not the case
 * some of the code would need to be expanded for such cases. - usually std::cout message should be displayed.
 * */


int main() {
	std::string PATH = "/home/matej/eclipse-workspace/IPLoop_ginac/Part_propto_tau/3-loop/";
	std::string currentPath;

	clock_t begin = clock();

	std::vector<Diagram> diags;
	std::string topologyT5 = "e12|34|35|e|55||:";
	Diagram T51("e12|34|35|e|55||:0P_mP_pP|mP_pP|pP_Pp|0p|mP_pP||", "y", "n", 1); //d1
	Diagram T52("e12|34|35|e|55||:0P_mP_pP|pP_mP|pP_Pp|0p|mP_pP||", "y", "n", 1); //d2
	Diagram T53("e12|34|35|e|55||:0P_mP_pP|pP_Pp|pP_mP|0p|Pm_Pp||", "y", "n", 1); //d3
	Diagram T54("e12|34|35|e|55||:0P_mP_pP|pP_Pp|mP_pP|0p|Pm_Pp||", "y", "n", 1); //d4

	//diags = {T51, T52, T53, T54}; //všetky
	diags = {T53, T52, T51}; // 1 časť
	//diags = {T54}; // 5 častí
	for(int i = 0; i < diags.size(); i++){

		currentPath = PATH + "T5_" + topologyT5 + "/" + diags.at(i).getName() + "/";
		//fs::create_directory(currentPath);
		findDivergentPartsPropToExtMom_3loop(diags.at(i), currentPath);

	}

	std::cout << "Done!";
	clock_t end = clock();
	std::cout << "Elapsed time: " << double(end - begin) / CLOCKS_PER_SEC
			<< "s." << std::endl;
	return 0;
}
