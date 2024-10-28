//============================================================================
// Name        : IPLoop.cpp
// Author      : M. Kecer
// Version     : 12.1
// Description : Code performs 3-loop calculations in dynamic isotropic percolation
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
 * this code cannot deal with non-causal props. as of yet (tadpoles)
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

	diags = {T51, T52, T53, T54}; // all

	for(int i = 0; i < diags.size(); i++){

		currentPath = PATH + "T5_" + topologyT5 + "/" + diags.at(i).getName() + "/";
		findDivergentPartsPropToExtMom_3loop(diags.at(i), currentPath);

	}

	std::cout << "Done!";
	clock_t end = clock();
	std::cout << "Elapsed time: " << double(end - begin) / CLOCKS_PER_SEC
			<< "s." << std::endl;
	return 0;
}
