//============================================================================
// Name        : IPLoop.cpp
// Author      : M. Kecer
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



int main() {
	// for CPU Cuba Vegas output
	//std::string CUBAPATH = "/home/matej/Documents/Cuba-4.2.2";
	//std::string PATH = "/home/matej/eclipse-workspace/IPLoop_ginac/2-loop/";

	// for GPU cuVegas output
	std::string PATH = "/home/matej/eclipse-workspace/IPLoop_ginac/2-loop_cuVegas/";

	std::string currentPath;

	clock_t begin = clock();

	std::vector<Diagram> diags;

	Diagram D2L1("e12|e3|33|:0P_mP_pP|0p_Pp|mP_pP|", "y", "n", 2);

	diags = {D2L1};

	for(int i = 0; i < diags.size(); i++){

		currentPath = PATH + "/" + diags.at(i).getName() + "/";
		//findDivergentPartsPropToExtFreq_2loop_timeOrderingwise(diags.at(i), currentPath, CUBAPATH);
		findDivergentPartsPropToExtFreq_2loop_freqPartwise_GPU(diags.at(i), currentPath);
	}

	std::cout << "Done!";
	clock_t end = clock();
	std::cout << "Elapsed time: " << double(end - begin) / CLOCKS_PER_SEC
			<< "s." << std::endl;
	return 0;
}

