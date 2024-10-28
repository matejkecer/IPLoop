/*
 * TimeOrderingsAndTimeCuts.cpp
 *
 *  Created on: Aug 10, 2023
 *      Author: matej
 */

#include"TimeOrderingsAndTimeCuts_ginac.hpp"
#include "FeynmanAndAlphaParam.hpp"
#include <ginac/ginac.h>

//===========================================================================
// code related to finding time orderings of diagram - implementation

void permute(std::vector<Vertex> _vertices, int l, int r,
		std::vector<std::vector<Vertex>> &_permutations,
		std::vector<Propagator> _propagators) {
	// Subprogram to make all permutations of elements of std::vector<Vertex> (vertices)
	// makes all possible permutations elements of std::vector, permuting from elem. at index l
	// to elem. at index r, then it saves the permutations into vector _permutations
	// Auxiliary program used in function "makeAllPermutationsOfVertices"

	if (l == r) {
		//cout << vertices << endl;
		if (isAllowedPermutation(_vertices, _propagators)) {
			_permutations.push_back(_vertices);
			//std::cout<<_permutations.size()<<std::endl;
		}
	} else {
		// Permutations made
		for (int i = l; i <= r; i++) {
			// Swapping done
			std::swap(_vertices.at(l), _vertices.at(i));

			// Recursion called
			permute(_vertices, l + 1, r, _permutations, _propagators);

			//backtrack
			std::swap(_vertices.at(l), _vertices.at(i));
		}
	}
	return;
}

bool isAllowedPermutation(std::vector<Vertex> _permutation,
		std::vector<Propagator> _propagators) {

	Vertex isLaterInPropag;
	Vertex isSoonerInPropag;
	int indexOfLaterInOrdering; // index of vertex-that-is-later-in-propag inside given ordering
	int indexOfSoonerInOrdering;

	for (int k = 0; k < _propagators.size(); k++) {
		// go through all propagators in set of propagators

		//determine, which letter of vertex is later in time in causal propagator
		isLaterInPropag = _propagators.at(k).getEndVert();
		isSoonerInPropag = _propagators.at(k).getStartVert();

		// find positions of given vertices in given ordering
		indexOfLaterInOrdering = findPositionOfVertexInOrdering(_permutation,
				isLaterInPropag);
		indexOfSoonerInOrdering = findPositionOfVertexInOrdering(_permutation,
				isSoonerInPropag);

		// if the relative position is wrong, throw out the ordering from results
		// NOTE: This works only for causal propagators - we don't allow closed loops
		// i.e. propagators type {A,A,momentum}
		if (!(indexOfLaterInOrdering < indexOfSoonerInOrdering)) {
			return false;
		}
	}

	return true;
}

std::vector<std::vector<Vertex>> makeAllAllowedPermutationsOfVertices(
		std::vector<Vertex> _vertices, std::vector<Propagator> _propagators) {
	// Makes all permutations of vertices and saves them into vector "permutations", uses function permute
	// by default makes all the possible permutations of all the vertices, from first to last

	int n = _vertices.size();
	std::vector<std::vector<Vertex>> permutations;

	permute(_vertices, 0, n - 1, permutations, _propagators);

	return permutations;
}

int findPositionOfVertexInOrdering(std::vector<Vertex> _permutation,
		Vertex _vert) {
	// function determines index of a particular vertex _vert within a time ordering (
	// some particular permutation of all vertices) _permutation
	// function returns this index
	// If the vertex is not there, it returns -1
	// Auxiliary function used in "throwOutForbiddenPermutations"

	int indexOfVertex = -1;
	std::string s;

	for (int j = 0; j < _permutation.size(); j++) {
		s = _permutation.at(j).getVertName();
		if (s.compare(_vert.getVertName()) == 0) {
			indexOfVertex = j;
		}
	}
	return indexOfVertex;
}

void throwOutForbiddenPermutations(
		std::vector<std::vector<Vertex>> &_allPermutations,
		std::vector<Propagator> _propagators) {
	// TODO test!!!!!

	// function checks all the permutations of vertices and gradually removes those,
	// which don't respect causality of propagators in _propagators
	// At the end, only allowed orderings remain in container allowedOrderings -> which is
	// what replaces _allPermutations after the function finishes the job
	// ASSUMES that _allPermutations contains all the possible permutations before the function is executed
	//         (those that are not allowed are gradually thrown out)

	// initialization of variables
	Vertex isLaterInPropag;
	Vertex isSoonerInPropag;
	int indexOfLaterInOrdering;	// index of vertex-that-is-later-in-propag inside given ordering
	int indexOfSoonerInOrdering;
	std::vector<Vertex> currentOrdering;

	std::vector<std::vector<Vertex>> allowedOrderings = _allPermutations;
	std::vector<std::vector<Vertex>> aux;	// auxiliary container used inside

	for (int k = 0; k < _propagators.size(); k++) {
		// go through all propagators in set of propagators

		//determine, which letter of vertex is later in time in causal propagator
		isLaterInPropag = _propagators.at(k).getEndVert();
		isSoonerInPropag = _propagators.at(k).getStartVert();

		for (int i = 0; i < allowedOrderings.size(); i++) {
			// go over all orderings
			currentOrdering = allowedOrderings.at(i);

			// find positions of given vertices in given ordering
			indexOfLaterInOrdering = findPositionOfVertexInOrdering(
					currentOrdering, isLaterInPropag);
			indexOfSoonerInOrdering = findPositionOfVertexInOrdering(
					currentOrdering, isSoonerInPropag);

			// if the relative position is wrong, throw out the ordering from results
			// NOTE: This works only for causal propagators - we don't allow closed loops
			// i.e. propagators type {A,A,momentum}
			if (indexOfLaterInOrdering < indexOfSoonerInOrdering) {
				aux.push_back(currentOrdering);
				//continue; // TODO - check it but i think it is just useless to put it here, it does nothing
			}

		}

		// now allowedOrderings holds all the orderings that are consistent with causal structure
		// of propagator _propagators.at(k)
		// we still have to check other propagators
		allowedOrderings.clear();
		allowedOrderings = aux;
		aux.clear();
	}
	_allPermutations.clear();
	_allPermutations = allowedOrderings;
	return;
}

std::vector<std::vector<Vertex>> Diagram::findAllPossibleTimeOrderings(
		std::string _writeData) {
	// Argument - writeData to output file? - choices "y" or "n"
	// Return - std::vector of all possible time orderings of vertices of diagram
	//
	// Function lists all possible time orderings of a diagram, first it makes all the possible permutations
	// of vertices, then filters only "allowed" orderings (ones that don't screw up the causality of propagators)

	// initializations
	std::vector<std::vector<Vertex>> allOrderings;

	// find all possible permutations
	allOrderings = makeAllAllowedPermutationsOfVertices(this->vertices, this->intPropags);

	// throws out all the time orderings that are not allowed
	//throwOutForbiddenPermutations(allOrderings, this->intPropags);

	// prints the number of allowed time orderings for the input diagram
	//std::cout << allOrderings.size() << std::endl;

	if (_writeData.compare("y") == 0) {
		writeTimeOrderings(*this, allOrderings);
	} else {
		if (_writeData.compare("n") == 0) {
			return allOrderings;
		} else {
			std::cout
					<< "Wrong input. Argument is either \"y\" or \"n\". The data was not written to output file.";
		}
	}

	return allOrderings;
}

void writeTimeOrderings(Diagram _diag,
		std::vector<std::vector<Vertex>> _orderings) {
	// Arguments - _diagram whose orderings we are writing, _orderings of the diagram
	// to be printed, file is written in the current folder (where project is)

	// Printwriter to write possible time-ordering data into output file located in project folder

	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << "Diagram_" << _diag.getName() << ".txt";

	std::fstream writer;

	//TODO - try catch
	writer.open(s.str(), std::fstream::app);

	writer << "Diagram: " << _diag.getName() << "\n";

	writer << "Number of Loops: " << _diag.getNumOfLoops();
	writer << "\n";

	writer << "Vertices: ";
	for (int i = 0; i < _diag.getVertices().size(); i++) {
		writer << _diag.getVertices().at(i).info();
		if (i != _diag.getVertices().size() - 1) {
			writer << ", ";
		}
	}

	writer << "\n";

	writer << "External point designation: ";
	writer << _diag.getExtPoint().info();
	writer << "\n";

	writer << "External propagators: ";
	for (int i = 0; i < _diag.getExtPropags().size(); i++) {
		writer << _diag.getExtPropags().at(i).info();
		if (i != _diag.getExtPropags().size() - 1) {
			writer << ", ";
		}
	}

	writer << "\n";

	writer << "Internal propagators: ";
	for (int i = 0; i < _diag.getIntPropags().size(); i++) {
		writer << _diag.getIntPropags().at(i).info();
		if (i != _diag.getIntPropags().size() - 1) {
			writer << ", ";
		}
	}

	writer << "\n";

	writer << "Internal momenta: {";
	for (int i = 0; i < _diag.getIntMomenta().size(); i++) {
		writer << _diag.getIntMomenta().at(i);
		if (i != _diag.getIntMomenta().size() - 1) {
			writer << ", ";
		}
	}
	writer << "}" << "\n";

	writer << "External momenta: {";
	for (int i = 0; i < _diag.getExtMomenta().size(); i++) {
		writer << _diag.getExtMomenta().at(i);
		if (i != _diag.getExtMomenta().size() - 1) {
			writer << ", ";
		}
	}
	writer << "}" << "\n";
	writer << "\n";

	writer << "Number of time-orderings: " << _orderings.size() << '\n' << '\n';

	for (int i = 0; i < _orderings.size(); i++) {
		for (int j = 0; j < _orderings.at(i).size(); j++) {
			writer << _orderings.at(i).at(j).getVertName();
		}
		writer << "\n";
	}

	writer << '\n';
	writer << "---------------------------------------------------------";
	writer << '\n';
	writer.close();
}

//===========================================================================

//===========================================================================
// code related to performing time cuts and making integrand expression after time integration - implementation
std::vector<Propagator> getPropagatorsWithEndVertex(Vertex _endVertex,
		std::vector<Propagator> _propagsToSearch) {
	// Auxilary function used in makeIntegrand function
	// I give it endVertex, say A, and it finds me all such propagators, that end with A (e.g. AX, AE, ... but not e.g. XA)
	// (end in a sense, that A is later in time in causal propagator)
	std::vector<Propagator> result;

	for (int i = 0; i < _propagsToSearch.size(); i++) {
		// go through all propagators
		if (_propagsToSearch.at(i).getEndVert().getVertName().compare(
				_endVertex.getVertName()) == 0) {
			// if the propagator starts with the right vertex, then add it to currentPropagators
			result.push_back(_propagsToSearch.at(i));
		}
	}
	return result;
}

std::vector<Propagator> getPropagatorsWithStartVertex(Vertex _startVertex,
		std::vector<Propagator> _propagsToSearch) {
	// Auxilary function used in makeIntegrand function
	// I give it endVertex, say A, and it finds me all such propagators, that end with A (e.g. AX, AE, ... but not e.g. XA)
	// (end in a sense, that A is later in time in causal propagator)
	std::vector<Propagator> result;

	for (int i = 0; i < _propagsToSearch.size(); i++) {
		// go through all propagators
		if (_propagsToSearch.at(i).getStartVert().getVertName().compare(
				_startVertex.getVertName()) == 0) {
			// if the propagator starts with the right vertex, then add it to currentPropagators
			result.push_back(_propagsToSearch.at(i));
		}
	}
	return result;
}

std::vector<Propagator> getPropagatorsStartingSooner(
		std::vector<Propagator> _propagsToSearch,
		std::vector<Vertex> _timeOrdering, Vertex _startPoint) {
	// Auxiliary function used in makeIntegrand. It checks if the propagators startVertex (the vertex in second place - sooner in time)
	// is sooner in time than startpoint of current interval where cut is being made (sooner or exactly it)
	//
	// example - propag {A,C, mom}, and cut being made at interval BC, inside ordering ABCDE -> index of C (from propag)
	// inside ordering is >= index of C (startpoint of interval), so this propagator will be selected
	// also {A,D,mom}, {B,E,mom} even {D,E,mom} would work and be selected
	//
	// We will be trying to find propagators that contribute to given cut in later functions,
	// things like {D,E} will be taken care of and filtered later (not in this function)

	int indexToBeCompared;
	std::vector<Propagator> result;

	for (int i = 0; i < _propagsToSearch.size(); i++) {
		// go through all propagators

		// find index of propagator startPoint within the given time ordering
		indexToBeCompared = findPositionOfVertexInOrdering(_timeOrdering,
				_propagsToSearch.at(i).getStartVert());

		if (indexToBeCompared
				>= findPositionOfVertexInOrdering(_timeOrdering, _startPoint)) {
			// if the index of propagator endpoint is >= index of interval endpoint,
			// then it contributes to the cut through given interval
			result.push_back(_propagsToSearch.at(i));
		}
	}
	return result;
}

std::string exToString(GiNaC::ex _ex) {
	// it is a bit of a cheat, but anyway
	std::stringstream a;
	a << _ex;

	return a.str();
}

int countTaus(std::vector<std::vector<Propagator>> _allCuts, int _cut) {
	// TODO - try to find out if comparing expressions wouldnt work -> don't use exToString
	// Auxiliary function to be used in wolframFormatIntegrand
	// is needed because if I simply put numOfMomenta in cut* tau I would have problem,
	// because it would count also momentum 0

	// Arguments: _allCuts is vector<vector<Propagator>> (one cut is vector<Propagator>
	// and it holds propagators contributing to given cut),
	// argument _cut is index of some particular cut in which I count taus

	// initialize as number of momenta that contribute to cut
	int numberOfTaus = _allCuts.at(_cut).size();

	// now subtract instances where the momenta are 0
	for (int i = 0; i < _allCuts.at(_cut).size(); i++) {
		if (exToString(_allCuts.at(_cut).at(i).getMomentum()).compare("0")
				== 0) {
			numberOfTaus -= 1;
		}
	}

	return numberOfTaus;
}

void appendVectors(std::vector<Propagator> &_a, std::vector<Propagator> _b) {
	// Aux function
	// function takes all the components of vector _b and appends it to _a
	// as a push_back basically
	// new _a is {old _a, _b}

	for (int i = 0; i < _b.size(); i++) {
		_a.push_back(_b.at(i));
	}

	return;
}

void appendVectors(std::vector<Vertex> &_a, std::vector<Vertex> _b) {
	// Aux function
	// function takes all the components of vector _b and appends it to _a
	// as a push_back basically
	// new _a is {old _a, _b}

	for (int i = 0; i < _b.size(); i++) {
		_a.push_back(_b.at(i));
	}

	return;
}

std::vector<std::vector<Propagator>> getCuts(std::vector<Vertex> _timeOrdering,
		Diagram _diag) {
	// Function finds and returns all the time cuts for given time ordering (we will use them to
	// build integrand later).


	// return vector<cuts> -> holding all the time cuts for given time ordering
	// one cut is type vector<Propagator> and holds propagators that contribute to the cut
	// through some particular interval within time ordering

	// argument _timeOrdering - ordering for which we are making cuts
	// argument _diag - diagram in question

	// initialize variables
	Vertex endPoint;
	Vertex startPoint;
	Vertex currentPoint;
	std::vector<Propagator> currentPropagators;
	std::vector<Propagator> contributingPropagsToCut = { };
	std::vector<Propagator> cut;
	std::vector<std::vector<Propagator>> allCuts;

	for (int i = 0; i < _timeOrdering.size() - 1; i++) {
		// for all intervals in given time ordering - do cut

		// define your interval
		endPoint = _timeOrdering.at(i);
		startPoint = _timeOrdering.at(i + 1);

		//a) -> look forward, find propagators that have endpoint timeOrdering.at(i)
		// say interval AX, we are looking for propagators type {"A", "sth"}
		// these propags. will be saved to currentPropagators
		currentPropagators = getPropagatorsWithEndVertex(endPoint,
				_diag.getIntPropags());

		// check if these propagators contribute to cuts through given interval
		// basically now if we are in interval AX and check the propag. {"A", "sth"}, we check
		// if in the particular time ordering "sth" is sooner in time then "X" (or "X" itself).
		// such propagators contribute and get saved to aux
		contributingPropagsToCut = getPropagatorsStartingSooner(
				currentPropagators, _timeOrdering, startPoint);
		appendVectors(cut, contributingPropagsToCut);
		//b) -> look behind, for all indices before ith, repeat the same procedure as in case a)
		// say we have ordering "AXBYCDZEF", we are in interval BY, we already checked a) - propags. that start with B
		// now we move on to X, we ask which propags. are of type {"X", "sth"}, and if "sth" is later/farther than Y (or Y itself)
		// then we move on to A and repeat the procedure

		for (int j = 0; j < i; j++) {
			// clear currentPropagators, so that we didn't account for the same thing multiple times
			currentPropagators.clear();
			contributingPropagsToCut.clear();

			// startpoint is now vertex at index j
			currentPoint = _timeOrdering.at(j);

			//rest of the procedure such as in a)
			currentPropagators = getPropagatorsWithEndVertex(currentPoint,
					_diag.getIntPropags());
			contributingPropagsToCut = getPropagatorsStartingSooner(
					currentPropagators, _timeOrdering, startPoint);
			appendVectors(cut, contributingPropagsToCut);
		}

		// again, clear the vectors in order not to count propagators in this iteration also in the next iteration
		currentPropagators.clear();
		contributingPropagsToCut.clear();

		// now vector cut contains all the propags contributing to cut in the given interval
		// add it to allCuts
		allCuts.push_back(cut);

		// at last, clear the current cut before moving on to next interval
		cut.clear();
	}

	// Now we have vector<vector<Propagator>> allCuts that contains all the cuts for given time ordering
	// (that means cuts through all the intervals inside the ordering)

	return allCuts;
}

GiNaC::ex makeIntegrandCorespToCuts(
		std::vector<std::vector<Propagator>> _cutsForGivenTimeOrdering,
		GiNaC::symbol &_tau) {
	// Function returns ginac expression for integrand after performing time cuts on given
	// time ordering.
	// Argument - cuts corresponding to given ordering (result of "getCuts" function)
	// return - ginac expression for integrand

	GiNaC::ex result = 1;
	GiNaC::ex mom = 0;
	int numOfTaus = 0;

	// go through all intervals of time ordering
	for (int i = 0; i < _cutsForGivenTimeOrdering.size(); i++) {

		// go through all propagators contributing to cut through given interval
		for (int j = 0; j < _cutsForGivenTimeOrdering.at(i).size(); j++) {
			mom += GiNaC::pow(
					_cutsForGivenTimeOrdering.at(i).at(j).getMomentum(), 2);
		}
		//std::cout<<mom<<"\n";
		numOfTaus = countTaus(_cutsForGivenTimeOrdering, i);
		result *= 1 / (mom + numOfTaus * _tau);

		// reinitialize mom to zero
		mom = 0;
	}
	return result;
}

GiNaC::ex Diagram::makeIntegrandAfterTimeCuts(
		std::vector<std::vector<Vertex>> _allTimeOrderings,
		std::string _writeData) {
	// Function to find and return ginac expression for integrand of the diagram after time cuts
	// Arguments: _allTimeOrderings - time orderings of the corresponding diagram
	//			  _tau - parameter that every propagator in IP theory contains
	//			  _writeData to output file? - choices "y" or "n"
	// Return - ginac expression for integranda fter time cuts

	GiNaC::ex result = 0;
	std::vector<std::vector<Propagator>> cutsForTheOrdering;

	for (int k = 0; k < _allTimeOrderings.size(); k++) {
		cutsForTheOrdering = getCuts(_allTimeOrderings.at(k), *this);
		result = result
				+ makeIntegrandCorespToCuts(cutsForTheOrdering, this->tau);
	}

	if (_writeData.compare("y") == 0) {
		writeIntegrandAfterTimeCuts(*this, result);
	} else {
		if (_writeData.compare("n") == 0) {
			return result;
		} else {
			std::cout
					<< "Wrong input. Argument is either \"y\" or \"n\". The data was not written to output file.";
		}
	}

	return result;
}

void writeIntegrandAfterTimeCuts(Diagram _diag, GiNaC::ex _integrand) {
	// Arguments - _diagram whose integrand we are writing, _integrand of the diagram after time cuts
	// to be printed, file is written in the current folder (where project is)

	// Printwriter to write integrand after time cuts into output file located in project folder

	std::stringstream s;

	//s << "/home/matej/eclipse-workspace/IPLoop/src/" << "Diagram_" << _diag.getName() <<".txt";
	s << "Diagram_" << _diag.getName() << ".txt";

	std::fstream writer;

	//TODO - try catch
	writer.open(s.str(), std::fstream::app);

	writer << "Integrand: " << "\n";

	writer << _integrand;

	writer << '\n';
	writer << "---------------------------------------------------------";
	writer << '\n';
	writer.close();
}

//===========================================================================

//===========================================================================
// some code variants useful in Feynman parametrization

int countTausInCut(std::vector<Propagator> _cut) {
	// Auxiliary function to be used in wolframFormatIntegrand
	// is needed because if I simply put numOfMomenta in cut* tau I would have problem,
	// because it would count also momentum 0

	// Arguments:  cut is vector<Propagator> and it holds propagators contributing
	//to given cut), where i count taus

	// initialize as number of momenta that contribute to cut
	int numberOfTaus = _cut.size();

	// now subtract instances where the momenta are 0
	for (int i = 0; i < _cut.size(); i++) {
		//if (exToString(_cut.at(i).getMomentum()).compare("0")== 0) {
		if (_cut.at(i).getMomentum() == 0) {
			numberOfTaus -= 1;
		}
	}

	return numberOfTaus;
}

GiNaC::ex findDenominatorForCut(std::vector<Propagator> _cut,
		GiNaC::symbol &_tau) {

	GiNaC::ex result = 1;
	GiNaC::ex mom = 0;
	int numOfTaus = 0;

	// go through all propagators contributing to cut _cut
	for (int j = 0; j < _cut.size(); j++) {
		mom += GiNaC::pow(_cut.at(j).getMomentum(), 2);
	}
	//std::cout<<mom<<"\n";
	numOfTaus = countTausInCut(_cut);
	result = (mom + numOfTaus * _tau);

	// reinitialize mom to zero
	mom = 0;

	return result;
}
//===========================================================================

