/*
 * TimeOrderingsAndTimeCuts.hpp
 *
 *  Created on: Aug 10, 2023
 *      Author: matej
 */

#ifndef TIMEORDERINGSANDTIMECUTS_GINAC_HPP_
#define TIMEORDERINGSANDTIMECUTS_GINAC_HPP_

#include<iostream>
#include<fstream>
#include "ElementsAndInput_ginac.hpp"

//===========================================================================
// code related to finding time orderings of diagram
// Implementation in TimeOrderingsAndTimeCuts.cpp

void permute(std::vector<Vertex> _vertices, int l, int r,
		std::vector<std::vector<Vertex>> &permutations,
		std::vector<Propagator> _propagators);

bool isAllowedPermutation(std::vector<Vertex> _permutation,
		std::vector<Propagator> _propagators);

std::vector<std::vector<Vertex>> makeAllAllowedPermutationsOfVertices(
		std::vector<Vertex> _vertices, std::vector<Propagator> _propagators);

int findPositionOfVertexInOrdering(std::vector<Vertex> _permutation,
		Vertex _vert);

void throwOutForbiddenPermutations(
		std::vector<std::vector<Vertex>> &_allPermutations,
		std::vector<Propagator> _propagators);

void writeTimeOrderings(Diagram _diag,
		std::vector<std::vector<Vertex>> _orderings);

//===========================================================================

//===========================================================================
// code related to performing time cuts and making integrand expression after time integration

std::vector<Propagator> getPropagatorsWithEndVertex(Vertex _endVertex,
		std::vector<Propagator> _propagsToSearch);
std::vector<Propagator> getPropagatorsWithStartVertex(Vertex _startVertex,
		std::vector<Propagator> _propagsToSearch);
std::vector<Propagator> getPropagatorsStartingSooner(
		std::vector<Propagator> _currentPropagators,
		std::vector<Vertex> _timeOrdering, Vertex _startPoint);

int countTaus(std::vector<std::vector<Propagator>> _allCuts, int _cut);

void appendVectors(std::vector<Propagator> &_a, std::vector<Propagator> _b);
void appendVectors(std::vector<Vertex> &_a, std::vector<Vertex> _b);

std::vector<std::vector<Propagator>> getCuts(std::vector<Vertex> _timeOrdering,
		Diagram _diag);

GiNaC::ex makeIntegrandCorespToCuts(
		std::vector<std::vector<Propagator>> _cutsForGivenTimeOrdering,
		GiNaC::symbol &_tau);

void writeIntegrandAfterTimeCuts(Diagram _diag, GiNaC::ex _integrand);

void updateFactorFromVertsAndPropsAfterCuts(Diagram &_diag);
//===========================================================================

//===========================================================================
// some code variants useful in Feynman parametrization
int countTausInCut(std::vector<Propagator> _cut);
GiNaC::ex findDenominatorForCut(std::vector<Propagator> _cut,
		GiNaC::symbol &_tau);
//===========================================================================

#endif /* TIMEORDERINGSANDTIMECUTS_GINAC_HPP_ */
