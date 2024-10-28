/*
 * ElementsAndInput.hpp
 *
 *  Created on: Aug 7, 2023
 *      Author: M. Kecer
 */

#ifndef ELEMENTSANDINPUT_GINAC_HPP_
#define ELEMENTSANDINPUT_GINAC_HPP_

#include<iostream>
#include<sstream>
#include<vector>
#include <ginac/ginac.h>

//===========================================================================
// class Vertex
class Vertex {
	/* class vertex holds vertices (coded by name of vertex, e.g. "A").
	 * TODO - implement type of propagator.
	 * */

protected:
	// attributes
	std::string vertName;
	std::string vertType;
	// in IP there are 3 types of fields P (psi'), p (psi), m(phi with derivative)
	// and possible types of vertices are "Ppm" and "PPp"
	// if we take this m creates new vertex and there is new propags, we will make it
	// every mP propag split into two of type pM, MP
	// example: {A,B,(k),"mP"} -> {A,X,(0), "pM"}, {X,B,(k), "MP"}
	// then there will be new vertices "MM"
	// external point is always vertex {"e","e"}
public:
	// constructors
	Vertex();
	Vertex(std::string _vertName);
	Vertex(std::string _vertName, std::string _type);

	// getters
	std::string getVertName() const;
	std::string getVertType() const;

	// setters
	void setVertName(std::string _name);
	void setVertType(std::string _name);

	// print
	std::string info() const;
	void print() const;

	// overloading operators == and != to compare vertices
	bool operator==(const Vertex &rhs) const;
	bool operator!=(const Vertex &rhs) const;

	// overloading assignment operator = for vertices
	Vertex& operator=(const Vertex &rhs);
};

//===========================================================================

//===========================================================================
// class Propagator

class Propagator {
	/* class propagator holds oriented causal propagators. Namely its startpoint,
	 * endpoint, and momentum. Startpoint is always sooner in time (endpoint is later).
	 * TODO - implement type of propagator.
	 * */
protected:
	// attributes
	Vertex endVert; // later in time
	Vertex startVert; // sooner in time
	GiNaC::ex momentum; // always points from startVert to endVert
	std::string propType;
	// types for IP "mP", "pP", after taking "m" as extra vertex
	// types will be "pM", "MP", "pP"
	//example: {A,B,(k),"mP"} -> {A,X,(0), "pM"}, {X,B,(k), "MP"}

public:
	// constructors
	Propagator();
	Propagator(std::string _end, std::string _start, std::string mom);
	Propagator(std::string _end, std::string _start, std::string _mom,
			std::string _type);
	Propagator(Vertex _end, Vertex _start, GiNaC::ex _mom);
	Propagator(Vertex _end, Vertex _start, GiNaC::ex _mom, std::string _type);

	// getters
	Vertex getStartVert() const;
	Vertex getEndVert() const;
	GiNaC::ex getMomentum() const; //TODO - will this do what i want it to do??
	std::string getPropType() const;

	// setters
	void setStartVert(Vertex _start);
	void setEndVert(Vertex _end);
	void setMomentum(GiNaC::ex _mom);
	void setPropType(std::string _type);

	// printer
	std::string info() const;
	void print() const;

	// overloading operators == and != to compare propagators
	bool operator==(const Propagator &rhs) const;
	bool operator!=(const Propagator &rhs) const;

	// overloading assignment operator = for propagators
	Propagator& operator=(const Propagator &rhs);
};

//===========================================================================

//===========================================================================
// class Diagram
class Diagram {
	/* class diagram holds Feynman graphs. Every graph has a name, std::vector
	 * of propagators that correspond to edges of graph, and std::vector of vertices
	 * which correspond to nodes of the graph.
	 * */
protected:
	// attributes
	std::string name;
	std::vector<Propagator> extPropags;
	std::vector<Propagator> intPropags;
	std::vector<Vertex> vertices;
	Vertex extPoint; // is always just one type "E" -> otherwise you have a problem
	int numOfLoops;
	std::vector<GiNaC::symbol> intMomenta;
	std::vector<GiNaC::symbol> extMomenta;
	GiNaC::symbol tau;
	GiNaC::symbol u_0; // inside factorsFromVertsAndProps
	GiNaC::symbol D_0; // inside factorsFromVertsAndProps
	GiNaC::symbol d;
	int symmetryFactor;
	GiNaC::ex factorsFromVertsAndProps; // here is u_0 and D_0 - every Ppp gets +D_0 u_0, every Ppm gets -D_0 u_0, every Pm prop gets D_0
	// (later when time cuts are made there will be 1/D_0 for every cut - but that is later not in constructor here)
	//
	//
	// part necessary for calculating contributions propto external momentum p^2
	GiNaC::ex createdNumeratorByDerivative;
	std::vector<int> indicesOfPropsWhichContrToNumerator;
	GiNaC::symbol y1;
	GiNaC::symbol y2;

public:
	// constructors
	Diagram();
	Diagram(std::string _name, std::vector<Propagator> _propags,
			std::vector<Vertex> _verts, int _symmetryFactor);
	Diagram(std::string _name, std::vector<Propagator> _propags,
			std::vector<Vertex> _verts, int _numOfLoops, int _symmetryFactor);
	Diagram(std::string _name, std::vector<std::string> _propags,
			std::string _verts, int _symmetryFactor); // TODO - can i delete this?
	Diagram(std::string _nickel, int _symmetryFactor);
	Diagram(std::string _nickel, std::string _IPMomRouting,
			std::string _zeroExtMomenta, int _symmetryFactor);

	// getters
	std::string getName() const;
	std::vector<Propagator> getExtPropags() const;
	std::vector<Propagator> getIntPropags() const;
	std::vector<Vertex> getVertices() const;
	int getNumOfLoops() const;
	std::vector<GiNaC::symbol> getIntMomenta() const;
	std::vector<GiNaC::symbol> getExtMomenta() const;
	Vertex getExtPoint() const;
	Propagator getIntPropagAtIndex(int _index) const;
	Propagator getExtPropagAtIndex(int _index) const;
	Vertex getVertexAtIndex(int _index) const;
	GiNaC::symbol getTau() const;
	GiNaC::symbol getU_0() const;
	GiNaC::symbol getD_0() const;
	int getSymmetryFactor() const;
	GiNaC::ex getFactorsFromVertsAndProps() const;
	GiNaC::symbol getD() const;
	GiNaC::ex getCreatedNumeratorByDerivative() const;
	std::vector<int> getIndicesOfPropsWhichContrToNumerator() const;
	GiNaC::symbol getY1() const;
	GiNaC::symbol getY2() const;

	// setters
	void setName(std::string _name);
	void setExtPropags(std::vector<Propagator> _props);
	void setIntPropags(std::vector<Propagator> _props);
	void setVertices(std::vector<Vertex> _verts);
	void setNumOfLoops(int _num);
	void setIntMomenta(std::vector<GiNaC::symbol> _momenta);
	void setExtMomenta(std::vector<GiNaC::symbol> _momenta);
	void setExtPoint(Vertex _extPoint);
	void setIntPropagAtIndex(Propagator _propToSet, int _index);
	void setExtPropagAtIndex(Propagator _propToSet, int _index);
	void setVertexAtIndex(Vertex _vertToSet, int _index);
	void setTau(GiNaC::symbol _tau);
	void setU_0(GiNaC::symbol _u_0);
	void setD_0(GiNaC::symbol _D_0);
	void setSymmetryFactor(int _symmetryFactor);
	void setFactorsFromVertsAndProps(GiNaC::ex _factor);
	void setD(GiNaC::symbol _d);
	void setCreatedNumeratorByDerivative(GiNaC::ex _numerator);
	void setIndicesOfPropsWhichContrToNumerator(std::vector<int> _indices);
	void setY1(GiNaC::symbol _y1);
	void setY2(GiNaC::symbol _y2);

	// print
	void print() const;

	// overloading assignment operator =
	Diagram& operator=(const Diagram &rhs);

	// find possible time orderings
	// implementation can be found in TimeOrderingsAndTimeCuts.cpp
	std::vector<std::vector<Vertex>> findAllPossibleTimeOrderings(
			std::string _writeData);
	GiNaC::ex makeIntegrandAfterTimeCuts(
			std::vector<std::vector<Vertex>> _allTimeOrderings,
			std::string _writeData);

	// find number of loops
	int findNumOfLoopsPhi3(); // implementation in ElementsAndInput
};
//===========================================================================

//===========================================================================
// functions to be used in diagram constructors
std::vector<Propagator> findExtPropags(std::vector<Propagator> _propToSearch);
std::vector<Propagator> findIntPropags(std::vector<Propagator> _propToSearch);
std::vector<GiNaC::symbol> defSymbolMomenta();
std::vector<Vertex> findVertices(std::vector<Vertex> _vertsToSearch);
Vertex findExtPoint(std::vector<Vertex> _vertsToSearch);
GiNaC::ex findFactorFromVertsAndProps(Diagram _diag);
//===========================================================================

//===========================================================================
// functions to make program compatible with previous code
// TODO - DELETE
// following are used in constructor of Diagram
std::vector<Propagator> makeDiagramsPropagatorsFromStrings(
		std::vector<std::string> _input);
std::vector<Vertex> makeDiagramsVerticesFromString(std::string _verts);
//===========================================================================

//===========================================================================
// Code for constructor of Diagram from Nickel indices
// NOTE: should be performed on phi3 type of diagram, with
// vertices type Ppm, PPp.

std::vector<std::vector<std::string>> divideNickel(std::string _nickelInput);
std::vector<Vertex> findNodesFromNickel(std::string _nickelInput);
std::vector<Vertex> cleanRepeatedVertices(std::vector<Vertex> _vertices);
Vertex pickVertexByNameInsideVector(std::string _name,
		std::vector<Vertex> _input);
int findIndexOfVertexInsideVector(Vertex _vert, std::vector<Vertex> _input);
int findIndexOfPropInsideVector(Propagator _prop,
		std::vector<Propagator> _input);
std::vector<Propagator> findEdgesFromNickel(std::string _nickel,
		std::vector<Vertex> _verts);
std::string propTypeFromDivNickel(
		std::vector<std::vector<std::string>> _divNickel, int i, int j);
std::vector<std::vector<std::string>> divideDividedNickel(
		std::vector<std::vector<std::string>> _divNickel);
std::vector<Propagator> repairCausality(std::vector<Propagator> _propags);
void mySwap(std::string &_one, std::string &_two);
std::vector<Vertex> assignVertexTypes(std::vector<Vertex> _vertsWithoutType,
		std::vector<Propagator> _allEdges);
std::vector<Propagator> addVertTypesToProps(
		std::vector<Propagator> _propsWOVertType,
		std::vector<Vertex> _vertsWithType);

// asigning momenta
std::vector<GiNaC::symbol> addGeneralIntMomenta(Diagram &_diag);
std::vector<GiNaC::symbol> addGeneralExtMomenta(Diagram &_diag,
		std::string _zeroExtMomenta);
std::vector<Propagator> getPropsContainingVert(Vertex _vert,
		std::vector<Propagator> _propsToSearch);
std::vector<GiNaC::symbol> assignExtMom_IPMomRouting(Diagram &_diag,
		std::string _zeroExtMomenta);
std::vector<GiNaC::symbol> assignIntMom_IPMomRouting(Diagram &_diag);
bool doesMomConsHold(Diagram _diag);
std::vector<Propagator> assignMomByConservation(Diagram &_diag,
		std::vector<Propagator> &_usedProps);
//===========================================================================

//===========================================================================
// Code for dividing propags type mP
void dividePropsType_mP(Diagram &_diag);
//===========================================================================

//===========================================================================
// Code for getting part proportional to ext. frequency
std::vector<Diagram> addVertCorrespToFrequencyDeriv(Diagram _parentDiag);
Diagram setExtMomToZero(Diagram _inputDiag);
//===========================================================================

//===========================================================================
// Code for getting part proportional to tau
std::vector<Diagram> addVertCorrespToTauDeriv(Diagram _parentDiag);
//===========================================================================

//===========================================================================
// Code for getting part proportional to ext. momentum (p^2)
std::vector<Diagram> addVertsCorrespToPDeriv(Diagram &_parentDiag);

std::vector<Diagram> addVertType1(Diagram &_parentDiag);
std::vector<Diagram> addVertsType2(Diagram &_parentDiag);
std::vector<Diagram> addVertsType3_crossTerms(Diagram &_parentDiag);
//===========================================================================

#endif /* ELEMENTSANDINPUT_GINAC_HPP_ */
