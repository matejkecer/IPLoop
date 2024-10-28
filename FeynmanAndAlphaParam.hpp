/*
 * FeynmanAndAlphaParam.hpp
 *
 *  Created on: Aug 22, 2023
 *      Author: matej
 */

#ifndef FEYNMANANDALPHAPARAM_HPP_
#define FEYNMANANDALPHAPARAM_HPP_

#include "ElementsAndInput_ginac.hpp"
#include "TimeOrderingsAndTimeCuts_ginac.hpp"

class FeynmanParam {
	// The feynman param in general formula looks like
	// 1/(A1^{a1}...An^{an}) = gamma (alpha)/ [gamma(a1)...gamma(an)] int dx1...dxn delta(1-x1-x2-...-xn) *  x1^{a1-1} ... xn^{an-1}/Q^{alpha}
	// where alpha = a1+a2+...+an
	// after integration over loop momenta as and using Vasiliev formula->
	// (OverallFactor)*Int dx_1 ... dx_n *deltaFunction * numeratorFactor * (U)^{alpha-(l+1)*d/2}/(F)^{alpha-l*d/2}
	// where numeratorFactor = x1^{a1-1} ... xn^{an-1}, U = det(V), F = det(V)[C-A_i*A_j* V^{-1}_ij],
	// and V and A are from factoring out Q = x1*A1+...+xn*An = V_ij k_i k_j + 2 A_i k_i + C

	// see Binoth
protected:
	//std::vector<Vertex> timeOrdering;

	// diagram of which feynman parametrization is to be made
	Diagram diag;
	std::vector<Vertex> timeOrdering;
	// ginac symbols that will appear
	std::vector<GiNaC::symbol> loopMomenta;
	std::vector<GiNaC::symbol> extMomenta;
	GiNaC::symbol tau;
	GiNaC::symbol u_0, D_0; // symbols from diagram -> maybe its overkill to reinvent these here, and we could have used e.g.
	// diag.getU_0(), the same thing with D_0, tau, etc.
	std::vector<GiNaC::symbol> params; // feynman parameters
	std::vector<GiNaC::ex> powers_a; // a1, a2, ... in 1/(A1^a1 * A2^a2 *...) in general formula
	// powers_a are powers of denominators (propag. denoms) coresponding to feyn. parameters in attribute params
	// see numerator of general formula for feynman parametrization for more
	GiNaC::symbol d; //dimension of space (in IP it will be 6-\epsilon)

	// expressions and elements that will appear in feynman parametrization
	GiNaC::ex Q;
	GiNaC::ex alpha;
	GiNaC::ex C;
	GiNaC::matrix V;
	GiNaC::matrix A;
	// parts that will make up final resulting expression after
	// performing all the steps
	GiNaC::ex overallNumFactor; // overall numerical factor
	GiNaC::ex overallFactor2; // here gamma(2*eps) will be stored, also 4pi^{...eps} and factored tau^{...}
	//GiNaC::ex numerator;
	//GiNaC::ex denominator;

	// new ones // TODO - make proper notes
	GiNaC::ex numeratorFactor; // numerator in feynman param. -> paramters to respective powers ^(a_i-1)
	GiNaC::ex U; // det(V)
	GiNaC::ex F; // det(V) * (C - A_i*A_j*V^{-1}_ij)

	// for calculation of part propto p^2
	GiNaC::symbol y1;
	GiNaC::symbol y2;
public:
	// constructors
	FeynmanParam();
	FeynmanParam(Diagram &_diag, std::vector<Vertex> &_timeOrdering);

	//TODO - remaining getters and setters

	// getters
	Diagram getDiag() const;
	std::vector<Vertex> getTimeOrdering() const;

	std::vector<GiNaC::symbol> getLoopMomenta() const;
	std::vector<GiNaC::symbol> getExtMomenta() const;
	GiNaC::symbol getTau() const;
	GiNaC::symbol getU_0() const;
	GiNaC::symbol getD_0() const;
	std::vector<GiNaC::symbol> getParams() const;
	std::vector<GiNaC::ex> getPowers_a() const;
	GiNaC::symbol getD() const;

	GiNaC::ex getQ() const;
	GiNaC::ex getAlpha() const;
	GiNaC::ex getC() const;
	GiNaC::matrix getV() const;
	GiNaC::matrix getA() const;

	GiNaC::ex getOverallNumFactor() const;
	GiNaC::ex getOverallFactor2() const;
	GiNaC::ex getNumeratorFactor() const;
	GiNaC::ex getU() const;
	GiNaC::ex getF() const;

	GiNaC::symbol getY1() const;
	GiNaC::symbol getY2() const;

	// setters
	void setDiag(Diagram _diag);
	void setTimeOrdering(std::vector<Vertex> _timeOrdering);

	void setLoopMomenta(std::vector<GiNaC::symbol> _loopMom);
	void setExtMomenta(std::vector<GiNaC::symbol> _extMom);
	void setTau(GiNaC::symbol _tau);
	void setU_0(GiNaC::symbol _u_0);
	void setD_0(GiNaC::symbol _D_0);
	void setParams(std::vector<GiNaC::symbol> _setOfNus);
	void setPowers_a(std::vector<GiNaC::ex> _powers);
	void setD(GiNaC::symbol _d);

	void setQ(GiNaC::ex _Q);
	void setAlpha(GiNaC::ex _alpha);
	void setC(GiNaC::ex _C);
	void setV(GiNaC::matrix _V);
	void setA(GiNaC::matrix _A);

	void setOverallNumFactor(GiNaC::ex _overallNumFactor);
	void setOverallFactor2(GiNaC::ex _overallFactor2);
	void setNumeratorFactor(GiNaC::ex _numeratorFactor);
	void setU(GiNaC::ex _U);
	void setF(GiNaC::ex _F);

	void setY1(GiNaC::symbol _y1);
	void setY2(GiNaC::symbol _y2);
	// print
	void print() const;
};

//===========================================================================
// code to be used in constructor of Feynman parametrization
GiNaC::ex findQ(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau, std::vector<GiNaC::symbol> &_params, GiNaC::symbol &_y1, GiNaC::symbol &_y2);
std::vector<GiNaC::symbol> findParamSymbols(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau);
std::vector<GiNaC::ex> findPowers(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau);
GiNaC::ex findOverallNumFactor(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau);
GiNaC::ex findOverallNumFactor_type2(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau);
GiNaC::ex findOverallFactor2(FeynmanParam &_feyn);
GiNaC::ex findAlpha(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau);
GiNaC::matrix findV(GiNaC::ex &_Q, std::vector<GiNaC::symbol> &_k,
		std::vector<GiNaC::symbol> &_x);
GiNaC::matrix findA(GiNaC::ex &_Q, std::vector<GiNaC::symbol> &_k,
		std::vector<GiNaC::symbol> &_x);
GiNaC::ex findC(GiNaC::ex &_Q, std::vector<GiNaC::symbol> &_k,
		std::vector<GiNaC::symbol> &_x);
GiNaC::ex findNumeratorFactor(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau);
GiNaC::ex findU(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau, GiNaC::symbol &_y1, GiNaC::symbol &_y2);
GiNaC::ex findF(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau, GiNaC::symbol &_y1, GiNaC::symbol &_y2);
bool propDenomMatches(GiNaC::ex _denomA, GiNaC::ex _denomB,
		FeynmanParam &_feyn);
std::vector<std::vector<GiNaC::ex>> getPropDenomsForFeynmanParam_forBasicFormula(
		std::vector<std::vector<Propagator>> _cutsForGivenTimeOrdering,
		FeynmanParam &_feyn, GiNaC::symbol &_tau);
std::vector<std::vector<GiNaC::ex>> getPropDenomsForFeynmanParam_forGeneralFormula(
		std::vector<std::vector<Propagator>> _cutsForGivenTimeOrdering,
		FeynmanParam &_feyn, GiNaC::symbol &_tau);
//===========================================================================

//===========================================================================
// using quadratic formula (Vasiliev) - integration over loop momenta
// Feynman parametrization
void useQuadraticFormulaToIntegrateOverLoopMomenta(FeynmanParam &_feyn);
void factorTauOutsideIntoOverallFactor(FeynmanParam &_feyn);
//===========================================================================

//===========================================================================
// ALPHA PARAMETRIZATION - TODO

//===========================================================================

#endif /* FEYNMANANDALPHAPARAM_HPP_ */
