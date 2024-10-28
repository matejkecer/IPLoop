/*
 * FeynmanAndAlphaParam.cpp
 *
 *  Created on: Aug 22, 2023
 *      Author: matej
 */

#include "FeynmanAndAlphaParam.hpp"

// constructors
FeynmanParam::FeynmanParam() {

	//this->timeOrdering = { };
	this->diag = Diagram();
	this->timeOrdering = { };
	this->loopMomenta = { };
	this->extMomenta = { };
	this->tau = this->diag.getTau();
	this->u_0 = this->diag.getU_0();
	this->D_0 = this->diag.getD_0();
	this->Q = 1;
	this->params = { }; // set of feynman parameters x1, x2, ...
	this->powers_a = { };

	this->alpha = this->params.size();
	this->C = 1;
	this->V = GiNaC::matrix();
	this->overallNumFactor = 1;
	this->overallFactor2 = 1;
	this->A = GiNaC::matrix();
	this->d = GiNaC::symbol("d");

	this->numeratorFactor = 1;
	this->U = 1;
	this->F = 1;
	this->y1 = GiNaC::symbol("y1");
	this->y2 = GiNaC::symbol("y2");
}
FeynmanParam::FeynmanParam(Diagram &_diag, std::vector<Vertex> &_timeOrdering) {

	// NOTE: _diag should have IP Momentum routing
	// NOTE: constructs Feynman param by general formula -> does not integrate over loop momenta
	// the format of constructed param is written thorugh Q, overall factor and numeratorFactor

	// symbols we will further use

	this->diag = _diag;
	this->timeOrdering = _timeOrdering;
	this->loopMomenta = _diag.getIntMomenta();
	this->extMomenta = _diag.getExtMomenta();
	this->tau = this->diag.getTau();
	this->u_0 = this->diag.getU_0();
	this->D_0 = this->diag.getD_0();
	this->y1 = diag.getY1();
	this->y2 = diag.getY2();

	this->d = _diag.getD();
	this->params = findParamSymbols(*this, _timeOrdering, this->tau); // must be run first
	this->powers_a = findPowers(*this, _timeOrdering, this->tau);
	this->Q = findQ(*this, _timeOrdering, this->tau, this->params, this->y1,
			this->y2);
	this->alpha = findAlpha(*this, _timeOrdering, this->tau); // must be run prior to findOverallFactor()

	this->C = findC(this->Q, this->loopMomenta, this->params);
	this->V = findV(this->Q, this->loopMomenta, this->params);
	this->A = findA(this->Q, this->loopMomenta, this->params);

	this->overallNumFactor = findOverallNumFactor(*this, _timeOrdering,
					this->tau);
	this->overallFactor2 = findOverallFactor2(*this); // it gets correction in useQuadraticFormula() and factorTau()
	this->numeratorFactor = findNumeratorFactor(*this, _timeOrdering,
			this->tau);


	this->U = findU(*this, _timeOrdering, this->tau, this->y1, this->y2);
	this->F = findF(*this, _timeOrdering, this->tau, this->y1, this->y2);
}
// getters
Diagram FeynmanParam::getDiag() const {
	return this->diag;
}
std::vector<Vertex> FeynmanParam::getTimeOrdering() const {
	return this->timeOrdering;
}

std::vector<GiNaC::symbol> FeynmanParam::getLoopMomenta() const {
	return this->loopMomenta;
}
std::vector<GiNaC::symbol> FeynmanParam::getExtMomenta() const {
	return this->extMomenta;
}
GiNaC::symbol FeynmanParam::getTau() const {
	return this->tau;
}

GiNaC::symbol FeynmanParam::getU_0() const {
	return this->u_0;
}
GiNaC::symbol FeynmanParam::getD_0() const {
	return this->D_0;
}

std::vector<GiNaC::symbol> FeynmanParam::getParams() const {
	return this->params;
}
std::vector<GiNaC::ex> FeynmanParam::getPowers_a() const {
	return this->powers_a;
}

GiNaC::symbol FeynmanParam::getD() const {
	return this->d;
}

GiNaC::ex FeynmanParam::getQ() const {
	return this->Q;
}
GiNaC::ex FeynmanParam::getAlpha() const {
	return this->alpha;
}
GiNaC::ex FeynmanParam::getC() const {
	return this->C;
}
GiNaC::matrix FeynmanParam::getV() const {
	return this->V;
}
GiNaC::matrix FeynmanParam::getA() const {
	return this->A;
}

GiNaC::ex FeynmanParam::getOverallNumFactor() const {
	return this->overallNumFactor;
}
GiNaC::ex FeynmanParam::getOverallFactor2() const {
	return this->overallFactor2;
}
GiNaC::ex FeynmanParam::getNumeratorFactor() const {
	return this->numeratorFactor;
}
GiNaC::ex FeynmanParam::getU() const {
	return this->U;
}
GiNaC::ex FeynmanParam::getF() const {
	return this->F;
}
GiNaC::symbol FeynmanParam::getY1() const {
	return this->y1;
}
GiNaC::symbol FeynmanParam::getY2() const {
	return this->y2;
}

// setters
void FeynmanParam::setDiag(Diagram _diag) {
	this->diag = _diag;
	return;
}
void FeynmanParam::setTimeOrdering(std::vector<Vertex> _timeOrdering) {
	this->timeOrdering = _timeOrdering;
	return;
}

void FeynmanParam::setLoopMomenta(std::vector<GiNaC::symbol> _loopMom) {
	this->loopMomenta = _loopMom;
	return;
}
void FeynmanParam::setExtMomenta(std::vector<GiNaC::symbol> _extMom) {
	this->extMomenta = _extMom;
	return;
}
void FeynmanParam::setTau(GiNaC::symbol _tau) {
	this->tau = _tau;
	return;
}

void FeynmanParam::setU_0(GiNaC::symbol _u_0) {
	this->u_0 = _u_0;
	return;
}

void FeynmanParam::setD_0(GiNaC::symbol _D_0) {
	this->D_0 = _D_0;
	return;
}

void FeynmanParam::setParams(std::vector<GiNaC::symbol> _setOfNus) {
	this->params = _setOfNus;
	return;
}
void FeynmanParam::setPowers_a(std::vector<GiNaC::ex> _powers) {
	this->powers_a = _powers;
	return;
}

void FeynmanParam::setD(GiNaC::symbol _d) {
	this->d = _d;
	return;
}

void FeynmanParam::setQ(GiNaC::ex _Q) {
	this->Q = _Q;
	return;
}
void FeynmanParam::setAlpha(GiNaC::ex _alpha) {
	this->alpha = _alpha;
	return;
}
void FeynmanParam::setC(GiNaC::ex _C) {
	this->C = _C;
	return;
}
void FeynmanParam::setV(GiNaC::matrix _V) {
	this->V = _V;
	return;
}
void FeynmanParam::setA(GiNaC::matrix _A) {
	this->A = _A;
	return;
}

void FeynmanParam::setOverallNumFactor(GiNaC::ex _overallNumFactor) {
	this->overallNumFactor = _overallNumFactor;
	return;
}
void FeynmanParam::setOverallFactor2(GiNaC::ex _overallFactor2) {
	this->overallFactor2 = _overallFactor2;
	return;
}
void FeynmanParam::setNumeratorFactor(GiNaC::ex _numeratorFactor) {
	this->numeratorFactor = _numeratorFactor;
	return;
}
void FeynmanParam::setU(GiNaC::ex _U) {
	this->U = _U;
	return;
}
void FeynmanParam::setF(GiNaC::ex _F) {
	this->F = _F;
	return;
}

void FeynmanParam::setY1(GiNaC::symbol _y1) {
	this->y1 = _y1;
	return;
}
void FeynmanParam::setY2(GiNaC::symbol _y2) {
	this->y2 = _y2;
	return;
}

// print
void FeynmanParam::print() const {
	std::cout << "Feynman parametrization" << "\n";
	std::cout << "tau = " << this->tau << "\n";
	std::cout << "u_0 = " << this->u_0 << "\n";
	std::cout << "D_0 = " << this->D_0 << "\n";

	std::cout << "Params = {";
	for (int i = 0; i < this->params.size(); i++) {
		std::cout << this->params.at(i);
		if (i != this->params.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}";
	std::cout << "\n";

	std::cout << "Powers a_i  = {";
	for (int i = 0; i < this->powers_a.size(); i++) {
		std::cout << this->powers_a.at(i);
		if (i != this->powers_a.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}";
	std::cout << "\n";

	std::cout << "Loop momenta = {";
	for (int i = 0; i < this->loopMomenta.size(); i++) {
		std::cout << this->loopMomenta.at(i);
		if (i != this->loopMomenta.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}";
	std::cout << "\n";

	std::cout << "Ext. momenta = {";
	for (int i = 0; i < this->extMomenta.size(); i++) {
		std::cout << this->extMomenta.at(i);
		if (i != this->extMomenta.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}";
	std::cout << "\n";

	std::cout << "Q = " << this->Q << "\n";
	std::cout << "Alpha = " << this->alpha << "\n";
	std::cout << "C = " << this->C << "\n";
	std::cout << "V = " << this->V << "\n";
	std::cout << "A = " << this->A << "\n";
	std::cout << "---------------------------------" << "\n";
	std::cout << "Overall numerical factor = " << this->overallNumFactor
			<< "\n";
	std::cout << "Overall factor = " << this->overallFactor2 << "\n";
	std::cout << "Numerator factor = " << this->numeratorFactor << "\n";
	std::cout << "U = " << this->U << "\n";
	std::cout << "F = " << this->F << "\n";
	std::cout << "" << "\n";

}

//===========================================================================
// code to be used in constructor

std::vector<GiNaC::symbol> findParamSymbols(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau) {
	// Function creates ginac symbols for feynman parameters that will be used
	// in feynman parametrization.

	std::vector<GiNaC::symbol> result;
	std::vector<std::vector<GiNaC::ex>> info =
			getPropDenomsForFeynmanParam_forGeneralFormula(
					getCuts(_timeOrdering, _feyn.getDiag()), _feyn, _tau);
	std::vector<GiNaC::ex> denoms = info.at(0);

	// for generating new ginac symbol
	GiNaC::symbol current;
	std::stringstream name;

	// we need one new parameter for every propDenominator
	// e.g. 1/(A1*A2*...*An) will produce n-many feynman parameters
	// x1, ..., xn
	for (int i = 0; i < denoms.size(); i++) {
		name.str("");
		name << "x" << i;
		current = GiNaC::symbol(name.str());
		result.push_back(current);
	}

	return result;
}

std::vector<GiNaC::ex> findPowers(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau) {
	// Function finds and returns the vector of powers, corresponding to powers
	// of individual propDenominators that appear in integral
	// i.e. in general formula 1/(A1^a1 * A2^a2 * ... * An^an) = ...
	// the return of this function is vector {a1, a2, ..., an}

	std::vector<std::vector<GiNaC::ex>> info =
			getPropDenomsForFeynmanParam_forGeneralFormula(
					getCuts(_timeOrdering, _feyn.getDiag()), _feyn, _tau);

	return info.at(1);
}

GiNaC::ex findQ(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau, std::vector<GiNaC::symbol> &_params,
		GiNaC::symbol &_y1, GiNaC::symbol &_y2) {
	// Function determines the parameter Q in general formula for feynman
	// parametrization - that is Q = [x0*A0 + x1*A1 * ...]
	//
	// NOTE: can be used only after feynman parameters have been defined,
	// i.e. after findParamSymbols();

	GiNaC::ex result;
	std::vector<std::vector<GiNaC::ex>> info =
			getPropDenomsForFeynmanParam_forGeneralFormula(
					getCuts(_timeOrdering, _feyn.getDiag()), _feyn, _tau);
	std::vector<GiNaC::ex> denoms = info.at(0);

	for (int i = 0; i < denoms.size(); i++) {
		result += (_params.at(i) * denoms.at(i));
	}

	// TODO
	std::string controlName = "";
	for (int s = 0; s < 5; s++) {
		controlName += _feyn.getDiag().getName().at(s);
	}

	// if type2
	if (controlName.compare("type2") == 0) {
		std::vector<int> indices =
				_feyn.getDiag().getIndicesOfPropsWhichContrToNumerator();

		result +=
				2*_y1
						* _feyn.getDiag().getIntPropagAtIndex(indices.at(0)).getMomentum();
	}
	// if type3
	if (controlName.compare("type3") == 0) {
		std::vector<int> indices =
				_feyn.getDiag().getIndicesOfPropsWhichContrToNumerator();
		GiNaC::symbol y1 = _feyn.getY1();
		GiNaC::symbol y2 = _feyn.getY2();
		std::vector<GiNaC::symbol> symbols = { _y1, _y2 };

		for (int i = 0; i < 2; i++) {
			result +=
					2*symbols.at(i)
							* _feyn.getDiag().getIntPropagAtIndex(indices.at(i)).getMomentum();
		}
	}

	//std::cout<<result;
	return result;
}

GiNaC::ex findAlpha(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau) {
	// Function determines parameter alpha from general formula for feynman
	// parametrization (Q^alpha - that is the power of Q).

	GiNaC::ex result = 0;
	std::vector<std::vector<GiNaC::ex>> info =
			getPropDenomsForFeynmanParam_forGeneralFormula(
					getCuts(_timeOrdering, _feyn.getDiag()), _feyn, _tau);
	for (int i = 0; i < info.at(1).size(); i++) {
		result += info.at(1).at(i);
	}

	return result;
}

GiNaC::ex findOverallNumFactor(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau) {
	// Function finds overall factor that will be muliplying final expression
	// for feynman parametrization. It contains factor from denominators
	// ( like factor 1/2 in 1/[ (denom)*(2*(denom) ] ), and it also
	// contains gamma functions from general expression for feynman parametrization
	//
	// NOTE: needs to be performed after alpha is known - after findAlpha()

	GiNaC::ex result = 1;

	// add symmetry factor of diagram
	result *= _feyn.getDiag().getSymmetryFactor();
	//std::cout<<_feyn.getDiag().getSymmetryFactor()<<std::endl;

	std::vector<std::vector<GiNaC::ex>> info =
			getPropDenomsForFeynmanParam_forGeneralFormula(
					getCuts(_timeOrdering, _feyn.getDiag()), _feyn, _tau);

	// add numerical factor that is factored from used denominators
	result *= info.at(2).at(0);

	// add gamma functions resulting from feynman parametrization
	// NOTE: this actually holds for both basic and general formula gamma(a1+a2+...)/( gamma(a1)gamma(a2)...)
	// in basic formula donominator will just be 1 and upstairs there will be gamma(alpha)

	// first gamma in numerator
	result *= GiNaC::tgamma(_feyn.getAlpha());
	// then gammas in denominator
	for (int i = 0; i < info.at(1).size(); i++) {
		result *= 1 / (GiNaC::tgamma(info.at(1).at(i)));
	}

	result*=_feyn.getDiag().getCreatedNumeratorByDerivative();//.subs(_feyn.getD() == 6 - 2*_feyn.getDiag().get);
	return result;
}

GiNaC::ex findOverallFactor2(FeynmanParam &_feyn) {
	// Function adds factorFromVertsAndProps (D_0's and U_0's) from diagram
	// and updates it for D_0's that come out of cuts -> every cut gets 1/D_0
	GiNaC::ex result = 1;

	std::vector<std::vector<Propagator>> allCuts;
	allCuts = getCuts(_feyn.getTimeOrdering(), _feyn.getDiag());
	int power = allCuts.size();

	result *= _feyn.getDiag().getFactorsFromVertsAndProps();

	result *= GiNaC::pow(_feyn.getD_0(), -power);
	return result;
}

GiNaC::matrix findV(GiNaC::ex &_Q, std::vector<GiNaC::symbol> &_k,
		std::vector<GiNaC::symbol> &_x) {

	// Function determines the matrix V - of coeff. corresponding
	// to terms quadratic in loop momenta inside Q parameter from
	// general formula for feynman parametrization

	// initialize variables
	GiNaC::ex Q = _Q.normal(); // normal simplifies it
	GiNaC::matrix V(_k.size(), _k.size());

	// highest possible power that can be inside Q of any momentum k is square
	// coefficients of ki^2 make up diagonal elements of V, while coefs of ki*kj
	// make up off diagonal elements symmetrically V_ij = V_ji = (coef/2) *ki*kj

	// find diagonal elements
	for (int i = 0; i < _k.size(); i++) {
		V(i, i) = Q.diff(_k[i], 2) / 2;
	}

	// find off diagonal elements
	for (int i = 0; i < _k.size(); i++) {
		for (int j = i; j < _k.size(); j++) {

			if (i == j) {
				// diagonal elements are already set
				continue;
			}
			// the coefficient half for symmetrization
			V(i, j) = Q.diff(_k[i], 1).diff(_k[j], 1) / 2;
			V(j, i) = Q.diff(_k[i], 1).diff(_k[j], 1) / 2;
		}
	}

	return V;
}

GiNaC::matrix findA(GiNaC::ex &_Q, std::vector<GiNaC::symbol> &_k,
		std::vector<GiNaC::symbol> &_x) {
	// Function determines a column vector A (1xl matrix) from
	// general formula for feynman parametrization
	// A holds coefficients with terms of Q that are linear in
	// loop momenta

	// initialize variables
	GiNaC::ex Q = _Q.normal(); // normal simplifies it
	GiNaC::matrix A(1, _k.size());

	// find elements of A
	for (int i = 0; i < _k.size(); i++) {
		// differentiating w.r.t. loop momentum will get parts prop to it
		A(0, i) = Q.diff(_k[i], 1);

		// there can still be parts quadratic in loop momenta inside, so
		// set all remaining loop momenta to zero - what remains is part
		// linear in _k[i]
		for (int j = 0; j < _k.size(); j++) {
			A(0, i) = A(0, i).subs(_k[j] == 0);
		}

		// at the end you need to divide it by factor 2 to fit general formula
		A(0, i) = A(0, i) / 2;
	}

	return A;
}

GiNaC::ex findC(GiNaC::ex &_Q, std::vector<GiNaC::symbol> &_k,
		std::vector<GiNaC::symbol> &_x) {
	// Function determines parameter C from general formula for feynman
	// parametrization

	// C is basically part of Q that is not dependent on loop momenta
	// so initialize it to Q
	GiNaC::ex C = _Q;

	// since no loop momenta there - set all of them to zero and
	// what remains is a constant term
	for (int i = 0; i < _k.size(); i++) {
		C = C.subs(_k[i] == 0);
	}

	return factor(C);
}

GiNaC::ex findNumeratorFactor(FeynmanParam &_feyn,
		std::vector<Vertex> _timeOrdering, GiNaC::symbol &_tau) {
	GiNaC::ex result = 1;
	std::vector<std::vector<GiNaC::ex>> info =
			getPropDenomsForFeynmanParam_forGeneralFormula(
					getCuts(_timeOrdering, _feyn.getDiag()), _feyn, _tau);
	std::vector<GiNaC::ex> powers = info.at(1);

	for (int i = 0; i < powers.size(); i++) {
		result *= GiNaC::pow(_feyn.getParams().at(i), powers.at(i) - 1);
	}

	return result;
}

GiNaC::ex findU(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau, GiNaC::symbol &_y1, GiNaC::symbol &_y2) {
	GiNaC::ex result;
	result = GiNaC::determinant(_feyn.getV());
	return result;
}
GiNaC::ex findF(FeynmanParam &_feyn, std::vector<Vertex> _timeOrdering,
		GiNaC::symbol &_tau, GiNaC::symbol &_y1, GiNaC::symbol &_y2) {

	GiNaC::ex result = _feyn.getC();
	GiNaC::ex l = _feyn.getLoopMomenta().size();
	GiNaC::ex d = _feyn.getD();
	GiNaC::ex alpha = _feyn.getAlpha();
	GiNaC::matrix A = _feyn.getA();
	GiNaC::matrix V = _feyn.getV();
	GiNaC::matrix zero(1, _feyn.getLoopMomenta().size()); // for comparison
	GiNaC::matrix Vinv = V.inverse();
	//std::cout<<A(0,0);
	if (A == zero) {
		result *= _feyn.getU();
	} else {
		// if A is nonzero
		//std::cout<<A(0,0);
		for (int i = 0; i < _feyn.getLoopMomenta().size(); i++) {
			for (int s = 0; s < _feyn.getLoopMomenta().size(); s++) {
				result -= A(0, s) * Vinv(i, s) * A(0, i);
			}
		}
		result *= _feyn.getU();
	}
	return result;
}

bool propDenomMatches(GiNaC::ex _denomA, GiNaC::ex _denomB,
		FeynmanParam &_feyn) {

	// Function determines if _denomA = const * _denomB
	// NOTE: works for IP theory (kind of propagators that arise in it
	// in particular the expression for propags in IP, that is (mom^2 + tau) type )

	GiNaC::ex current;
	current = _denomA / _denomB;

	for (int i = 0; i < _feyn.getLoopMomenta().size(); i++) {
		if (GiNaC::has(current, _feyn.getLoopMomenta().at(i))) {
			return false;
		}
	}

	// if there are some external momenta
	if (_feyn.getExtMomenta().size() > 0) {
		for (int i = 0; i < _feyn.getExtMomenta().size(); i++) {
			if (GiNaC::has(current, _feyn.getExtMomenta().at(i))) {
				return false;
			}
		}
	}

	if (GiNaC::has(current, _feyn.getTau())) {
		return false;
	}

	// if I got to here, the expression current does not contain any
	// symbol that can be possibly present inside the expression for denominator
	// of whatever cut (or expression for propagator of the theory in general)
	// therefore it is mere number - a constant
	// therefore _denomA = const * _denomB
	// therefor denominators match
	return true;
}

std::vector<std::vector<GiNaC::ex>> getPropDenomsForFeynmanParam_forGeneralFormula(
		std::vector<std::vector<Propagator>> _cutsForGivenTimeOrdering,
		FeynmanParam &_feyn, GiNaC::symbol &_tau) {

	// return vector<vector<ginac::ex>> = { {denoms}, {their powers in integrand}, {overall num factor} }

	// Is used in Feynman parametrisation.
	// returns a vector of ginac expressions that correspond to denominators in integrand after time cuts.
	// that is propagator denominators (propDenoms) -> say if propag. is 1/A1 -> propDenom is A1
	//
	// e.g. integrand is  1/(A1*A2*A3) then return will hold {A1, A2, A3} - ginac expressions
	// e.g. integrand is 1/(A1*A1*A2) then return will hold {A1, A2} -> this is in order to use
	// the general formula for feynman param with powers i.e. 1/(A1^a1*...*An^an) = ...
	//
	// the function thus collects all the like propDenominators and extracts possible num factor
	// that will contribute to overall factor in front of integrals (there will be this numerical * symbols
	// that come out of verts, etc.)
	//
	// if there are two propDenominators which are merely multiples of each other,
	// e.g. (k^2+t) and (2k^2+2t) = 2*(k^2+t), then the one with lowest coefficient is put into return
	// vector i.e. (k^2+t) will get pushed back, the overall factor will be updated by 1/2
	//
	// however if there is only (2k^2+2t) = 2*(k^2+t) term in integrand and no (k^2+t), then the full
	// 2*(k^2+t) will put as propDenom into return, and overall factor will not be updated

	std::vector<GiNaC::ex> result = { };
	std::vector<int> powers = { };
	std::vector<GiNaC::ex> ginacPowers;

	GiNaC::ex overallNumFactor = 1;
	GiNaC::ex currentDenom;
	GiNaC::ex currentNumFactor = 1;
	int indexOfCurrentInResult = -1; //TODO - dont forget to reset it later in for cycle

	std::vector<GiNaC::ex> allDenoms;
	std::vector<GiNaC::ex> usedDenoms;

	// find propDenominators for all the cuts, without differentiating if two entries are
	// the same or multiples
	for (int i = 0; i < _cutsForGivenTimeOrdering.size(); i++) {
		allDenoms.push_back(
				findDenominatorForCut(_cutsForGivenTimeOrdering.at(i), _tau));
	}

// basiacally go through all cuts (propDenominators)
	for (int i = 0; i < allDenoms.size(); i++) {
		currentDenom = allDenoms.at(i);

		// if there is already something in final result vector (return)
		if (result.size() != 0) {

			// find out if currentDenom is already in there (in result)
			// and get its index in result
			for (int j = 0; j < result.size(); j++) {

				if (propDenomMatches(result.at(j), currentDenom, _feyn)) {
					indexOfCurrentInResult = j;
				}

			}
		}

		// if it already is in result
		if (indexOfCurrentInResult != -1) {
			// if yes

			// general working principle:
			// check if it isn't smaller multiple
			// i.e. now in result there can be 2*(k^2+t) and let currentDenom be (k^2+t)
			// currentNumFactor = (k^2+t)/(2*(k^2+t)) < 1
			// in that case you replace the spot in result by currentDenom
			// and update the overall num factor correspondingly

			currentNumFactor = currentDenom / result.at(indexOfCurrentInResult);

			// TODO - check if this gives the correct ginac comparison
			if (currentNumFactor < 1) {
				// if currentNumFactor < 1 - i don't have correct expression in the
				// result vector, I have to replace the expression in result by currentDenom,
				// multiply it by itself and update numerical factor
				// correspondingly, power of denom will be increased in vector of powers

				// how to update num factor -> e.g. currentDenom = (...), inResult = 2*(...)
				// say power of inResult is N. - i.e. integrand looks like 1/[ (2*(...))^N (...)]
				// =  1/[ (inRes)^N (current)] = (1/2)^N * 1/ (current)^(N+1) = (1/2)^N * 1/(newInResult)^(N+1)
				// thereforenumFactor should be updated by (1/2)^power i.e. by (currentDenom/inResult)^power
				// and only after this i replace inResult by currentDenom to form newInResult

				// update numerical factor
				overallNumFactor *= GiNaC::pow(currentNumFactor,
						powers.at(indexOfCurrentInResult));

				// then replace the expression in result by currentDenom
				result.at(indexOfCurrentInResult) = currentDenom;

				// add extra power corresponding to adding currentDenom
				powers.at(indexOfCurrentInResult)++;

				}
			else {
				// if currentNumFactor >= 1 - i have correct expression in the
				// result, I don't change anything in result but i update power in vector of powers
				// also i update numerical factor correspondingly.

				// how to update num factor -> e.g. currentDenom = 2*(...), inResult = 1*(...)
				// 1/(2*(...)) = ( 1/(inResult)^2 )*(1/2) -> numFactor should be updated
				// by 1/2 i.e. by inResult/currentDenom

				// TODO - check
				// this is always the case with -> you update with num factor by const. to
				// power 1, only when you don't have correct expression in result
				// then do you need to put some power to it
				overallNumFactor *= result.at(indexOfCurrentInResult)
						/ currentDenom;

				// multiply the given denominator by itself -> basically in result
				// you don't need to do anything - correct expression is there already
				// but its power will go up
				powers.at(indexOfCurrentInResult)++;}

			}
		else {
			// denom is not in result yet -> add it with power 1
			result.push_back(currentDenom);
			powers.push_back(1);
		}

		// reset aux variable
		indexOfCurrentInResult = -1;

	}

	GiNaC::ex aux;

	for (int i = 0; i < powers.size(); i++) {
		aux = powers.at(i);
		ginacPowers.push_back(aux);
	}

	return {result, ginacPowers, {overallNumFactor}};
}

std::vector<std::vector<GiNaC::ex>> getPropDenomsForFeynmanParam_forBasicFormula(
		std::vector<std::vector<Propagator>> _cutsForGivenTimeOrdering,
		FeynmanParam &_feyn, GiNaC::symbol &_tau) {

	// TODO - will this give the same result as general formula? will it even work?

	// return vector<vector<ginac::ex>> = { {denoms}, {their powers in integrand}, {overall num factor} }
	// actually - it being a basic formula the powers will all be one, and overall num factor too
	// the format of return is only such to be compatible with other functions such as find Q,C,V etc.

	// Is used in Feynman parametrisation.
	// returns a vector of ginac expressions that correspond to denominators in integrand after time cuts.
	// that is propagator denominators (propDenoms) -> say if propag. is 1/A1 -> propDenom is A1
	//
	// e.g. integrand is  1/(A1*A2*A3) then return will hold {A1, A2, A3} - ginac expressions
	// e.g. integrand is 1/(A1*A1*A2) then return will hold {A1, A1, A2} -> i.e. it doesn't
	// put common factors together and treats every denominator by itself as if it was different denom
	// even if the actuall expressions match,
	// the basic formula for feynman param is used i.e. 1/(A1*A2*...*An) = ...
	//
	// if there are two denominators which are merely multiples of each other,
	// e.g. (k^2+t) and (2k^2+2t) = 2*(k^2+t), then both will be treated as separate distinct
	// denoms, also the factor 2 from 2*(k^2+t) will not be factored out

	// TODO - test

	std::vector<GiNaC::ex> allDenoms = { };
	std::vector<GiNaC::ex> ginacPowers;

	GiNaC::ex overallNumFactor = 1;
	GiNaC::ex currentDenom;
	GiNaC::ex currentNumFactor = 1;
	int indexOfCurrentInResult = -1; //TODO - dont forget to reset it later in for cycle

	// find denominators for all the cuts, without differentiating if two entries are
	// the same or multiples
	for (int i = 0; i < _cutsForGivenTimeOrdering.size(); i++) {
		allDenoms.push_back(
				findDenominatorForCut(_cutsForGivenTimeOrdering.at(i), _tau));
	}

	ginacPowers = std::vector<GiNaC::ex>(allDenoms.size(), 1);

	return {allDenoms, ginacPowers, {overallNumFactor}};
}
//===========================================================================

//===========================================================================
// using quadratic formula (Vasiliev) - integration over loop momenta

void useQuadraticFormulaToIntegrateOverLoopMomenta(FeynmanParam &_feyn) {

	// quadratic formula says
	// 1/(2 \pi)^{l*d} int d^d k1 ... int d^d kl 1/(V_{is} k_i k_s + 2 A_i k_i + C)^Alpha =
	// 1/(4 \pi)^{l*d/2} [gamma(Alpha - l*d/2)/gamma(Alpha)] * [(det V)^{-d/2}/(C - V^{-1}_{is} A_i A_s)^{Alpha-l*d/2}]

	GiNaC::ex l = _feyn.getLoopMomenta().size();
	GiNaC::ex d = _feyn.getD();
	GiNaC::ex alpha = _feyn.getAlpha();
	GiNaC::matrix A = _feyn.getA();
	GiNaC::matrix V = _feyn.getV();


	// a) update overallFactor by (4 \pi) thing and gamma functions
	GiNaC::ex newOverallNumFactor = _feyn.getOverallNumFactor();
	GiNaC::ex newOverallFactor2 = _feyn.getOverallFactor2();

	newOverallFactor2 *= 1 / GiNaC::pow((4 * GiNaC::Pi), l * d / 2);
	//newOverallFactor2 *= GiNaC::tgamma(alpha - l * d / 2) / GiNaC::tgamma(alpha);
	newOverallFactor2 *= GiNaC::tgamma(alpha - l * d / 2);
	newOverallNumFactor /= GiNaC::tgamma(alpha);

	_feyn.setOverallFactor2(newOverallFactor2);
	_feyn.setOverallNumFactor(newOverallNumFactor);

	// b) find determinant of V
	GiNaC::ex detV = V.determinant();
	_feyn.setU(detV);
	//GiNaC::ex newNumerator = _feyn.getNumerator();
	//newNumerator *= GiNaC::pow(detV, -d / 2);
	//_feyn.setNumerator(newNumerator);

	// c) find inverse of V - only if A is nonzero vector
	// now A is matrix of size 1 x loopMomenta.size() - this is used in comparison, after that
	// set FeynmanParam attribute denominator - denominator on rhs i.e. (C - V^{-1}_{is} A_i A_s)^{Alpha-l*d/2}
	// if A is zero matrix - only C will contribute

	GiNaC::matrix zero(1, _feyn.getLoopMomenta().size());	// for comparison
	GiNaC::matrix Vinv = V.inverse();

	GiNaC::ex newDenom = _feyn.getC();

	if (A == zero) {
		_feyn.setF(newDenom * detV);
		//_feyn.setDenominator(GiNaC::pow(_feyn.getC(), alpha - l * d / 2));
	} else {
		// if A is nonzero

		for (int i = 0; i < _feyn.getLoopMomenta().size(); i++) {
			for (int s = 0; s < _feyn.getLoopMomenta().size(); s++) {
				newDenom -= Vinv(i, s) * A(0, i) * A(0, s);
			}
		}

		_feyn.setF(newDenom * detV);

		//newDenom = GiNaC::pow(newDenom, alpha - l * d / 2);
		//_feyn.setDenominator(newDenom);
	}

	return;
}

void factorTauOutsideIntoOverallFactor(FeynmanParam &_feyn) {
	// NOTE: works for IP for diagrams proportional to external frequency where
	// external momenta are set to 0. In that case A = 0 and after momentum integration
	// by help of quadratic formula, tau pops out of denominator as tau^{-(Alpha-l*d/2)}.
	// What remains in the denominator are some linear combination of feynman parameters
	// that is (C/tau)^{alpha - l*d/2}.

	// Besides that, after this step only things that remain are numerator where there is the determinant
	// and delta function, and integrals over feynman parameters

	GiNaC::ex l = _feyn.getLoopMomenta().size();
	GiNaC::ex d = _feyn.getD();
	GiNaC::ex alpha = _feyn.getAlpha();
	GiNaC::symbol tau = _feyn.getTau();
	GiNaC::matrix A = _feyn.getA();
	GiNaC::matrix zero(1, _feyn.getLoopMomenta().size());	// for comparison

	// a) message - in case diag. prop to ext momentum is given as input
	if ((A != zero)) {
		std::cout
				<< "Diagram is proportional to external momenta - vector A is non-zero."
						"Tau cannot be straightforwardly factored out."
						"-factorTauOutsideIntoOverallFactor()" << std::endl;
		return;
	}
	// b) update overall Factor
	GiNaC::ex newOverallFactor2 = _feyn.getOverallFactor2();
	newOverallFactor2 *= GiNaC::pow(tau, -alpha + l * d / 2);
	_feyn.setOverallFactor2(newOverallFactor2);

	// c) update FeynmanParam attribute denominator

	//GiNaC::ex newDenominator = _feyn.getDenominator();
	//newDenominator = newDenominator.subs(tau == 1);
	//_feyn.setDenominator(newDenominator);

	// d) update F
	GiNaC::ex newF = _feyn.getF();
	newF = newF.subs(tau == 1);
	_feyn.setF(newF);

	return;
}

//===========================================================================

