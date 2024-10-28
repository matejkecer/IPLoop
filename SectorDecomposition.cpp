/*
 * SectorDecomposition.cpp
 *
 *  Created on: 25. 8. 2023
 *      Author: matej
 */

#include "SectorDecomposition.hpp"

//===========================================================================
// Class Sector
std::vector<GiNaC::symbol> Sector::getExtMomenta() const {
	return this->extMomenta;
}
GiNaC::symbol Sector::getTau() const {
	return this->tau;
}
GiNaC::symbol Sector::getU0() const {
	return this->u0;
}
GiNaC::symbol Sector::getD0() const {
	return this->D0;
}
GiNaC::symbol Sector::getY1() const {
	return this->y1;
}
GiNaC::symbol Sector::getY2() const {
	return this->y2;
}
std::vector<GiNaC::symbol> Sector::getSectorVars() const {
	return this->sectorVars;
}
std::vector<GiNaC::ex> Sector::getPowers_a_bar() const {
	return this->powers_a_bar;
}
GiNaC::symbol Sector::getD() const {
	return this->d;
}
GiNaC::symbol Sector::getEps() const {
	return this->eps;
}

GiNaC::ex Sector::getAlpha() const {
	return this->alpha;
}
GiNaC::ex Sector::getL() const {
	return this->l;
}
GiNaC::ex Sector::getC_bar() const {
	return this->C_bar;
}
GiNaC::matrix Sector::getV_bar() const {
	return this->V_bar;
}
GiNaC::matrix Sector::getA_bar() const {
	return this->A_bar;
}

GiNaC::ex Sector::getOverallNumFactor_bar() const {
	return this->overallNumFactor_bar;
}
GiNaC::ex Sector::getOverallFactor2_bar() const {
	return this->overallFactor2_bar;
}
GiNaC::ex Sector::getNumeratorFactor_bar() const {
	return this->numeratorFactor_bar;
}
GiNaC::ex Sector::getU_bar() const {
	return this->U_bar;
}
GiNaC::ex Sector::getF_bar() const {
	return this->F_bar;
}

FeynmanParam Sector::getParentFeyn() const {
	return this->parentFeyn;
}

// setters
void Sector::setExtMomenta(std::vector<GiNaC::symbol> _extMomenta) {
	this->extMomenta = _extMomenta;
	return;
}
void Sector::setTau(GiNaC::symbol _tau) {
	this->tau = _tau;
	return;
}
void Sector::setU0(GiNaC::symbol _u0) {
	this->u0 = _u0;
	return;
}
void Sector::setD0(GiNaC::symbol _D0) {
	this->D0 = _D0;
	return;
}
void Sector::setY1(GiNaC::symbol _y1) {
	this->y1 = _y1;
	return;
}
void Sector::setY2(GiNaC::symbol _y2) {
	this->y2 = _y2;
	return;
}
void Sector::setSectorVars(std::vector<GiNaC::symbol> _sectorVars) {
	this->sectorVars = _sectorVars;
	return;
}
void Sector::setPowers_a_bar(std::vector<GiNaC::ex> _powers_a_bar) {
	this->powers_a_bar = _powers_a_bar;
	return;
}
void Sector::setD(GiNaC::symbol _d) {
	this->d = _d;
	return;
}
void Sector::setEps(GiNaC::symbol _eps) {
	this->eps = _eps;
	return;
}

void Sector::setAlpha(GiNaC::ex _alpha) {
	this->alpha = _alpha;
	return;
}
void Sector::setL(GiNaC::ex _l) {
	this->l = _l;
	return;
}
void Sector::setC_bar(GiNaC::ex _C_bar) {
	this->C_bar = _C_bar;
	return;
}
void Sector::setV_bar(GiNaC::matrix _V_bar) {
	this->V_bar = _V_bar;
	return;
}
void Sector::setA_bar(GiNaC::matrix _A_bar) {
	this->A_bar = _A_bar;
	return;
}

void Sector::setOverallNumFactor_bar(GiNaC::ex _overallNumFactor_bar) {
	this->overallNumFactor_bar = _overallNumFactor_bar;
	return;
}
void Sector::setOverallFactor2_bar(GiNaC::ex _overallFactor2_bar) {
	this->overallFactor2_bar = _overallFactor2_bar;
	return;
}
void Sector::setNumeratorFactor_bar(GiNaC::ex _numerator_bar) {
	this->numeratorFactor_bar = _numerator_bar;
	return;
}
void Sector::setU_bar(GiNaC::ex _U_bar) {
	this->U_bar = _U_bar;
	return;
}
void Sector::setF_bar(GiNaC::ex _F_bar) {
	this->F_bar = _F_bar;
	return;
}

void Sector::setParentFeyn(FeynmanParam _parentFeyn) {
	this->parentFeyn = _parentFeyn;
	return;
}

void Sector::print() const {

	//std::cout << "Sector" << "\n";

	std::cout << "Sector variables = {";
	for (int i = 0; i < this->sectorVars.size(); i++) {
		std::cout << this->sectorVars.at(i);
		if (i != this->sectorVars.size() - 1) {
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

	std::cout << "tau = " << this->tau << "\n";
	std::cout << "Alpha = " << this->alpha << "\n";
	std::cout << "C_bar = " << this->C_bar << "\n";
	std::cout << "V_bar = " << this->V_bar << "\n";
	std::cout << "A_bar = " << this->A_bar << "\n";
	std::cout << "---------------------------------" << "\n";
	std::cout << "OverallNumFactor_bar = " << this->overallNumFactor_bar
			<< "\n";
	std::cout << "OverallFactor2_bar = " << this->overallFactor2_bar << "\n";
	std::cout << "NumeratorFactor_bar = " << this->numeratorFactor_bar << "\n";
	std::cout << "U_bar = " << this->U_bar << "\n";
	std::cout << "F_bar = " << this->F_bar << "\n";

	std::cout << "" << "\n";
}

//===========================================================================

//===========================================================================
// Class PrimarySector

PrimarySector::PrimarySector() {
	this->parentFeyn = FeynmanParam();
	this->paramIntegratedOut = GiNaC::symbol();

	this->extMomenta = { };
	this->tau = GiNaC::symbol();
	this->sectorVars = { };
	this->powers_a_bar = { };

	this->d = GiNaC::symbol();
	this->eps = GiNaC::symbol();
	this->alpha = 0;
	this->l = 0;
	this->u0 = GiNaC::symbol();
	this->D0 = GiNaC::symbol();

	this->C_bar = 0;
	this->V_bar = GiNaC::matrix();
	this->A_bar = GiNaC::matrix();

	// parts that will make up final resulting expression after
	// performing all the steps
	this->overallNumFactor_bar = 0;
	this->overallFactor2_bar = 0;
	this->numeratorFactor_bar = 0;
	this->U_bar = 0;
	this->F_bar = 0;
	this->y1 = GiNaC::symbol();
	this->y2 = GiNaC::symbol();
}
PrimarySector::PrimarySector(FeynmanParam &_feyn,
		GiNaC::symbol &_paramIntegratedOut_xl) {
	this->parentFeyn = _feyn;
	this->paramIntegratedOut = _paramIntegratedOut_xl;

	this->extMomenta = _feyn.getExtMomenta();
	this->tau = _feyn.getTau();
	this->d = _feyn.getD();
	this->eps = GiNaC::symbol("eps");
	this->alpha = _feyn.getAlpha();
	this->u0 = _feyn.getU_0();
	this->D0 = _feyn.getD_0();
	this->l = _feyn.getLoopMomenta().size();

	this->sectorVars = findSectorVars(this->parentFeyn,
			this->paramIntegratedOut);
	this->powers_a_bar = findPowers_a_bar(this->parentFeyn,
			this->paramIntegratedOut);
	this->C_bar = findC_bar(this->parentFeyn, this->paramIntegratedOut, *this);
	this->V_bar = findV_bar(this->parentFeyn, this->paramIntegratedOut, *this);
	this->A_bar = findA_bar(this->parentFeyn, this->paramIntegratedOut, *this);

	// parts that will make up final resulting expression after
	// performing all the steps
	this->overallNumFactor_bar = findOverallNumFactor_bar(this->parentFeyn);//
	//.subs(this->d==6-2*this->eps);
	this->overallFactor2_bar = findOverallFactor2_bar(this->parentFeyn);
	this->numeratorFactor_bar = findNumeratorFactor_bar(this->parentFeyn,
			this->paramIntegratedOut, *this);
	this->U_bar = findU_bar(this->parentFeyn, this->paramIntegratedOut, *this);
	this->F_bar = findF_bar(this->parentFeyn, this->paramIntegratedOut, *this);

	this->y1 = _feyn.getY1();
	this->y2 = _feyn.getY2();
}

// getters
//FeynmanParam PrimarySector::getParenFeyn() const {
//	return this->parentFeyn;
//}
GiNaC::symbol PrimarySector::getParamIntegratedOut() const {
	return this->paramIntegratedOut;
}

// setters
//void PrimarySector::setParenFeyn(FeynmanParam _parentFeyn) {
//	this->parentFeyn = _parentFeyn;
//	return;
//}
void PrimarySector::setParamIntegratedOut(GiNaC::symbol _param) {
	this->paramIntegratedOut = _param;
	return;
}

void PrimarySector::print() const {

	std::cout << "-----------------------------" << "\n";
	std::cout << "Parent ";
	this->parentFeyn.print();
	std::cout << "-----------------------------" << "\n";

	std::cout << "Primary sector" << "\n";
	std::cout << "Parameter integrated out = {" << this->paramIntegratedOut
			<< "} \n";

	std::cout << "Sector variables = {";
	for (int i = 0; i < this->sectorVars.size(); i++) {
		std::cout << this->sectorVars.at(i);
		if (i != this->sectorVars.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}";
	std::cout << "\n";

	std::cout << "Powers a_i_bar  = {";
	for (int i = 0; i < this->powers_a_bar.size(); i++) {
		std::cout << this->powers_a_bar.at(i);
		if (i != this->powers_a_bar.size() - 1) {
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

	std::cout << "tau = " << this->tau << "\n";
	std::cout << "Alpha = " << this->alpha << "\n";
	std::cout << "C_bar = " << this->C_bar << "\n";
	std::cout << "V_bar = " << this->V_bar << "\n";
	std::cout << "A_bar = " << this->A_bar << "\n";
	std::cout << "---------------------------------" << "\n";
	std::cout << "OverallNumFactor_bar = " << this->overallNumFactor_bar
			<< "\n";
	std::cout << "OverallFactor2_bar = " << this->overallFactor2_bar << "\n";
	std::cout << "NumeratorFactor_bar = " << this->numeratorFactor_bar << "\n";
	std::cout << "U_bar = " << this->U_bar << "\n";
	std::cout << "F_bar = " << this->F_bar << "\n";

	std::cout << "" << "\n";

}

//===========================================================================
// Code used in constructor
std::vector<GiNaC::symbol> PrimarySector::findSectorVars(FeynmanParam &_feyn,
		GiNaC::symbol &_xl) {
	// Function creates and returns ginac symbols for sector variables in primary sector
	// Argument _feyn - feynman parametrization from which we are making primary sectors
	// Argument _xl - feynman parameter that is integrated out in the process of making primary
	// 			sectors

	std::vector<GiNaC::symbol> sectorVars;
	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	bool laterThanL = false;

	// for generating new ginac symbol
	GiNaC::symbol current;
	std::stringstream name;

	for (int i = 0; i < parentFeynParams.size(); i++) {

		// we need to leave out the parameter that is integrated out
		if (parentFeynParams.at(i) == _xl) {
			laterThanL = true;
			continue;
		}

		// if i < l
		if (!laterThanL) {
			// otherwise create new symbol corresponding to sector variable
			name.str("");
			name << "t" << i;
			current = GiNaC::symbol(name.str());
			sectorVars.push_back(current);
		} else {
			// if i > l
			// otherwise create new symbol corresponding to sector variable
			name.str("");
			name << "t" << i - 1;
			current = GiNaC::symbol(name.str());
			sectorVars.push_back(current);
		}

	}

	return sectorVars;
}

std::vector<GiNaC::ex> PrimarySector::findPowers_a_bar(FeynmanParam &_feyn,
		GiNaC::symbol &_xl) {
	std::vector<GiNaC::ex> powers_bar;
	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::ex> parentPowers = _feyn.getPowers_a();

	for (int i = 0; i < parentFeynParams.size(); i++) {

		// we need to leave out the parameter that is integrated out
		if (parentFeynParams.at(i) == _xl) {
			continue;
		}

		// this readily makes powers a_i_bar corresponding to new sector vars
		// it is consistent with the substitution made in generation of primary
		// sectors
		powers_bar.push_back(parentPowers.at(i));
	}

	return powers_bar;
}

GiNaC::ex PrimarySector::findC_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
		PrimarySector &_sect) {
	// Function rewrites C of original Feynman parametrization into new sector
	// variables. Basically it has three steps:
	// 1.) if i < l then x_i -> t_i
	// 2.) subs x_l == 1
	// 3.) if i > l then x_i -> t_(i-1)

	GiNaC::ex result = _feyn.getC();

	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::symbol> sectorVars = _sect.getSectorVars();
	bool laterThanL = false;

	for (int i = 0; i < parentFeynParams.size(); i++) {

		// we need to leave out the parameter that is integrated out
		if (parentFeynParams.at(i) == _xl) {
			// since we factored it out 1 should be substituted for it
			result = result.subs(parentFeynParams.at(i) == 1);
			laterThanL = true;
			continue;
		}

		// if i < l
		if (!laterThanL) {
			result = result.subs(parentFeynParams.at(i) == sectorVars.at(i));
		} else {
			// if i > l
			result = result.subs(
					parentFeynParams.at(i) == sectorVars.at(i - 1));
		}
	}

	return result;
}

GiNaC::matrix PrimarySector::findV_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
		PrimarySector &_sect) {
	// Function rewrites V of original Feynman parametrization into new sector
	// variables. Basically it has three steps:
	// 1.) if i < l then x_i -> t_i
	// 2.) subs x_l == 1
	// 3.) if i > l then x_i -> t_(i-1)

	GiNaC::matrix result = _feyn.getV(); //

	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::symbol> sectorVars = _sect.getSectorVars();
	bool laterThanL = false;

	// for every parameter in original Feynman parametrization
	for (int k = 0; k < parentFeynParams.size(); k++) {

		// if we have xl
		if (parentFeynParams.at(k) == _xl) {

			// go throgh all elements of the matrix
			for (int i = 0; i < result.rows(); i++) {
				for (int j = 0; j < result.cols(); j++) {

					// since we factored it out 1 should be substituted for it
					result(i, j) = result(i, j).subs(
							parentFeynParams.at(k) == 1);

				}
			}

			laterThanL = true;
			continue;
		}

		if (!laterThanL) {
			// if i < l
			// go throgh all elements of the matrix
			for (int i = 0; i < result.rows(); i++) {
				for (int j = 0; j < result.cols(); j++) {
					// substitute in new sector variables
					result(i, j) = result(i, j).subs(
							parentFeynParams.at(k) == sectorVars.at(k));
				}
			}

		} else {
			// if i > l
			// go throgh all elements of the matrix
			for (int i = 0; i < result.rows(); i++) {
				for (int j = 0; j < result.cols(); j++) {
					// substitute in new sector variables
					result(i, j) = result(i, j).subs(
							parentFeynParams.at(k) == sectorVars.at(k - 1));
				}
			}
		}

	}
	return result;
}

GiNaC::matrix PrimarySector::findA_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
		PrimarySector &_sect) {
	// Function rewrites V of original Feynman parametrization into new sector
	// variables. Basically it has three steps:
	// 1.) if i < l then x_i -> t_i
	// 2.) subs x_l == 1
	// 3.) if i > l then x_i -> t_(i-1)

	GiNaC::matrix result = _feyn.getA(); //

	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::symbol> sectorVars = _sect.getSectorVars();
	bool laterThanL = false;

	// for every parameter in original Feynman parametrization
	for (int k = 0; k < parentFeynParams.size(); k++) {
		// if we have xl
		if (parentFeynParams.at(k) == _xl) {

			// go throgh all elements of the matrix
			for (int i = 0; i < result.cols(); i++) {
				// since we factored it out 1 should be substituted for it
				result(0, i) = result(0, i).subs(parentFeynParams.at(k) == 1);
			}

			laterThanL = true;
			continue;
		}

		if (!laterThanL) {
			// if i < l
			// go throgh all elements of the matrix
			for (int i = 0; i < result.cols(); i++) {
				// substitute in new sector variables
				result(0, i) = result(0, i).subs(
						parentFeynParams.at(k) == sectorVars.at(k));

			}

		} else {
			// if i > l
			// go throgh all elements of the matrix
			for (int i = 0; i < result.cols(); i++) {
				// substitute in new sector variables
				result(0, i) = result(0, i).subs(
						parentFeynParams.at(k) == sectorVars.at(k - 1));

			}
		}

	}
	return result;
}

GiNaC::ex PrimarySector::findOverallNumFactor_bar(FeynmanParam &_feyn) {
	// overall factor does not change
	return _feyn.getOverallNumFactor();
}
GiNaC::ex PrimarySector::findOverallFactor2_bar(FeynmanParam &_feyn) {
	// overall factor does not change
	return _feyn.getOverallFactor2();
}
GiNaC::ex PrimarySector::findNumeratorFactor_bar(FeynmanParam &_feyn,
		GiNaC::symbol &_xl, PrimarySector &_sect) {
	// one way to do it is to reproduce the steps as we have done
	// for FeynmanParam but now starting with V_bar, A_bar, etc.,
	// easier way is to just substitute in sector variables and
	// sector powers a_bar into final expression for numerator from
	// FeynmanParam

	GiNaC::ex result = _feyn.getNumeratorFactor();
	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::symbol> sectorVars = _sect.getSectorVars();
	bool laterThanL = false;

	for (int i = 0; i < parentFeynParams.size(); i++) {

		// we need to leave out the parameter that is integrated out
		if (parentFeynParams.at(i) == _xl) {
			// since we factored it out 1 should be substituted for it
			result = result.subs(parentFeynParams.at(i) == 1);
			laterThanL = true;
			continue;
		}

		// if i < l
		if (!laterThanL) {
			result = result.subs(parentFeynParams.at(i) == sectorVars.at(i));
		} else {
			// if i > l
			result = result.subs(
					parentFeynParams.at(i) == sectorVars.at(i - 1));
		}
	}
	return result;
}
GiNaC::ex PrimarySector::findU_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
		PrimarySector &_sect) {
	// one way to do it is to reproduce the steps as we have done
	// for FeynmanParam but now starting with V_bar, A_bar, etc.,
	// easier way is to just substitute in sector variables and
	// sector powers a_bar into final expression for numerator from
	// FeynmanParam

	GiNaC::ex result = _feyn.getU();
	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::symbol> sectorVars = _sect.getSectorVars();
	bool laterThanL = false;

	for (int i = 0; i < parentFeynParams.size(); i++) {

		// we need to leave out the parameter that is integrated out
		if (parentFeynParams.at(i) == _xl) {
			// since we factored it out 1 should be substituted for it
			result = result.subs(parentFeynParams.at(i) == 1);
			laterThanL = true;
			continue;
		}

		// if i < l
		if (!laterThanL) {
			result = result.subs(parentFeynParams.at(i) == sectorVars.at(i));
		} else {
			// if i > l
			result = result.subs(
					parentFeynParams.at(i) == sectorVars.at(i - 1));
		}
	}
	result = result.normal();
	return result;
}
GiNaC::ex PrimarySector::findF_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
		PrimarySector &_sect) {
	// one way to do it is to reproduce the steps as we have done
	// for FeynmanParam but now starting with V_bar, A_bar, etc.,
	// easier way is to just substitute in sector variables and
	// sector powers a_bar into final expression for numerator from
	// FeynmanParam

	GiNaC::ex result = _feyn.getF();
	std::vector<GiNaC::symbol> parentFeynParams = _feyn.getParams();
	std::vector<GiNaC::symbol> sectorVars = _sect.getSectorVars();
	bool laterThanL = false;

	for (int i = 0; i < parentFeynParams.size(); i++) {

		// we need to leave out the parameter that is integrated out
		if (parentFeynParams.at(i) == _xl) {
			// since we factored it out 1 should be substituted for it
			result = result.subs(parentFeynParams.at(i) == 1);
			laterThanL = true;
			continue;
		}

		// if i < l
		if (!laterThanL) {
			result = result.subs(parentFeynParams.at(i) == sectorVars.at(i));
		} else {
			// if i > l
			result = result.subs(
					parentFeynParams.at(i) == sectorVars.at(i - 1));
		}
	}
	result = result.normal();
	return result;
}
//===========================================================================

//===========================================================================
// Class Subsector
SubSector::SubSector() {
	this->parent = Sector();
	this->largestParam = GiNaC::symbol();

	this->extMomenta = { };
	this->tau = GiNaC::symbol();
	this->sectorVars = { };
	this->powers_a_bar = { };

	this->d = GiNaC::symbol();
	this->u0 = GiNaC::symbol();
	this->D0 = GiNaC::symbol();
	this->eps = GiNaC::symbol();
	this->alpha = 0;
	this->l = 0;

	this->C_bar = 0;
	this->V_bar = GiNaC::matrix();
	this->A_bar = GiNaC::matrix();

	// parts that will make up final resulting expression after
	// performing all the steps
	this->overallNumFactor_bar = 0;
	this->overallFactor2_bar = 0;
	this->numeratorFactor_bar = 0;
	this->U_bar = 0;
	this->F_bar = 0;

	this->y1 = GiNaC::symbol();
	this->y2 = GiNaC::symbol();
}
SubSector::SubSector(Sector &_parent, std::vector<GiNaC::symbol> _params,
		GiNaC::symbol _largestParam_tl) {

	this->parent = _parent;
	this->largestParam = _largestParam_tl; // was used for subsector var substitution

	this->extMomenta = _parent.getExtMomenta();
	this->tau = _parent.getTau();
	this->d = _parent.getD();
	this->u0 = _parent.getU0();
	this->D0 = _parent.getD0();
	this->eps = GiNaC::symbol("eps");
	this->alpha = _parent.getAlpha();
	this->l = _parent.getL();

	this->sectorVars = findSectorVars(this->parent, this->largestParam);
	this->powers_a_bar = findPowers_a_bar(this->parent, this->largestParam);
	this->C_bar = findC_bar(this->parent, _params, this->largestParam, *this);
	this->V_bar = findV_bar(this->parent, _params, this->largestParam, *this);
	this->A_bar = findA_bar(this->parent, _params, this->largestParam, *this);

	// parts that will make up final resulting expression after
	// performing all the steps
	this->overallNumFactor_bar = findOverallNumFactor_bar(this->parent);
	this->overallFactor2_bar = findOverallFactor2_bar(this->parent);
	this->numeratorFactor_bar = findNumeratorFactor_bar(this->parent, _params,
			this->largestParam, *this);
	this->U_bar = findU_bar(this->parent, _params, this->largestParam, *this);
	this->F_bar = findF_bar(this->parent, _params, this->largestParam, *this);

	factor_tl_ifPossible(*this, this->largestParam);

	this->y1 = _parent.getY1();
	this->y2 = _parent.getY2();
}

// getters
Sector SubSector::getParent() const {
	return this->parent;
}
GiNaC::symbol SubSector::getLargestParam() const {
	return this->largestParam;
}

// setters
void SubSector::setParent(Sector _parent) {
	this->parent = _parent;
	return;
}
void SubSector::setLargestParam(GiNaC::symbol _param) {
	this->largestParam = _param;
	return;
}

void SubSector::print() const {

	std::cout << "-----------------------------" << "\n";
	std::cout << "Parent sector" << "\n";
	this->parent.print();
	std::cout << "-----------------------------" << "\n";

	std::cout << "SubSector" << "\n";
	std::cout << "Largest parameter = {" << this->largestParam << "} \n";

	std::cout << "Sector variables = {";
	for (int i = 0; i < this->sectorVars.size(); i++) {
		std::cout << this->sectorVars.at(i);
		if (i != this->sectorVars.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}";
	std::cout << "\n";

	std::cout << "Powers a_i_bar  = {";
	for (int i = 0; i < this->powers_a_bar.size(); i++) {
		std::cout << this->powers_a_bar.at(i);
		if (i != this->powers_a_bar.size() - 1) {
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

	std::cout << "tau = " << this->tau << "\n";
	std::cout << "Alpha = " << this->alpha << "\n";
	std::cout << "C_bar = " << this->C_bar << "\n";
	std::cout << "V_bar = " << this->V_bar << "\n";
	std::cout << "A_bar = " << this->A_bar << "\n";
	std::cout << "---------------------------------" << "\n";
	std::cout << "OverallNumFactor_bar = " << this->overallNumFactor_bar
			<< "\n";
	std::cout << "OverallFactor2_bar = " << this->overallFactor2_bar << "\n";
	std::cout << "NumeratorFactor_bar = " << this->numeratorFactor_bar << "\n";
	std::cout << "U_bar = " << this->U_bar << "\n";
	std::cout << "F_bar = " << this->F_bar << "\n";

	std::cout << "" << "\n";

}

// Code used in constructor
std::vector<GiNaC::symbol> SubSector::findSectorVars(Sector _parent,
		GiNaC::symbol _tl) {
	// The function should find are return new subsector variables

	// old sector vars were (t0,t1,...,tn), now we have in argument _params
	// the minimal set of parameters that need to further worked at within
	// sector decomposition, out of which _tl has been chosen as largest
	// new sector variables are again labeled with (t0, t1, ...)

	// e.g. lets have sector vars (t0,t1,t2), out of which (t1,t2) is minimal subset
	// and t1 is taken as largest. Then t1 = t1_prime, t2 = t1_prime*t2_prime
	// new set of parameters - new sector variables is (t0, t1_prime, t2_prime)
	// however we change the names from primes back so new vars are again (t0,t1,t2)

	// thus - since no change in number of sector variables or their names occur
	// this function just returns parent sectors variables

	return _parent.getSectorVars();
}

std::vector<GiNaC::ex> SubSector::findPowers_a_bar(Sector _parent,
		GiNaC::symbol _tl) {
	// powers don't change
	return _parent.getPowers_a_bar();
}

GiNaC::ex SubSector::findC_bar(Sector _parent,
		std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
		SubSector _sector) {
	// The function finds and returns C in new sector variables
	// i.e. after performing the substitution.

	// NOTE: Does not factor out possible _tl even if possible. This will be done in
	// another function (with F and U).

	GiNaC::ex result = _parent.getC_bar();
	std::vector<GiNaC::symbol> minimalSet = _params;

	// do the substitution into sector variables
	for (int i = 0; i < minimalSet.size(); i++) {

		if (minimalSet.at(i) == _tl) {
			continue;
		}

		result = result.subs(minimalSet.at(i) == minimalSet.at(i) * _tl);
	}

	return result;

}

GiNaC::matrix SubSector::findV_bar(Sector _parent,
		std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
		SubSector _sector) {

	GiNaC::matrix result = _parent.getV_bar();
	std::vector<GiNaC::symbol> minimalSet = _params;

	// do the substitution into sector variables
	for (int k = 0; k < minimalSet.size(); k++) {

		if (minimalSet.at(k) == _tl) {
			continue;
		}

		for (int i = 0; i < result.rows(); i++) {
			for (int j = 0; j < result.cols(); j++) {
				result(i, j) = result(i, j).subs(
						minimalSet.at(k) == minimalSet.at(k) * _tl);
			}
		}

	}

	return result;

}

GiNaC::matrix SubSector::findA_bar(Sector _parent,
		std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
		SubSector _sector) {

	GiNaC::matrix result = _parent.getA_bar();
	std::vector<GiNaC::symbol> minimalSet = _params;

	// do the substitution into sector variables
	for (int k = 0; k < minimalSet.size(); k++) {

		if (minimalSet.at(k) == _tl) {
			continue;
		}

		for (int j = 0; j < result.cols(); j++) {
			result(0, j) = result(0, j).subs(
					minimalSet.at(k) == minimalSet.at(k) * _tl);
		}

	}

	return result;
}

GiNaC::ex SubSector::findOverallNumFactor_bar(Sector _parent) {
	// no change there
	return _parent.getOverallNumFactor_bar();
}
GiNaC::ex SubSector::findOverallFactor2_bar(Sector _parent) {
	// no change there
	return _parent.getOverallFactor2_bar();
}

GiNaC::ex SubSector::findNumeratorFactor_bar(Sector _parent,
		std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
		SubSector _sector) {

	GiNaC::ex result = _parent.getNumeratorFactor_bar();
	std::vector<GiNaC::symbol> minimalSet = _params;

	// do the substitution into sector variables
	for (int i = 0; i < minimalSet.size(); i++) {

		if (minimalSet.at(i) == _tl) {
			continue;
		}

		result = result.subs(minimalSet.at(i) == minimalSet.at(i) * _tl);
	}

	// take care of factor that comes up from measure
	result *= GiNaC::pow(_tl, _params.size() - 1);

	return result;
}

GiNaC::ex SubSector::findU_bar(Sector _parent,
		std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
		SubSector _sector) {
	// NOTE: Does not factor out possible _tl even if possible. This will be done in
	// another function (with F and U).

	GiNaC::ex result = _parent.getU_bar();
	std::vector<GiNaC::symbol> minimalSet = _params;

	// do the substitution into sector variables
	for (int i = 0; i < minimalSet.size(); i++) {

		if (minimalSet.at(i) == _tl) {
			continue;
		}

		result = result.subs(minimalSet.at(i) == minimalSet.at(i) * _tl);
	}

	result = result.normal();
	return result;
}

GiNaC::ex SubSector::findF_bar(Sector _parent,
		std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
		SubSector _sector) {

	GiNaC::ex result = _parent.getF_bar();
	std::vector<GiNaC::symbol> minimalSet = _params;

	// do the substitution into sector variables
	for (int i = 0; i < minimalSet.size(); i++) {

		if (minimalSet.at(i) == _tl) {
			continue;
		}

		result = result.subs(minimalSet.at(i) == minimalSet.at(i) * _tl);
	}
	result = result.normal();
	return result;
}

bool SubSector::param_tl_CanBeFactored(GiNaC::ex &_ex_tl_toBeFactoredFrom,
		GiNaC::symbol &_tl) {

	GiNaC::ex aux = _ex_tl_toBeFactoredFrom;
	aux = aux.subs(_tl == 0);
	if (aux == 0) {
		return true;
	}

	return false;
}

int SubSector::findPowerOf_tl_inExpression(GiNaC::ex _ex, GiNaC::symbol _tl) {
	// assumes that tl can be factored out and finds which power of tl can be factored

	int counter = 0;
	GiNaC::ex ex = _ex;

	if (param_tl_CanBeFactored(ex, _tl)) {

		while (param_tl_CanBeFactored(ex, _tl)) {
			ex = ex.diff(_tl, 1);
			counter++;
		}
		return counter;
	} else {
		std::cout
				<< "Actually it is not possible to factor _tl from the expression you put in. "
						"findPowerOf_tl_inExpression()" << std::endl;
		return -1;
	}
}

void SubSector::factor_tl_ifPossible(SubSector &_sector, GiNaC::symbol _tl) {

	GiNaC::ex U_bar = _sector.getU_bar();
	GiNaC::ex F_bar = _sector.getF_bar();
	GiNaC::ex powerOf_tl;
	GiNaC::ex fullPower;

	GiNaC::ex l = _sector.getL();
	GiNaC::ex d = _sector.getD();
	GiNaC::ex alpha = _sector.getAlpha();

	GiNaC::ex newNumeratorFactor = _sector.getNumeratorFactor_bar();

	// try to see if you can factor out tl from U_bar
	if (param_tl_CanBeFactored(U_bar, _tl)) {
		// update numerator Factor
		powerOf_tl = findPowerOf_tl_inExpression(U_bar, _tl);
		fullPower = powerOf_tl * (alpha - (l * d) / 2 - d / 2);

		//newNumeratorFactor *= GiNaC::pow(_tl, fullPower);
		//_sector.setNumeratorFactor_bar(newNumeratorFactor);

		// update U_bar
		_sector.setU_bar((U_bar / (GiNaC::pow(_tl, powerOf_tl))).normal());
	}

	// try to see if you can factor out tl from F_bar
	if (param_tl_CanBeFactored(F_bar, _tl)) {
		// update numerator Factor
		powerOf_tl = findPowerOf_tl_inExpression(F_bar, _tl);
		fullPower += powerOf_tl * (-alpha + l * d / 2);

		_sector.setNumeratorFactor_bar(
				newNumeratorFactor * GiNaC::pow(_tl, fullPower));

		// update F_bar
		_sector.setF_bar((F_bar / (GiNaC::pow(_tl, powerOf_tl))).normal());
	}

	// finally - simplify the new numeratorFactor
	// _sector.setNumeratorFactor_bar(_sector.getNumeratorFactor_bar().normal());

	return;

}

//===========================================================================

//===========================================================================
// Code related to generating primary sectors
std::vector<PrimarySector> generatePrimarySectors(FeynmanParam &_feyn) {
	// Function generates primary sectors for given FeynmanParam given as
	// input argument _feyn.
	// return - vector of primary sectors,
	// Original integral that needs to be solved is now divided into
	// sum over integrals corresponding to individual primary sectors.

	std::vector<PrimarySector> result;
	PrimarySector current;

	for (int i = 0; i < _feyn.getParams().size(); i++) {
		current = PrimarySector(_feyn, _feyn.getParams().at(i));
		result.push_back(current);
	}

	return result;
}

//===========================================================================

//===========================================================================
// class Final Integral

FinalIntegral::FinalIntegral() {

	this->overallNumFactor = 0;
	this->overallFactor2 = 0;
	this->numeratorFactor = 0;
	this->I = 0; // integrand
	this->orderOfHighestPole = 0;
	this->eps = GiNaC::symbol();
	this->tau = GiNaC::symbol();
	this->u0 = GiNaC::symbol();
	this->D0 = GiNaC::symbol();
	this->integVars = { };
	this->y1 = GiNaC::symbol();
	this->y2 = GiNaC::symbol();
	this->alpha = 0;

}
FinalIntegral::FinalIntegral(Sector &_sec) {

	this->integVars = _sec.getSectorVars();
	this->eps = _sec.getEps();
	this->tau = _sec.getTau();
	this->u0 = _sec.getU0();
	this->D0 = _sec.getD0();
	this->orderOfHighestPole = _sec.getL();
	substituteInEpsilon(_sec);
	this->overallNumFactor = _sec.getOverallNumFactor_bar();
	this->overallFactor2 = _sec.getOverallFactor2_bar();
	this->numeratorFactor = _sec.getNumeratorFactor_bar();
	this->I = findI(_sec);
	this->y1 = _sec.getY1();
	this->y2 = _sec.getY2();
	this->alpha = _sec.getAlpha(); //careful, if we have type 2 or 3 then it is already N-2 !!!
}

GiNaC::ex FinalIntegral::findI(Sector _sec) {

	GiNaC::ex result;
	GiNaC::ex numer = _sec.getNumeratorFactor_bar();
	GiNaC::ex U = _sec.getU_bar();
	GiNaC::ex F = _sec.getF_bar();
	GiNaC::ex l = _sec.getL();
	GiNaC::ex alpha = _sec.getAlpha();
	GiNaC::ex eps = this->eps;

	result = GiNaC::pow(U, alpha - (l + 1) * (6 - 2 * eps) / 2)
			* GiNaC::pow(F, -alpha + l * (6 - 2 * eps) / 2);

	return result;
}

// getters
GiNaC::symbol FinalIntegral::getEps() const {
	return this->eps;
}
GiNaC::symbol FinalIntegral::getTau() const {
	return this->tau;
}
GiNaC::symbol FinalIntegral::getU0() const {
	return this->u0;
}
GiNaC::symbol FinalIntegral::getD0() const {
	return this->D0;
}
GiNaC::ex FinalIntegral::getOverallNumFactor() const {
	return this->overallNumFactor;
}
GiNaC::ex FinalIntegral::getOverallFactor2() const {
	return this->overallFactor2;
}
GiNaC::ex FinalIntegral::getNumeratorFactor() const {
	return this->numeratorFactor;
}
GiNaC::ex FinalIntegral::getI() const {
	return this->I;
}
std::vector<GiNaC::symbol> FinalIntegral::getIntegVars() const {
	return this->integVars;
}
GiNaC::ex FinalIntegral::getOrderOfHighestPole() const {
	return this->orderOfHighestPole;
}

GiNaC::symbol FinalIntegral::getY1() const {
	return this->y1;
}
GiNaC::symbol FinalIntegral::getY2() const {
	return this->y2;
}

GiNaC::ex FinalIntegral::getAlpha() const {
	return this->alpha;
}
// setters
void FinalIntegral::setEps(GiNaC::symbol _eps) {
	this->eps = _eps;
	return;
}
void FinalIntegral::setTau(GiNaC::symbol _tau) {
	this->tau = _tau;
	return;
}
void FinalIntegral::setU0(GiNaC::symbol _u0) {
	this->u0 = _u0;
	return;
}
void FinalIntegral::setD0(GiNaC::symbol _D0) {
	this->D0 = _D0;
	return;
}
void FinalIntegral::setOverallNumFactor(GiNaC::ex _overallNumFactor) {
	this->overallNumFactor = _overallNumFactor;
	return;
}
void FinalIntegral::setOverallFactor2(GiNaC::ex _overallFactor2) {
	this->overallFactor2 = _overallFactor2;
	return;
}
void FinalIntegral::setNumeratorFactor(GiNaC::ex _numFact) {
	this->numeratorFactor = _numFact;
	return;
}
void FinalIntegral::setI(GiNaC::ex _I) {
	this->I = _I;
	return;
}
void FinalIntegral::setIntegVars(std::vector<GiNaC::symbol> _vars) {
	this->integVars = _vars;
	return;
}
void FinalIntegral::setOrderOfHighestPole(GiNaC::ex _order) {
	this->orderOfHighestPole = _order;
	return;
}

void FinalIntegral::setY1(GiNaC::symbol _y1) {
	this->y1 = _y1;
	return;
}
void FinalIntegral::setY2(GiNaC::symbol _y2) {
	this->y2 = _y2;
	return;
}
void FinalIntegral::setAlpha(GiNaC::ex _alpha) {
	this->alpha = _alpha;
	return;
}

// print
void FinalIntegral::print() {
	std::cout << "---------------------------------" << "\n";
	std::cout << "Integration variables = {";
	for (int i = 0; i < this->integVars.size(); i++) {
		std::cout << this->integVars.at(i);
		if (i != this->integVars.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}" << "\n";
	std::cout << "OverallNumFactor = " << this->overallNumFactor << "\n";
	std::cout << "OverallFactor2 = " << this->overallFactor2 << "\n";
	std::cout << "NumeratorFactor = " << this->numeratorFactor << "\n";
	std::cout << "I = " << this->I << "\n";

	std::cout << "" << "\n";
	return;
}

//===========================================================================

//===========================================================================
// Code related to finding subsectors (step II.)

bool paramUOrFIsZeroAtLowerIntBoundary(Sector _sector,
		std::vector<GiNaC::symbol> _vars) {
// Auxiliary function. It takes as argument a sector _sector
// and set (actually vector) of parameters _vars. It checks
// if when _vars are set to zero U_bar and F_bar of _sector
// become zero. If so return true. If not zero - return false.

	GiNaC::ex U = _sector.getU_bar().normal();
	GiNaC::ex F = _sector.getF_bar().normal();

	for (int i = 0; i < _vars.size(); i++) {
		U = U.subs(_vars.at(i) == 0);
		F = F.subs(_vars.at(i) == 0);
	}

	if (U == 0 || F == 0) {
		return true;
	}

	return false;
}

std::vector<std::vector<GiNaC::symbol>> findAllCombinationsOfParams(
		std::vector<GiNaC::symbol> _params) {
// Returns vector of all possible subsets of set given by parameters
// e.g. {t0,t1,t2} -> return {{t0}, {t1}, {t2}, {t0,t1}, {t0,t2}, {t1,t2}, {t0,t1,t2}}
// Basically creates power set of set _params that is given as argument.
// This power set is returned in the form of vector.

// NOTE: I found this on the internet, its quite clever

	std::vector<std::vector<GiNaC::symbol>> result = { { } };
	for (int i = 0; i < _params.size(); i++) {
		std::vector<std::vector<GiNaC::symbol>> copy = result;

		for (int j = 0; j < copy.size(); j++) {
			copy.at(j).push_back(_params.at(i));

			result.push_back(copy.at(j));
		}

	}

	std::vector<std::vector<GiNaC::symbol>> finalResult;
// and now get rid of empty set from result
	for (int i = 0; i < result.size(); i++) {
		if (result.at(i).size() == 0) {
			continue;
		}

		finalResult.push_back(result.at(i));
	}

	return finalResult;
}

std::vector<GiNaC::symbol> findMinimalSetOfParams(Sector _sector) {
// Function finds and returns minimal set of sector variables, for which
// when these are set to zero - either U or F become zero

// NOTE: this is not unique. There may be many sets of equal size of params which U or F vanish
// when set to zero, I choose one of those basically at random

	std::vector<std::vector<GiNaC::symbol>> allPossibleSetsOfParams =
			findAllCombinationsOfParams(_sector.getSectorVars());
	GiNaC::ex U = _sector.getU_bar();
	GiNaC::ex F = _sector.getF_bar();

	std::vector<std::vector<GiNaC::symbol>> setsThatMakeUorFZero;

// check which of all the possible sets of params are such, that when all set to zero
// either U_bar or F_bar goes to zero

// for all possible sets
	for (int i = 0; i < allPossibleSetsOfParams.size(); i++) {

		// reset U and F
		U = _sector.getU_bar();
		F = _sector.getF_bar();

		// for all parameters in given set of parameters - set them to zero
		for (int j = 0; j < allPossibleSetsOfParams.at(i).size(); j++) {
			U = U.subs(allPossibleSetsOfParams.at(i).at(j) == 0);
			F = F.subs(allPossibleSetsOfParams.at(i).at(j) == 0);
		}

		// if this did the trick then push back this set of params
		// into vector setsThatMakeUorFZero
		if (U == 0 || F == 0) {
			setsThatMakeUorFZero.push_back(allPossibleSetsOfParams.at(i));
		}
	}

// now find the smallest possible set from setsThatMakeUorFZero
	std::vector<GiNaC::symbol> smallestSet;
	for (int i = 0; i < setsThatMakeUorFZero.size(); i++) {

		if (i == 0) {
			smallestSet = setsThatMakeUorFZero.at(i);
		}

		if (setsThatMakeUorFZero.at(i).size() < smallestSet.size()) {
			smallestSet = setsThatMakeUorFZero.at(i);
		}
	}

	return smallestSet;
}

std::vector<SubSector> generateSubsectors(Sector _sector) {
// Function takes sector as an argument, generates minimal set of parameters
// from which subsectors are constructed. Return - vector of subsectors
// created by this procedure

// NOTE: it is assumed that this function is only performed if _sector needs
// to be further split. (if this was not the case then the original sector
// should to be returned - this is done in function generateFullSetOfSectors()
// that is higher in hierarchy)

// NOTE: this is not final decomposition. Returned subsectors may need to
// be further decomposed.

	std::vector<SubSector> result;
	SubSector current;

// a) find minimal set of parameters
	std::vector<GiNaC::symbol> minimalSet = findMinimalSetOfParams(_sector);

// b) split the sector into subsectors working with sector variables
// corresponding to found minimal set
	for (int i = 0; i < minimalSet.size(); i++) {
		current = SubSector(_sector, minimalSet, minimalSet.at(i));
		result.push_back(current);
	}

	return result;
}

void fullyDecomposeSector(Sector _sector, std::vector<Sector> &_subsectors) {
// Function decomposes sector _sector from argument into
// subsectors such that each subsector from this full decomposition
// does not need further decomposition (both U_bar and F_bar are
// finite as sector variables are set to zero.
// all the generated subsectors are pushed back to vector
// _subsectors (argument of function)

	std::vector<SubSector> subs;

// if sector does not need to be decomposed further then push
// it back to vector _subsectors
	if (!paramUOrFIsZeroAtLowerIntBoundary(_sector, _sector.getSectorVars())) {
		_subsectors.push_back(_sector);
		return;
	}

// if i got here, then further decomposition is needed, so
// generate subsectors
	subs = generateSubsectors(_sector);

// go over all the generated subsectors
	for (int i = 0; i < subs.size(); i++) {

		// if the subsector does not need to be decomposed further
		// push it back to result (vector _subsectors)
		if (!paramUOrFIsZeroAtLowerIntBoundary(subs.at(i),
				subs.at(i).getSectorVars())) {
			_subsectors.push_back(subs.at(i));
		} else {
			// if it is decomposable - use recursion to fully decomopose
			// it and push back to result
			fullyDecomposeSector(subs.at(i), _subsectors);
		}
	}

	return;
}

std::vector<Sector> generateFullSetOfSectors(FeynmanParam _feyn) {
// Function performs sector decomposition. First it divides FeynmanParam
// _feyn into primary sectors, then it further decomposes these (if this is
// needed) into subsectors. If for some sector no further decomposition is
// needed, the sector as such gets pushed back into resulting vector.
// return - vector of all the sectors corresponding to given Feynman Parametrization
// all these are then ready for extraction of poles and numerical calculation.

	std::vector<Sector> result;
	std::vector<PrimarySector> primSectors;
	std::vector<Sector> needFurtherDecomp;

// a) generate primary sectors from feynman parametrization
	primSectors = generatePrimarySectors(_feyn);

// b) go over primary sectors - check if they need to be further
// decomposed. If no - push back to result. If yes - push back to
// needFurtherDecomp.
	for (int i = 0; i < primSectors.size(); i++) {
		if (paramUOrFIsZeroAtLowerIntBoundary(primSectors.at(i),
				primSectors.at(i).getSectorVars())) {
			// if further decompostition needed
			needFurtherDecomp.push_back(primSectors.at(i));
		} else {
			// if no further decomposition needed
			result.push_back(primSectors.at(i));
		}
	}

// c) go over needFurtherDecomp and fully decompose
	for (int i = 0; i < needFurtherDecomp.size(); i++) {

		fullyDecomposeSector(needFurtherDecomp.at(i), result);

	}

	return result;
}

std::vector<FinalIntegral> generateSetOfFinalIntegrals(
		std::vector<Sector> _allSectors) {

	std::vector<FinalIntegral> result;

	for (int i = 0; i < _allSectors.size(); i++) {
		result.push_back(FinalIntegral(_allSectors.at(i)));
	}

	return result;
}

// just some tests
void testFullyDecomposeSector(int _num, std::vector<int> &_result) {

//std::vector<Sector> result;
	std::vector<SubSector> subs;
//bool stillNeedDecomp = true;

// if sector does not need to be decomposed further
	if (_num % 2 != 0) {
		_result.push_back(_num);
		return;
	}

//subs = generateSubsectors(_sector);

	_num = _num / 2;

// if the subsector is does not need to be decomposed further
// push it back to result
	if (_num % 2 != 0) {
		_result.push_back(_num);
	} else {
		// if it is decomposable - fully decomopose it and push back to result
		testFullyDecomposeSector(_num, _result);
	}

// TODO test properly

	return;
}

std::vector<int> testTest(std::vector<int> _numbers) {

	std::vector<int> result;

	for (int i = 0; i < _numbers.size(); i++) {
		testFullyDecomposeSector(_numbers.at(i), result);
	}

	return result;

}
//===========================================================================

//===========================================================================
// Code related to extracting poles

void substituteInEpsilon(Sector &_sect) {
	GiNaC::symbol d = _sect.getD();
	GiNaC::symbol eps = _sect.getEps();

	GiNaC::ex overallNumFactor = _sect.getOverallNumFactor_bar();
	GiNaC::ex overallFactor2 = _sect.getOverallFactor2_bar();
	GiNaC::ex numeratorFactor = _sect.getNumeratorFactor_bar();
	GiNaC::ex U = _sect.getU_bar();
	GiNaC::ex F = _sect.getF_bar();

	// actually in overallNumFactor there should be no epsilons from construction
	// this does not hold for part propto ext momentum - type2 and type3
	//if (overallNumFactor.has(eps)) {
	//	std::cout << "Something wrong. Your overallNumFactor has epsilon in it."
	//			<< std::endl;
	//}

	overallNumFactor = overallNumFactor.subs(d == 6 - 2 * eps);
	_sect.setOverallNumFactor_bar(overallNumFactor);

	overallFactor2 = overallFactor2.subs(d == 6 - 2 * eps);
	_sect.setOverallFactor2_bar(overallFactor2);

	numeratorFactor = numeratorFactor.subs(d == 6 - 2 * eps);
	_sect.setNumeratorFactor_bar(numeratorFactor);

	U = U.subs(d == 6 - 2 * eps);
	_sect.setU_bar(U);

	F = F.subs(d == 6 - 2 * eps);
	_sect.setF_bar(F);

	return;
}

GiNaC::ex findPowerOfVar(GiNaC::ex _ex, GiNaC::symbol _var,
		FinalIntegral _finInt) {
	// NOTE: assumes the ex. is in form _ex = _var^{power}
	// retruns ginac expression power
	//
	// CAREFUL: used to return nonsensical things if form of input is not _ex = _var^{power}
	// e.g. take _ex = _var^{power}*otherVar^{power2}
	// so I do power(var)^{power-1}*otherVar^{power2}
	// and then power*otherVar^{power2} != power
	// THIS HAS BEEN SOLVED BY SUBSTITUTING 1 FOR ALL VARIABLES IN RESULT AFTER
	// DIFFERENTIATION - WORKS IF _ex ONLY CONTAINS TERMS OF TYPE (_var_i)^{power_i}

	GiNaC::ex result = _ex;

	result = result.diff(_var, 1);

	for (int i = 0; i < _finInt.getIntegVars().size(); i++) {
		result = result.subs(_finInt.getIntegVars().at(i) == 1);
	}

	return result;
}

GiNaC::ex findAj(GiNaC::ex _ex, GiNaC::symbol _var, GiNaC::symbol _eps,
		FinalIntegral _finInt) {
	// NOTE: assumes expression _ex can be written as _ex = _var^{Aj - Bj*_eps}.
	// Function returns Aj.
	// CAREFULL: it returns nonsensical things if input is not of form _ex = _var^{Aj - Bj*_eps}.
	// see also note and example of what can go wrong in function findPowerOfVar()

	GiNaC::ex result = findPowerOfVar(_ex, _var, _finInt);

	result = result.subs(_eps == 0);

	return result;

}

GiNaC::ex findBj(GiNaC::ex _ex, GiNaC::symbol _var, GiNaC::symbol _eps,
		FinalIntegral _finInt) {
	// NOTE: assumes expression _ex can be written as _ex = _var^{Aj - Bj*_eps}.
	// Function returns Bj.
	// CAREFULL: it returns nonsensical things if input is not of form _ex = _var^{Aj - Bj*_eps}.
	// see also note and example of what can go wrong in function findPowerOfVar()

	GiNaC::ex result = findPowerOfVar(_ex, _var, _finInt);

	result = -result.diff(_eps, 1);

	return result;

}

void doSeriesInVarUpToPower(FinalIntegral &_finInt, GiNaC::symbol _var,
		GiNaC::ex _upTo) {
	// Auxiliary function used in extractPoles(). Performs partialy step III.
	// of article by Binoth - extracting poles.
	//
	// Assumes integrand of the type _var^{Aj-Bj*eps} * I
	// it does taylor expansion up to order _var^(_upTo) of
	// this integrand so as to extract poles. (the _upTo should be
	// of course |Aj|-1) to fully extract poles

	// It is void function but what it does is it directly
	// changes final sector(FinalIntegral object) _finInt from argument
	// to implement the changes I -> new I set, as well as
	// numeratorFactor is updated
	//
	// Definitions used - article by Binoth:
	// R = I(t,eps) - sum_{p=0}^{_upTo} I^{(p)}(0,eps)*(_var^p) /p!
	// newI = sum_{p=0}^{_upTo} 1/(Aj + p+1 - Bj*eps) * I^{(p)}(0,eps)/p!
	// + _var^{Aj-Bj*eps} R

	GiNaC::symbol eps = _finInt.getEps();
	GiNaC::ex numerFactor = _finInt.getNumeratorFactor();

	GiNaC::ex nf = _finInt.getNumeratorFactor();

	GiNaC::ex R = _finInt.getI();
	GiNaC::ex newI = 0;
	GiNaC::ex newNumerFactor = 0;
	GiNaC::ex currentTerm;
	GiNaC::ex currentCoefInExpansion;

	GiNaC::ex Aj = findAj(nf, _var, eps, _finInt);
	GiNaC::ex Bj = findBj(nf, _var, eps, _finInt);

	int counter = 0;
	GiNaC::ex ginacCount = 0;

	while (ginacCount <= _upTo) {
		currentCoefInExpansion = _finInt.getI().diff(_var, counter);
		currentCoefInExpansion = currentCoefInExpansion.subs(_var == 0);

		newI += (1 / (Aj + ginacCount + 1 - Bj * eps))
				* (currentCoefInExpansion / (GiNaC::factorial(counter)));
		R -= (currentCoefInExpansion * (GiNaC::pow(_var, counter)))
				/ (GiNaC::factorial(counter));

		counter++;
		ginacCount = counter;
	}

	newI += ((GiNaC::pow(_var, Aj - Bj * eps)) * R);
	newNumerFactor = numerFactor.subs(_var == 1);

	_finInt.setI(newI);
	_finInt.setNumeratorFactor(newNumerFactor);

	return;

}

void extractPoles(FinalIntegral &_finInt) {
	// Function extracts poles from final sectors - FinalIntegral objects.
	// It basically performs step III. of article by Binoth

	std::vector<GiNaC::symbol> integrationVars = _finInt.getIntegVars();
	GiNaC::symbol eps = _finInt.getEps();
	GiNaC::ex numeratorFactor = _finInt.getNumeratorFactor();
	GiNaC::ex currentAj;
	GiNaC::ex currentBj;

	// go through all integration variables
	for (int i = 0; i < integrationVars.size(); i++) {
		//std::cout<<"im here "<< integrationVars.size() << std::endl;
		//std::cout<<"working on integr. var "<< integrationVars.at(i);
		currentAj = findAj(numeratorFactor, integrationVars.at(i), eps,
				_finInt);

		// of Aj >= 0 - no poles - nothing to extract
		if (currentAj >= 0) {
			continue;
		} else {
			// if Aj < 0 - extract poles up to order |Aj|-1
			// but since Aj<0 => |Aj| = -Aj

			doSeriesInVarUpToPower(_finInt, integrationVars.at(i),
					-currentAj - 1);
		}
	}

	// once I got here I have extracted all the poles
	return;
}

void changeOrderInVector(std::vector<GiNaC::ex> &_vector) {

	std::vector<GiNaC::ex> newVector;

	for (int i = _vector.size() - 1; i >= 0; i--) {
		newVector.push_back(_vector.at(i));
	}

	_vector = newVector;
	return;
}

std::vector<GiNaC::ex> findPoleCoefficients(FinalIntegral _finInt) {
	// Function finds and returns a vector of pole coefficients
	// (these need to be further integrated using some package for numerical integration)
	//
	// NOTE: the function should be performed after extractPoles()
	//
	// NOTE: (TODO) so far it does not expand overallFactor (involving gamma function)
	// that is why constant term is also in return vector - it is actually eps^-1 pole
	// if gamma function was also expanded
	//
	// the format of return is as follows, e.g. for integral with highest pole
	// of order eps^{-3} the return will be:
	// {coeff of eps^{-3} , coeff of eps^{-2}, coeff of eps^{-1}}

	std::vector<GiNaC::ex> result;
	GiNaC::symbol eps = _finInt.getEps();
	GiNaC::ex numeratorFactor = _finInt.getNumeratorFactor();
	GiNaC::ex I = _finInt.getI();
	GiNaC::ex currentCoefficient;

	// NOTE: (TODO) -1 because i am not expanding gamma function in overallFactor
	GiNaC::ex orderOfHighestPole = _finInt.getOrderOfHighestPole() - 1;

	// expand in eps
	GiNaC::ex expandedIntegral = (numeratorFactor * I).series(eps, 1);
	// get rid of order term (series + O(n)) - gets rid of O(n)
	expandedIntegral = GiNaC::series_to_poly(expandedIntegral);
	//std::cout<<"im here - finding poles "<<_finInt.getOrderOfHighestPole() << std::endl;
	// isolate pole coefficients
	// to do this we first need to multiply the series by
	// highest pole power, then differentiate appropriate number
	// of times in eps, divide by appropriate factorial and finally set eps=0)
	// e.g. (A_{-2} eps^{-2} + A_{-1} eps^{-1} + A_{0} eps^{0})
	// we multiply by eps^{2} and get
	// (A_{-2} + A_{-1} eps^{+1} + A_{0} eps^{+2}), if we want A_{-2}
	// diff. 0 times, divide by 0! and set eps->0, if we want A_{-1}
	// then diff. 1 times, divide by 1! and set eps-> 0

	//expandedIntegral.coeff(eps,-2);
	//std::cout<<"poly coeff -2" << std::endl;
	//expandedIntegral.coeff(eps,-1);
	//std::cout<<"poly coeff -1" << std::endl;
	//expandedIntegral.coeff(eps,0);
	//std::cout<<"poly coeff 0" << std::endl;
	//expandedIntegral *= GiNaC::pow(eps, orderOfHighestPole);

	int counter = 0;
	GiNaC::ex ginacCount = 0;

	while (ginacCount <= orderOfHighestPole) {

		result.push_back(expandedIntegral.coeff(eps, -counter));

		counter++;
		ginacCount = counter;
	}

	changeOrderInVector(result);

	return result;
}

std::vector<GiNaC::ex> findFinitePartCoefficients(FinalIntegral _finInt,
		int _upToOrder) {
	// Function finds and returns a vector of finite part coefficients from overall
	// eps^0 up to eps^{_upToOrder}
	// (these need to be further integrated using some package for numerical integration)
	//
	// NOTE: the function should be performed after extractPoles()
	//
	// NOTE: (TODO) so far it does not expand overallFactor (involving gamma function)
	// that is why constant term is also in return vector - it is actually eps^-1 pole
	// if gamma function was also expanded
	//
	// the format of return is as follows, e.g. for integral with _upToOrder = 2
	// the return will be:
	// {coeff of eps^{0} , coeff of eps^{1}, coeff of eps^{2}}

	std::vector<GiNaC::ex> result;
	GiNaC::symbol eps = _finInt.getEps();
	GiNaC::ex numeratorFactor = _finInt.getNumeratorFactor();
	GiNaC::ex I = _finInt.getI();
	GiNaC::ex currentCoefficient;

	// expand in eps
	GiNaC::ex expandedIntegral = (numeratorFactor * I).series(eps,
			_upToOrder + 2);
	//std::cout<<expandedIntegral<<std::endl;
	// get rid of order term (series + O(n)) - gets rid of O(n)
	expandedIntegral = GiNaC::series_to_poly(expandedIntegral);

	int counter = 0;
	GiNaC::ex ginacCount = 0;

	while (ginacCount <= _upToOrder) {
		// if these are not poles then push them back until order _upToOrder

		// counter+1 because I am not expanding gamma function in overall factor,
		// so actually eps^1 order in expanded integral is overall eps^0 order etc.
		result.push_back(expandedIntegral.coeff(eps, counter + 1));
		counter++;
		ginacCount = counter;
	}

	return result;
}

//===========================================================================

//===========================================================================
// Aux functions
bool isInsideVector(GiNaC::symbol _toBeCompared,
		std::vector<GiNaC::symbol> _vector) {

	// auxiliary function. Returns true if ginac symbol _toBeCompared
	// is within the vector of ginac symbols _vector
	// false if otherwise

	bool result;

	for (int i = 0; i < _vector.size(); i++) {
		if (_vector.at(i) == _toBeCompared) {
			return true;
		}
	}

	return false;
}

void makeProperOverallFactor2_1loop(FinalIntegral &_int) {
	GiNaC::symbol eps = _int.getEps();
	GiNaC::symbol tau = _int.getTau();
	GiNaC::symbol u0 = _int.getU0();
	GiNaC::symbol D0 = _int.getD0();
	GiNaC::ex oldOvFact = _int.getOverallFactor2();
	GiNaC::ex newOvFact = -1 * GiNaC::pow(4 * GiNaC::Pi, -3 + eps)
			* GiNaC::pow(tau, -eps) * GiNaC::tgamma(eps) * GiNaC::pow(u0, 2);
	_int.setI(((oldOvFact / newOvFact) * _int.getI()).subs(tau == 0));
	_int.setOverallFactor2(newOvFact);
	return;
}

void makeProperOverallFactor2_2loop(FinalIntegral &_int) {
	GiNaC::symbol eps = _int.getEps();
	GiNaC::symbol tau = _int.getTau();
	GiNaC::symbol u0 = _int.getU0();
	GiNaC::symbol D0 = _int.getD0();
	GiNaC::ex oldOvFact = _int.getOverallFactor2();
	GiNaC::ex newOvFact = GiNaC::pow(4 * GiNaC::Pi, -6 + 2 * eps)
			* GiNaC::pow(tau, -2 * eps) * GiNaC::tgamma(2 * eps)
			* GiNaC::pow(u0, 4);
	_int.setI(((oldOvFact / newOvFact) * _int.getI()).subs(tau == 0));
	_int.setOverallFactor2(newOvFact);
	return;
}

void makeProperOverallFactor2_3loop(FinalIntegral &_int) {
	GiNaC::symbol eps = _int.getEps();
	GiNaC::symbol tau = _int.getTau();
	GiNaC::symbol u0 = _int.getU0();
	GiNaC::symbol D0 = _int.getD0();
	GiNaC::ex oldOvFact = _int.getOverallFactor2();

	GiNaC::ex newOvFact = -1 * GiNaC::pow(4 * GiNaC::Pi, -9 + 3 * eps)
			* GiNaC::pow(tau, -3 * eps) * GiNaC::tgamma(3 * eps)
			* GiNaC::pow(u0, 6);
	_int.setI(((oldOvFact / newOvFact) * _int.getI()).subs(tau == 0));
	_int.setOverallFactor2(newOvFact);
	return;
}
//===========================================================================

//===========================================================================
// Putting it all together

std::vector<GiNaC::ex> listCoefsOfDivergentPartPropToExtFreq(Diagram _diag) {
	// Function returns a vector of ginac expressions corresponding to
	// coefficients of poles inside diagram _diag (these need to be further integrated using
	// some numerical method - i.e. gives back integrands!! not calculated integrals.)
	//
	// format of return: e.g. for 3-loop graph returns
	// {coef. of eps^{-2} pole, coef. of eps^{-1} pole, coef. of eps^{0} pole}
	// it is like this because there is gamma function in overall factor ~ eps^{-1}
	// so in reality it is:
	// {coef. of eps^{-3} pole, coef. of eps^{-2} pole, coef. of eps^{-1} pole}

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// derivative w.r.t. external frequency produces another vertex -> get all possible insertions of it
	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> freqPart;
	freqPart = addVertCorrespToFrequencyDeriv(diag);

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of freq. part diagrams: " << freqPart.size()
			<< std::endl;
	// for all frequency parts
	for (int i = 0; i < freqPart.size(); i++) {

		// find all possible time orderings
		orderings = freqPart.at(i).findAllPossibleTimeOrderings("y");

		std::cout << "Num. of orderings in freq. part " << i << ": "
				<< orderings.size() << std::endl;

		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			//std::cout<<"working on ordering " << j<<std::endl;
			// construct feynman parametrization
			param = FeynmanParam(freqPart.at(i), orderings.at(j));

			// useQuadraticFormula to integrate over loop momenta in a
			useQuadraticFormulaToIntegrateOverLoopMomenta(param);
			// factor tau
			factorTauOutsideIntoOverallFactor(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);

			// TODO - do it better
			//overallFactor2 = finalInts.at(0).getOverallFactor2();
			//std::cout<<overallFactor2<<std::endl;

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				//std::cout<<"working on final int " << k << "out of "<< finalInts.size() << std::endl;
				currentInt = finalInts.at(k);
				currentNumFactor = finalInts.at(k).getOverallNumFactor();

				extractPoles(currentInt);
				currentPoleCoefs = findPoleCoefficients(currentInt);

				// add coefs to result
				for (int m = 0; m < currentPoleCoefs.size(); m++) {
					result.at(m) += currentNumFactor * currentPoleCoefs.at(m);
				}

			}

		}

	}

	return result;
}

void findDivergentPartsPropToExtFreq_3loop(Diagram _diag, std::string _path) {
	// Function is void type, very similar to findDivergentPartPropToExtFreq()
	// instead of returning coefficients of poles in a vector of ginac expressions
	// this function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside one numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// derivative w.r.t. external frequency produces another vertex -> get all possible insertions of it
	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> freqPart;
	freqPart = addVertCorrespToFrequencyDeriv(diag);

	//just a check
	//for(int i = 0; i< freqPart.size(); i++){
	//freqPart.at(i).print();
	//}

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of freq. part diagrams: " << freqPart.size()
			<< std::endl;
	// for all frequency parts
	for (int i = 0; i < freqPart.size(); i++) {

		// find all possible time orderings
		orderings = freqPart.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in freq. part " << i << ": "
				<< orderings.size() << std::endl;

		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(freqPart.at(i), orderings.at(j));

			// useQuadraticFormula to integrate over loop momenta in a
			useQuadraticFormulaToIntegrateOverLoopMomenta(param);
			// factor tau
			factorTauOutsideIntoOverallFactor(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);
			//std::cout<<"Num. of final ints in ordering "<<j<<": "<< finalInts.size()<<std::endl;

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);
				currentNumFactor = finalInts.at(k).getOverallNumFactor(); //numerical factor

				//std::cout<<currentInt.getOverallFactor2()<<std::endl;

				// this is for check if all the overallFactor2 are the same in all finalInts
				overallFactors.push_back(currentInt.getOverallFactor2());

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentPoleCoefs = findPoleCoefficients(currentInt);

				// add coefs to result
				for (int m = 0; m < currentPoleCoefs.size(); m++) {
					result.at(m) += currentNumFactor * currentPoleCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}

			}

		}

	}

	//  together wls
	//writeWLSCombinedVegasAndEx_separately(diag, integVariables, result,
	//			finalInts.at(0).getOverallFactor2(), _path);

	// Nparts wls
	//writeWLSCombinedExVegas_Nparts(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);

	// together c vegas
	//writeCVegas_eps0(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path);

	// ostrý
	//-----------------
	writeCombinedCVegasAndEx_3loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);
	//writeCVegas_Eps0_Nparts(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path, 5);
	//-----------------
	return;
}

void findDivergentPartsPropToTau_3loop(Diagram _diag, std::string _path) {
	// Function is void type, very similar to findDivergentPartPropToExtFreq().
	// This function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// derivative w.r.t. tau produces another vertex -> get all possible insertions of it
	// should be done with IPMomRouting and nonzero ext. momenta by analogy to FreqDeriv parts.

	std::vector<Diagram> tauParts;
	tauParts = addVertCorrespToTauDeriv(diag);

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of tau part diagrams: " << tauParts.size() << std::endl;
	// for all tau parts
	for (int i = 0; i < tauParts.size(); i++) {

		// just for a check
		//tauParts.at(i).print();

		// find all possible time orderings
		orderings = tauParts.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in tau part " << i << ": "
				<< orderings.size() << std::endl;

		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(tauParts.at(i), orderings.at(j));

			// useQuadraticFormula to integrate over loop momenta in a
			useQuadraticFormulaToIntegrateOverLoopMomenta(param);
			// factor tau
			factorTauOutsideIntoOverallFactor(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);
			//std::cout<<"Num. of final ints in ordering "<<j<<": "<< finalInts.size()<<std::endl;

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);
				currentNumFactor = finalInts.at(k).getOverallNumFactor(); //numerical factor

				//std::cout<<currentInt.getOverallFactor2()<<std::endl;

				// this is for check if all the overallFactor2 are the same in all finalInts
				overallFactors.push_back(currentInt.getOverallFactor2());

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentPoleCoefs = findPoleCoefficients(currentInt);

				// add coefs to result
				for (int m = 0; m < currentPoleCoefs.size(); m++) {
					result.at(m) += currentNumFactor * currentPoleCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}

			}

		}

	}

	//  together wls
	//writeWLSCombinedVegasAndEx_separately(diag, integVariables, result,
	//			finalInts.at(0).getOverallFactor2(), _path);

	// Nparts wls
	//writeWLSCombinedExVegas_Nparts(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);

	// together c vegas
	//writeCVegas_eps0(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path);

	// ostrý
	//-----------------
	//writeWLSVegas(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);
	writeCombinedCVegasAndEx_3loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);
	//writeCVegas_Eps0_Nparts(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path, 5);
	//-----------------
	return;
}

void findFinitePartsPropToExtFreq(Diagram _diag, std::string _path,
		int _orderOfEpsilon) {

	// Function is void type, very similar to findDivergentPartPropToExtFreq()
	// instead of returning coefficients of poles in a vector of ginac expressions
	// this function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside one numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_orderOfEpsilon + 1, 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentFinitePartCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// derivative w.r.t. external frequency produces another vertex -> get all possible insertions of it
	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> freqPart;
	freqPart = addVertCorrespToFrequencyDeriv(diag);

	//just a check
	//for(int i = 0; i< freqPart.size(); i++){
	//freqPart.at(i).print();
	//}

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of freq. part diagrams: " << freqPart.size()
			<< std::endl;
	// for all frequency parts
	for (int i = 0; i < freqPart.size(); i++) {

		// find all possible time orderings
		orderings = freqPart.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in freq. part " << i << ": "
				<< orderings.size() << std::endl;

		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(freqPart.at(i), orderings.at(j));

			// useQuadraticFormula to integrate over loop momenta in a
			useQuadraticFormulaToIntegrateOverLoopMomenta(param);
			// factor tau
			factorTauOutsideIntoOverallFactor(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);
			//std::cout<<"Num. of final ints in ordering "<<j<<": "<< finalInts.size()<<std::endl;

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);
				currentNumFactor = finalInts.at(k).getOverallNumFactor(); //numerical factor

				//std::cout<<currentInt.getOverallFactor2()<<std::endl;

				// this is for check if all the overallFactor2 are the same in all finalInts
				overallFactors.push_back(currentInt.getOverallFactor2());

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentFinitePartCoefs = findFinitePartCoefficients(currentInt,
						_orderOfEpsilon);

				// add coefs to result
				for (int m = 0; m < currentFinitePartCoefs.size(); m++) {
					result.at(m) += currentNumFactor
							* currentFinitePartCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}

			}

		}

	}

	//  together wls
	//writeWLSCombinedVegasAndEx_separately(diag, integVariables, result,
	//			finalInts.at(0).getOverallFactor2(), _path);

	// Nparts wls
	//writeWLSCombinedExVegas_Nparts(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);

	// together c vegas
	//writeCVegas_eps0(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path);

	// ostrý
	//-----------------
	writeCVegas_Finite_2loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path, _orderOfEpsilon);
	//writeCVegas_Eps0_Nparts(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path, 5);
	//-----------------
	return;
}

void findFinitePartsPropToTau(Diagram _diag, std::string _path,
		int _orderOfEpsilon) {

	std::vector<GiNaC::ex> result(_orderOfEpsilon + 1, 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentFinitePartCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// derivative w.r.t. external frequency produces another vertex -> get all possible insertions of it
	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> tauParts;
	tauParts = addVertCorrespToTauDeriv(diag);

	//just a check
	//for(int i = 0; i< freqPart.size(); i++){
	//freqPart.at(i).print();
	//}

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of tau part diagrams: " << tauParts.size() << std::endl;
	// for all tau parts
	for (int i = 0; i < tauParts.size(); i++) {

		// find all possible time orderings
		orderings = tauParts.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in tau part " << i << ": "
				<< orderings.size() << std::endl;

		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(tauParts.at(i), orderings.at(j));

			// useQuadraticFormula to integrate over loop momenta in a
			useQuadraticFormulaToIntegrateOverLoopMomenta(param);
			// factor tau
			factorTauOutsideIntoOverallFactor(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);
			//std::cout<<"Num. of final ints in ordering "<<j<<": "<< finalInts.size()<<std::endl;

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);
				currentNumFactor = finalInts.at(k).getOverallNumFactor(); //numerical factor

				//std::cout<<currentInt.getOverallFactor2()<<std::endl;

				// this is for check if all the overallFactor2 are the same in all finalInts
				overallFactors.push_back(currentInt.getOverallFactor2());

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentFinitePartCoefs = findFinitePartCoefficients(currentInt,
						_orderOfEpsilon);

				// add coefs to result
				for (int m = 0; m < currentFinitePartCoefs.size(); m++) {
					result.at(m) += currentNumFactor
							* currentFinitePartCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}

			}

		}

	}

	//  together wls
	//writeWLSCombinedVegasAndEx_separately(diag, integVariables, result,
	//			finalInts.at(0).getOverallFactor2(), _path);

	// Nparts wls
	//writeWLSCombinedExVegas_Nparts(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);

	// together c vegas
	//writeCVegas_eps0(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path);

	// ostrý
	//-----------------
	writeCVegas_Finite_2loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path, _orderOfEpsilon);
	//writeCVegas_Eps0_Nparts(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path, 5);
	//-----------------
	return;
}

void findDivergentParts_3pt_2loop(Diagram _diag, std::string _path) {
	// Function is void type, very similar to findDivergentPartPropToExtFreq(). But now it calculates
	// divergent part of 3pt functions in 2 loop.
	// This function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	std::vector<std::vector<Vertex>> orderings;

	// just for a check
	//tauParts.at(i).print();

	// find all possible time orderings
	orderings = diag.findAllPossibleTimeOrderings("n");
	std::cout << "Num. of orderings: " << orderings.size() << std::endl;

	// for every time ordering
	for (int j = 0; j < orderings.size(); j++) {

		// construct feynman parametrization
		param = FeynmanParam(diag, orderings.at(j));

		// useQuadraticFormula to integrate over loop momenta in a
		useQuadraticFormulaToIntegrateOverLoopMomenta(param);
		// factor tau
		factorTauOutsideIntoOverallFactor(param);

		// generate set of sectors corresp. to a which are fully decomposed
		sectors = generateFullSetOfSectors(param);

		// convert final set of sectors into objects FinalIntegral
		finalInts = generateSetOfFinalIntegrals(sectors);
		//std::cout<<"Num. of final ints in ordering "<<j<<": "<< finalInts.size()<<std::endl;

		// for every final integral
		for (int k = 0; k < finalInts.size(); k++) {
			currentInt = finalInts.at(k);
			currentNumFactor = finalInts.at(k).getOverallNumFactor(); //numerical factor

			//std::cout<<currentInt.getOverallFactor2()<<std::endl;

			// this is for check if all the overallFactor2 are the same in all finalInts
			overallFactors.push_back(currentInt.getOverallFactor2());

			// extract poles and find pole coefficients
			extractPoles(currentInt);
			currentPoleCoefs = findPoleCoefficients(currentInt);

			// add coefs to result
			for (int m = 0; m < currentPoleCoefs.size(); m++) {
				result.at(m) += currentNumFactor * currentPoleCoefs.at(m);
			}

			// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
			// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
			if (currentInt.getIntegVars().size() > numOfIntegVars) {
				numOfIntegVars = currentInt.getIntegVars().size();
				integVariables = currentInt.getIntegVars();
			}

		}

	}

	//  together wls
	//writeWLSCombinedVegasAndEx_separately(diag, integVariables, result,
	//			finalInts.at(0).getOverallFactor2(), _path);

	// Nparts wls
	//writeWLSCombinedExVegas_Nparts(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);

	// together c vegas
	//writeCVegas_eps0(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path);

	// ostrý
	//-----------------
	//writeWLSVegas(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);
	writeCombinedCVegasAndEx_2loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);
	//writeCVegas_Eps0_Nparts(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path, 5);
	//-----------------
	return;
}

void findDivergentParts_3pt_3loop(Diagram _diag, std::string _path) {
	// Function is void type, very similar to findDivergentPartPropToExtFreq(). But now it calculates
	// divergent part of 3pt functions in 2 loop.
	// This function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;

	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	std::vector<std::vector<Vertex>> orderings;

	// just for a check
	//tauParts.at(i).print();

	// find all possible time orderings
	orderings = diag.findAllPossibleTimeOrderings("n");
	std::cout << "Num. of orderings: " << orderings.size() << std::endl;

	// for every time ordering
	for (int j = 0; j < orderings.size(); j++) {

		// construct feynman parametrization
		param = FeynmanParam(diag, orderings.at(j));

		// useQuadraticFormula to integrate over loop momenta in a
		useQuadraticFormulaToIntegrateOverLoopMomenta(param);
		// factor tau
		factorTauOutsideIntoOverallFactor(param);

		// generate set of sectors corresp. to a which are fully decomposed
		sectors = generateFullSetOfSectors(param);

		// convert final set of sectors into objects FinalIntegral
		finalInts = generateSetOfFinalIntegrals(sectors);
		//std::cout<<"Num. of final ints in ordering "<<j<<": "<< finalInts.size()<<std::endl;

		// for every final integral
		for (int k = 0; k < finalInts.size(); k++) {
			currentInt = finalInts.at(k);
			currentNumFactor = finalInts.at(k).getOverallNumFactor(); //numerical factor

			//std::cout<<currentInt.getOverallFactor2()<<std::endl;

			// this is for check if all the overallFactor2 are the same in all finalInts
			overallFactors.push_back(currentInt.getOverallFactor2());

			// extract poles and find pole coefficients
			extractPoles(currentInt);
			currentPoleCoefs = findPoleCoefficients(currentInt);

			// add coefs to result
			for (int m = 0; m < currentPoleCoefs.size(); m++) {
				result.at(m) += currentNumFactor * currentPoleCoefs.at(m);
			}

			// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
			// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
			if (currentInt.getIntegVars().size() > numOfIntegVars) {
				numOfIntegVars = currentInt.getIntegVars().size();
				integVariables = currentInt.getIntegVars();
			}

		}

	}

	//  together wls
	//writeWLSCombinedVegasAndEx_separately(diag, integVariables, result,
	//			finalInts.at(0).getOverallFactor2(), _path);

	// Nparts wls
	//writeWLSCombinedExVegas_Nparts(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);

	// together c vegas
	//writeCVegas_eps0(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path);

	// ostrý
	//-----------------
	//writeWLSVegas(diag, integVariables, result, finalInts.at(0).getOverallFactor2(), _path);
	writeCombinedCVegasAndEx_3loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);
	//writeCVegas_Eps0_Nparts(diag, integVariables, result,
	//		finalInts.at(0).getOverallFactor2(), _path, 5);
	//-----------------
	return;
}
//===========================================================================

//===========================================================================
// Part propto external momentum

void findDivergentPartsPropToExtMom_1loop(Diagram _diag, std::string _path) {
	// Function is void type, very similar to findDivergentPartPropToExtFreq()
	// instead of returning coefficients of poles in a vector of ginac expressions
	// this function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside one numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;
	GiNaC::symbol y1;
	GiNaC::symbol y2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;
	//diag.print();
	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// when doing part propto ext. momentum, we don't add any new vertex
	// instead we keep ext. momenta in diagram, integrate over loop
	// momenta using quadratic formula, and then only F is dependent on ext.
	// momentum p0 - actually on p0^2, so you can straightforwardly diff.
	// w.r.t. p0^2 and then set p0 = 0

	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> pPart = { diag };
	std::cout << pPart.size() << std::endl;
	// we are not adding any extra vertex so there is only one part

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of p part diagrams: " << pPart.size() << std::endl;
	for (int i = 0; i < pPart.size(); i++) {

		// find all possible time orderings
		orderings = pPart.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in p part " << i << ": "
				<< orderings.size() << std::endl;

		//std::cout<<"i am here";
		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(pPart.at(i), orderings.at(j));

			useQuadraticFormulaToIntegrateOverLoopMomenta(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);

				// dealing with derivative w.r.t. p0^2
				GiNaC::symbol p = param.getDiag().getExtMomenta().at(0);
				GiNaC::symbol p_sq("p_sq");
				GiNaC::ex newI = currentInt.getI().subs(
						GiNaC::pow(p, 2) == p_sq);
				newI = newI.diff(p_sq, 1);
				newI = newI.subs(p_sq == 0);

				currentInt.setI(newI.normal());

				// now we couldn't have factored tau outside after using
				// quadratic formula straightforwardly, because we had
				// nonzero A matrix in feynman param.
				// Now however, I of FinalIntegral is propto tau^(-numOfLoops*eps)
				// so we can factor it out
				factorTauOutsideIntoOverallFactor_proptoExtMom_1loop(
						currentInt);

				// also, since no vertex was added, the overallFactor2
				// is now dependent on tgamma(-1+numOfLoops*eps) rather than tgamma(numOfLoops*eps)
				// however, the derivative w.r.t. p0^2 created factor (-1+numOfLoops*eps) that
				// is multiplying I of FinalIntegral, so together with gamma function in overallFactor2
				// creates proper tgamma as there should be - following function does that
				fixGamma_proptoExtMom_1loop(currentInt);

				// at this point there should be no nontrivial (1) numFactor
				// however we kept it here just in case
				combineNumFactWithI_proptoExtMom(currentInt);

				overallFactors.push_back(currentInt.getOverallFactor2());
				currentNumFactor = currentInt.getOverallNumFactor(); //numerical factor

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentPoleCoefs = findPoleCoefficients(currentInt);

				// add coefs to result
				for (int m = 0; m < currentPoleCoefs.size(); m++) {
					result.at(m) += currentPoleCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}
				//currentInt.print();
				finalInts.at(k) = currentInt;
			}

		}
	}

	writeCombinedCVegasAndEx_1loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);

	return;
}

void findDivergentPartsPropToExtMom_2loop(Diagram _diag, std::string _path) {

	// Function is void type, very similar to findDivergentPartPropToExtFreq()
	// instead of returning coefficients of poles in a vector of ginac expressions
	// this function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside one numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;
	GiNaC::symbol y1;
	GiNaC::symbol y2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;
	//diag.print();
	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// when doing part propto ext. momentum, we don't add any new vertex
	// instead we keep ext. momenta in diagram, integrate over loop
	// momenta using quadratic formula, and then only F is dependent on ext.
	// momentum p0 - actually on p0^2, so you can straightforwardly diff.
	// w.r.t. p0^2 and then set p0 = 0

	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> pPart = { diag };
	std::cout << pPart.size() << std::endl;
	// we are not adding any extra vertex so there is only one part

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of p part diagrams: " << pPart.size() << std::endl;
	for (int i = 0; i < pPart.size(); i++) {

		// find all possible time orderings
		orderings = pPart.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in p part " << i << ": "
				<< orderings.size() << std::endl;

		//std::cout<<"i am here";
		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(pPart.at(i), orderings.at(j));

			useQuadraticFormulaToIntegrateOverLoopMomenta(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);

				// dealing with derivative w.r.t. p0^2
				GiNaC::symbol p = param.getDiag().getExtMomenta().at(0);
				GiNaC::symbol p_sq("p_sq");
				GiNaC::ex newI = currentInt.getI().subs(
						GiNaC::pow(p, 2) == p_sq);
				newI = newI.diff(p_sq, 1);
				newI = newI.subs(p_sq == 0);

				currentInt.setI(newI.normal());

				// now we couldn't have factored tau outside after using
				// quadratic formula straightforwardly, because we had
				// nonzero A matrix in feynman param.
				// Now however, I of FinalIntegral is propto tau^(-numOfLoops*eps)
				// so we can factor it out
				factorTauOutsideIntoOverallFactor_proptoExtMom_2loop(
						currentInt);

				// also, since no vertex was added, the overallFactor2
				// is now dependent on tgamma(-1+numOfLoops*eps) rather than tgamma(numOfLoops*eps)
				// however, the derivative w.r.t. p0^2 created factor (-1+numOfLoops*eps) that
				// is multiplying I of FinalIntegral, so together with gamma function in overallFactor2
				// creates proper tgamma as there should be - following function does that
				fixGamma_proptoExtMom_2loop(currentInt);

				// at this point there should be no nontrivial (1) numFactor
				// however we kept it here just in case
				combineNumFactWithI_proptoExtMom(currentInt);

				overallFactors.push_back(currentInt.getOverallFactor2());
				currentNumFactor = currentInt.getOverallNumFactor(); //numerical factor

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentPoleCoefs = findPoleCoefficients(currentInt);

				// add coefs to result
				for (int m = 0; m < currentPoleCoefs.size(); m++) {
					result.at(m) += currentPoleCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}
				//currentInt.print();
				finalInts.at(k) = currentInt;
			}

		}
	}

	writeCombinedCVegasAndEx_2loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);

	return;
}

void findDivergentPartsPropToExtMom_3loop(Diagram _diag, std::string _path) {
	// Function is void type, very similar to findDivergentPartPropToExtFreq()
	// instead of returning coefficients of poles in a vector of ginac expressions
	// this function writes output file to directory given by _path
	// the output file is in .wls/.cpp format for further use of wolfram mathematica/CUBA
	// All coefficients from all finalIntegrals are
	// summed and one gigantic integrand is put inside one numerical integration procedure (of course
	// pole of every degree gets one NIntegrate[]/Vegas_calculate;

	std::vector<GiNaC::ex> result(_diag.getNumOfLoops(), 0);
	Diagram diag = _diag;
	FeynmanParam param;
	std::vector<Sector> sectors;
	std::vector<FinalIntegral> finalInts;
	FinalIntegral currentInt;
	std::vector<GiNaC::ex> currentPoleCoefs;
	GiNaC::ex currentNumFactor;
	GiNaC::ex overallFactor2;
	GiNaC::symbol y1;
	GiNaC::symbol y2;

	std::vector<GiNaC::ex> overallFactors;

	std::vector<GiNaC::symbol> integVariables;
	int numOfIntegVars = 0;
	//diag.print();
	// treat propagators type "mP" - new vertices
	dividePropsType_mP(diag);

	// when doing part propto ext. momentum, we don't add any new vertex
	// instead we keep ext. momenta in diagram, integrate over loop
	// momenta using quadratic formula, and then only F is dependent on ext.
	// momentum p0 - actually on p0^2, so you can straightforwardly diff.
	// w.r.t. p0^2 and then set p0 = 0

	// can be only done with IPMomRouting and nonzero ext. momenta
	std::vector<Diagram> pPart = { diag };
	std::cout << pPart.size() << std::endl;
	// we are not adding any extra vertex so there is only one part

	std::vector<std::vector<Vertex>> orderings;

	std::cout << "Num. of p part diagrams: " << pPart.size() << std::endl;
	for (int i = 0; i < pPart.size(); i++) {

		// find all possible time orderings
		orderings = pPart.at(i).findAllPossibleTimeOrderings("n");
		std::cout << "Num. of orderings in p part " << i << ": "
				<< orderings.size() << std::endl;

		//std::cout<<"i am here";
		// for every time ordering
		for (int j = 0; j < orderings.size(); j++) {

			// construct feynman parametrization
			param = FeynmanParam(pPart.at(i), orderings.at(j));

			useQuadraticFormulaToIntegrateOverLoopMomenta(param);

			// generate set of sectors corresp. to a which are fully decomposed
			sectors = generateFullSetOfSectors(param);

			// convert final set of sectors into objects FinalIntegral
			finalInts = generateSetOfFinalIntegrals(sectors);

			// for every final integral
			for (int k = 0; k < finalInts.size(); k++) {
				currentInt = finalInts.at(k);

				// dealing with derivative w.r.t. p0^2
				GiNaC::symbol p = param.getDiag().getExtMomenta().at(0);
				GiNaC::symbol p_sq("p_sq");
				GiNaC::ex newI = currentInt.getI().subs(
						GiNaC::pow(p, 2) == p_sq);
				newI = newI.diff(p_sq, 1);
				newI = newI.subs(p_sq == 0);

				currentInt.setI(newI.normal());

				// now we couldn't have factored tau outside after using
				// quadratic formula straightforwardly, because we had
				// nonzero A matrix in feynman param.
				// Now however, I of FinalIntegral is propto tau^(-numOfLoops*eps)
				// so we can factor it out
				factorTauOutsideIntoOverallFactor_proptoExtMom_3loop(
						currentInt);

				// also, since no vertex was added, the overallFactor2
				// is now dependent on tgamma(-1+numOfLoops*eps) rather than tgamma(numOfLoops*eps)
				// however, the derivative w.r.t. p0^2 created factor (-1+numOfLoops*eps) that
				// is multiplying I of FinalIntegral, so together with gamma function in overallFactor2
				// creates proper tgamma as there should be - following function does that
				fixGamma_proptoExtMom_3loop(currentInt);

				// at this point there should be no nontrivial (1) numFactor
				// however we kept it here just in case
				combineNumFactWithI_proptoExtMom(currentInt);

				overallFactors.push_back(currentInt.getOverallFactor2());
				currentNumFactor = currentInt.getOverallNumFactor(); //numerical factor

				// extract poles and find pole coefficients
				extractPoles(currentInt);
				currentPoleCoefs = findPoleCoefficients(currentInt);

				// add coefs to result
				for (int m = 0; m < currentPoleCoefs.size(); m++) {
					result.at(m) += currentPoleCoefs.at(m);
				}

				// update integration variables - these coeffs need to be further integrated over the integration variables from 0, 1
				// this finds the set of int. variables that need to be fed into mathematica or other numerical integration package
				if (currentInt.getIntegVars().size() > numOfIntegVars) {
					numOfIntegVars = currentInt.getIntegVars().size();
					integVariables = currentInt.getIntegVars();
				}
				//currentInt.print();
				finalInts.at(k) = currentInt;
			}

		}
	}

	writeCombinedCVegasAndEx_3loop(diag, integVariables, result,
			finalInts.at(0).getOverallFactor2(), _path);

	return;
}

void factorTauOutsideIntoOverallFactor_proptoExtMom_1loop(FinalIntegral &_fin) {
// the point is overall factor will eventually have tau^(-numOfLoops*eps)
// however, if we are calculating part proportional to external frequency
// especially type2 and type3 - tau cannot be straightforwardly factored out of
// feynman parametrization as now vector A in feynman param is nonzero
// so some taus remain in I of final integral

	GiNaC::ex newI = _fin.getI().subs(_fin.getTau() == 1);
	_fin.setI(newI);
	_fin.setOverallFactor2(
			_fin.getOverallFactor2()
					* GiNaC::pow(_fin.getTau(), -_fin.getEps()));
	return;
}

void factorTauOutsideIntoOverallFactor_proptoExtMom_2loop(FinalIntegral &_fin) {
// the point is overall factor will eventually have tau^(-numOfLoops*eps)
// however, if we are calculating part proportional to external frequency
// especially type2 and type3 - tau cannot be straightforwardly factored out of
// feynman parametrization as now vector A in feynman param is nonzero
// so some taus remain in I of final integral

	GiNaC::ex newI = _fin.getI().subs(_fin.getTau() == 1);

//std::cout <<_fin.getI().coeff(
//		GiNaC::pow(_fin.getTau(), 1 - 2 * _fin.getEps()), 1)<<std::endl;
// in one loop there appears combination tau^(-1)*tau(1-2*eps)
// which ginac cannot combine so first take care of tau(1-2*eps)

	_fin.setI(newI);
	_fin.setOverallFactor2(
			_fin.getOverallFactor2()
					* GiNaC::pow(_fin.getTau(), -2 * _fin.getEps()));
	return;
}

void factorTauOutsideIntoOverallFactor_proptoExtMom_3loop(FinalIntegral &_fin) {
// the point is overall factor will eventually have tau^(-numOfLoops*eps)
// however, if we are calculating part proportional to external frequency
// especially type2 and type3 - tau cannot be straightforwardly factored out of
// feynman parametrization as now vector A in feynman param is nonzero
// so some taus remain in I of final integral

	GiNaC::ex newI = _fin.getI().subs(_fin.getTau() == 1);
	_fin.setI(newI);
	_fin.setOverallFactor2(
			_fin.getOverallFactor2()
					* GiNaC::pow(_fin.getTau(), -3 * _fin.getEps()));
	return;
}

void fixGamma_proptoExtMom_1loop(FinalIntegral &_fin) {
// the point is there will be tgamma(numOfLoops*eps) in overallFactor2ň
// however, in type2 and type3 because of different feynman param. owing
// to nonzero numerator, there will be tgamma(-1+numOfLoops*eps) in overallFactor2
// and (-1+numOfLoops*eps) factor in I of final integral which together create
// the overall gamma function
// this function combines the two and writes it in proper form into overallFactor2

	GiNaC::ex newI = _fin.getI() / (-1 + _fin.getEps());
// in one loop there appears combination tau^(-1)*tau(1-eps)
// which ginac cannot combine so first take care of tau(1-eps)

	_fin.setI(newI);
	_fin.setOverallFactor2(
			(_fin.getOverallFactor2() / GiNaC::tgamma(-1 + _fin.getEps()))
					* GiNaC::tgamma(_fin.getEps()));
	return;
}

void fixGamma_proptoExtMom_2loop(FinalIntegral &_fin) {
// the point is there will be tgamma(numOfLoops*eps) in overallFactor2ň
// however, in type2 and type3 because of different feynman param. owing
// to nonzero numerator, there will be tgamma(-1+numOfLoops*eps) in overallFactor2
// and (-1+numOfLoops*eps) factor in I of final integral which together create
// the overall gamma function
// this function combines the two and writes it in proper form into overallFactor2

	GiNaC::ex newI = _fin.getI() / (-1 + 2 * _fin.getEps());

	_fin.setI(newI);
	_fin.setOverallFactor2(
			(_fin.getOverallFactor2() / GiNaC::tgamma(-1 + 2 * _fin.getEps()))
					* GiNaC::tgamma(2 * _fin.getEps()));
	return;
}

void fixGamma_proptoExtMom_3loop(FinalIntegral &_fin) {
// the point is there will be tgamma(numOfLoops*eps) in overallFactor2ň
// however, in type2 and type3 because of different feynman param. owing
// to nonzero numerator, there will be tgamma(-1+numOfLoops*eps) in overallFactor2
// and (-1+numOfLoops*eps) factor in I of final integral which together create
// the overall gamma function
// this function combines the two and writes it in proper form into overallFactor2

	GiNaC::ex newI = _fin.getI() / (-1 + 3 * _fin.getEps());

	_fin.setI(newI);
	_fin.setOverallFactor2(
			(_fin.getOverallFactor2() / GiNaC::tgamma(-1 + 3 * _fin.getEps()))
					* GiNaC::tgamma(3 * _fin.getEps()));
	return;
}

void factorD0IntoOverallFactor2_proptoExtMom(FinalIntegral &_fin) {
// there is some factor of D_0 in overallNumFactor or I, we want to get it into
// overallFactor2
//std::cout<<_fin.getI()<<std::endl;
//std::cout<<_fin.getI().subs(D0==0)<<std::endl;
	GiNaC::ex numFactor = _fin.getOverallNumFactor();
	numFactor.expand(); // ginac function used when working with polynomials
	int degree = numFactor.degree(_fin.getD0());
//std::cout<<degree<<std::endl;
	_fin.setOverallNumFactor(
			_fin.getOverallNumFactor().subs(_fin.getD0() == 1));
//_fin.setI(_fin.getI().subs(_fin.getD0()==1));
	_fin.setOverallFactor2(
			_fin.getOverallFactor2() * GiNaC::pow(_fin.getD0(), degree));
	return;
}

void combineNumFactWithI_proptoExtMom(FinalIntegral &_fin) {
//now there can be factor like -(-3+eps)^(-1) in numer factor
// this should be put in I and then poles need be extracted

	_fin.setI(_fin.getI() * _fin.getOverallNumFactor());
	_fin.setOverallNumFactor(1);
	return;
}
//===========================================================================

