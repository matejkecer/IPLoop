/*
 * SectorDecomposition.hpp
 *
 *  Created on: 25. 8. 2023
 *      Author: matej
 */

#ifndef SECTORDECOMPOSITION_HPP_
#define SECTORDECOMPOSITION_HPP_

#include "ElementsAndInput_ginac.hpp"
//#include "SpTreesAndForests.hpp"
#include "TimeOrderingsAndTimeCuts_ginac.hpp"
#include "FeynmanAndAlphaParam.hpp"
#include "OutputWriter.hpp"

class Sector {
	// parent class to primarySector and SubSector
protected:
	FeynmanParam parentFeyn;

	std::vector<GiNaC::symbol> extMomenta;
	GiNaC::symbol tau;
	GiNaC::symbol u0;
	GiNaC::symbol D0;
	GiNaC::symbol y1;
	GiNaC::symbol y2;
	std::vector<GiNaC::symbol> sectorVars; // new sector variables t_i (again labeled t)
	std::vector<GiNaC::ex> powers_a_bar; // powers a_i corresp. to t_i in numerator

	GiNaC::symbol d; //dimension of space (in IP it will be 6-\epsilon)
	GiNaC::symbol eps; // 2*eps = 6-d (in IP)
	GiNaC::ex alpha;
	GiNaC::ex l; // number of loops of original diagram

	GiNaC::ex C_bar; // C in new sector variables t_i
	GiNaC::matrix V_bar; // V in new sector variables t_i
	GiNaC::matrix A_bar; // A in new sector variables t_i

	// parts that will make up final resulting expression after
	// performing all the steps
	GiNaC::ex overallNumFactor_bar;
	GiNaC::ex overallFactor2_bar;
	GiNaC::ex numeratorFactor_bar;
	GiNaC::ex U_bar;
	GiNaC::ex F_bar;

public:
	std::vector<GiNaC::symbol> getExtMomenta() const;
	GiNaC::symbol getTau() const;
	GiNaC::symbol getU0() const;
	GiNaC::symbol getD0() const;
	GiNaC::symbol getY1() const;
	GiNaC::symbol getY2() const;
	std::vector<GiNaC::symbol> getSectorVars() const;
	std::vector<GiNaC::ex> getPowers_a_bar() const;
	GiNaC::symbol getD() const;
	GiNaC::symbol getEps() const;

	GiNaC::ex getAlpha() const;
	GiNaC::ex getL() const;
	GiNaC::ex getC_bar() const;
	GiNaC::matrix getV_bar() const;
	GiNaC::matrix getA_bar() const;

	GiNaC::ex getOverallNumFactor_bar() const;
	GiNaC::ex getOverallFactor2_bar() const;
	GiNaC::ex getNumeratorFactor_bar() const;
	GiNaC::ex getU_bar() const;
	GiNaC::ex getF_bar() const;

	FeynmanParam getParentFeyn() const;


	void setExtMomenta(std::vector<GiNaC::symbol> _extMomenta);
	void setTau(GiNaC::symbol _tau);
	void setU0(GiNaC::symbol _u0);
	void setD0(GiNaC::symbol _D0);
	void setY1(GiNaC::symbol _y1);
	void setY2(GiNaC::symbol _y2);
	void setSectorVars(std::vector<GiNaC::symbol> _sectorVars);
	void setPowers_a_bar(std::vector<GiNaC::ex> _powers_a_bar);
	void setD(GiNaC::symbol _d);
	void setEps(GiNaC::symbol _eps);

	void setAlpha(GiNaC::ex _alpha);
	void setL(GiNaC::ex _l);
	void setC_bar(GiNaC::ex _C_bar);
	void setV_bar(GiNaC::matrix _V_bar);
	void setA_bar(GiNaC::matrix _A_bar);

	void setOverallNumFactor_bar(GiNaC::ex _overallNumFactor_bar);
	void setOverallFactor2_bar(GiNaC::ex _overallFactor2_bar);
	void setNumeratorFactor_bar(GiNaC::ex _numerator_bar);
	void setU_bar(GiNaC::ex _U_bar);
	void setF_bar(GiNaC::ex _F_bar);

	void setParentFeyn(FeynmanParam _parentFeyn);
	// print
	void print() const;

};

class PrimarySector: public Sector {
protected:
	// i will use _bar in attribute names to make clear these are attributes in sector variables
	// without _bar are attributes of FeynmanParam in original Feynman parameters

	//FeynmanParam parentFeyn;
	GiNaC::symbol paramIntegratedOut; // Feynman parameter of parentFeyn that has been
	// integrated out when constructing primary sector (use of delta function)

public:
	PrimarySector();
	PrimarySector(FeynmanParam &_feyn, GiNaC::symbol &paramIntegratedOut_xl);

	// getters
	//FeynmanParam getParenFeyn() const;
	GiNaC::symbol getParamIntegratedOut() const;

	// setters
	//void setParenFeyn(FeynmanParam _parentFeyn);
	void setParamIntegratedOut(GiNaC::symbol _param);

	// print
	void print() const;
private:
	//===========================================================================
	// Code used in constructor
	std::vector<GiNaC::symbol> findSectorVars(FeynmanParam &_feyn,
			GiNaC::symbol &_xl);
	std::vector<GiNaC::ex> findPowers_a_bar(FeynmanParam &_feyn,
			GiNaC::symbol &_xl);
	GiNaC::ex findC_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
			PrimarySector &_sect);
	GiNaC::matrix findV_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
			PrimarySector &_sect);
	GiNaC::matrix findA_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
			PrimarySector &_sect);
	GiNaC::ex findOverallNumFactor_bar(FeynmanParam &_parent);
	GiNaC::ex findOverallFactor2_bar(FeynmanParam &_parent);
	GiNaC::ex findNumeratorFactor_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
			PrimarySector &_sect);
	GiNaC::ex findU_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
			PrimarySector &_sect);
	GiNaC::ex findF_bar(FeynmanParam &_feyn, GiNaC::symbol &_xl,
			PrimarySector &_sect);
	//===========================================================================

};

class SubSector: public Sector {
protected:
	Sector parent; // from which I am making subsectors
	GiNaC::symbol largestParam; // sector parameter of parent that has been
	// taken as largest (substitution to subsector vars was based on it)

public:
	SubSector();
	SubSector(Sector &_parent, std::vector<GiNaC::symbol> _params,
			GiNaC::symbol _largestParam_tl);
	// _params is minimal set of parameters that need be further taken care of
	// within sector decomposition

	// getters
	Sector getParent() const;
	GiNaC::symbol getLargestParam() const;

	// setters
	void setParent(Sector _parent);
	void setLargestParam(GiNaC::symbol _param);

	// print
	void print() const;

private:
	//===========================================================================
	// Code used in constructor
	std::vector<GiNaC::symbol> findSectorVars(Sector _parent,
			GiNaC::symbol _tl);
	std::vector<GiNaC::ex> findPowers_a_bar(Sector _parent, GiNaC::symbol _tl);
	GiNaC::ex findC_bar(Sector _parent, std::vector<GiNaC::symbol> _params,
			GiNaC::symbol _tl, SubSector _sector);
	GiNaC::matrix findV_bar(Sector _parent, std::vector<GiNaC::symbol> _params,
			GiNaC::symbol _tl, SubSector _sector);
	GiNaC::matrix findA_bar(Sector _parent, std::vector<GiNaC::symbol> _params,
			GiNaC::symbol _tl, SubSector _sector);
	GiNaC::ex findOverallNumFactor_bar(Sector _parent);
	GiNaC::ex findOverallFactor2_bar(Sector _parent);
	GiNaC::ex findNumeratorFactor_bar(Sector _parent,
			std::vector<GiNaC::symbol> _params, GiNaC::symbol _tl,
			SubSector _sector);
	GiNaC::ex findU_bar(Sector _parent, std::vector<GiNaC::symbol> _params,
			GiNaC::symbol _tl, SubSector _sector);
	GiNaC::ex findF_bar(Sector _parent, std::vector<GiNaC::symbol> _params,
			GiNaC::symbol _tl, SubSector _sector);
	bool param_tl_CanBeFactored(GiNaC::ex &_ex_tl_toBeFactoredFrom,
			GiNaC::symbol &_tl);
	int findPowerOf_tl_inExpression(GiNaC::ex _ex, GiNaC::symbol _tl);
	void factor_tl_ifPossible(SubSector &_sector, GiNaC::symbol _tl);
	//===========================================================================

};

class FinalIntegral {
	// Is basically an assembled integrand corresponding to a particular fully decomposed sector
	// along with all the integration variables and prefactors
protected:

	GiNaC::ex overallNumFactor;
	GiNaC::ex overallFactor2;
	GiNaC::ex numeratorFactor;
	GiNaC::ex I;
	GiNaC::symbol eps;
	GiNaC::symbol tau;
	GiNaC::symbol u0;
	GiNaC::symbol D0;
	std::vector<GiNaC::symbol> integVars;
	GiNaC::ex orderOfHighestPole;
	//GiNaC::ex integrand;
	//GiNaC::ex R;

	GiNaC::ex alpha;
	GiNaC::symbol y1;
	GiNaC::symbol y2;


public:
	FinalIntegral();
	FinalIntegral(Sector &_sector);

	// getters
	GiNaC::ex getOverallNumFactor() const;
	GiNaC::ex getOverallFactor2() const;
	GiNaC::ex getNumeratorFactor() const;
	GiNaC::ex getI() const;
	GiNaC::symbol getEps() const;
	GiNaC::symbol getTau() const;
	GiNaC::symbol getU0() const;
	GiNaC::symbol getD0() const;
	std::vector<GiNaC::symbol> getIntegVars() const;
	GiNaC::ex getOrderOfHighestPole() const;
	GiNaC::symbol getY1() const;
	GiNaC::symbol getY2() const;
	GiNaC::ex getAlpha() const;

	// setters
	void setOverallNumFactor(GiNaC::ex _overallNumFact);
	void setOverallFactor2(GiNaC::ex _overallFact2);
	void setNumeratorFactor(GiNaC::ex _numFact);
	void setI(GiNaC::ex _I);
	void setEps(GiNaC::symbol _eps);
	void setTau(GiNaC::symbol _tau);
	void setIntegVars(std::vector<GiNaC::symbol> _vars);
	void setOrderOfHighestPole(GiNaC::ex _order);
	void setY1(GiNaC::symbol _y1);
	void setY2(GiNaC::symbol _y2);
	void setAlpha(GiNaC::ex _alpha);
	void setU0(GiNaC::symbol _u0);
	void setD0(GiNaC::symbol _D0);

	// print
	void print();

private:
	//GiNaC::ex findIntegrand(Sector _sec);
	GiNaC::ex findI(Sector _sec);
	//GiNaC::ex findR(Sector _sec);
};

//===========================================================================
// Code related to generating primary sectors (step I.)
std::vector<PrimarySector> generatePrimarySectors(FeynmanParam &_feyn);
//===========================================================================

//===========================================================================
// Code related to finding subsectors (step II.)
// TODO

bool paramUOrFIsZeroAtLowerIntBoundary(Sector _sector,
		std::vector<GiNaC::symbol> _vars);
std::vector<std::vector<GiNaC::symbol>> findAllCombinationsOfParams(
		std::vector<GiNaC::symbol> _params);
std::vector<GiNaC::symbol> findMinimalSetOfParams(Sector _sector);
std::vector<SubSector> generateSubsectors(Sector _sector);
void fullyDecomposeSector(Sector _sector, std::vector<SubSector> &_subsectors);

std::vector<Sector> generateFullSetOfSectors(FeynmanParam _feyn); // works

std::vector<FinalIntegral> generateSetOfFinalIntegrals(
		std::vector<Sector> _allSectors); //test

// just some test functions
void testFullyDecomposeSector(int _num, std::vector<int> &_result);
std::vector<int> testTest(std::vector<int> _numbers);
//===========================================================================

//===========================================================================
// Code related to extracting poles

void substituteInEpsilon(Sector &_sect);
GiNaC::ex findPowerOfVar(GiNaC::ex _ex, GiNaC::symbol _var,
		FinalIntegral _finInt);
GiNaC::ex findAj(GiNaC::ex _ex, GiNaC::symbol _var, GiNaC::symbol _eps,
		FinalIntegral _finInt);
GiNaC::ex findBj(GiNaC::ex _ex, GiNaC::symbol _var, GiNaC::symbol _eps,
		FinalIntegral _finInt);

void doSeriesInVarUpToPower(FinalIntegral &_finInt, GiNaC::symbol _var,
		GiNaC::ex _upTo);

void extractPoles(FinalIntegral &_finInt);

std::vector<GiNaC::ex> findPoleCoefficients(FinalIntegral _finInt);

std::vector<GiNaC::ex> findFinitePartCoefficients(FinalIntegral _finInt, int _upToOrder);

void changeOrderInVector(std::vector<GiNaC::ex> &_vector);
//===========================================================================

//===========================================================================
// Aux functions
bool isInsideVector(GiNaC::symbol _toBeCompared,
		std::vector<GiNaC::symbol> _vector);

void makeProperOverallFactor2_1loop(FinalIntegral &_int);
void makeProperOverallFactor2_2loop(FinalIntegral &_int);
void makeProperOverallFactor2_3loop(FinalIntegral &_int);
//===========================================================================

//===========================================================================
// Putting it all together

std::vector<GiNaC::ex> listCoefsOfDivergentPartPropToExtFreq(Diagram _diag);

//bool OvFactorsAreTheSameEverywhere(std::vector<GiNaC::ex> _ovFactors);

void findDivergentPartsPropToExtFreq_3loop(Diagram _diag,
		std::string _path);

void findDivergentPartsPropToTau_3loop(Diagram _diag,
		std::string _path);

void findFinitePartsPropToExtFreq(Diagram _diag,
		std::string _path, int _orderOfEpsilon);
void findFinitePartsPropToTau(Diagram _diag,
		std::string _path, int _orderOfEpsilon);

void findDivergentParts_3pt_2loop(Diagram _diag, std::string _path);
void findDivergentParts_3pt_3loop(Diagram _diag, std::string _path);

void findDivergentPartsPropToExtMom_1loop(Diagram _diag, std::string _path);
void findDivergentPartsPropToExtMom_2loop(Diagram _diag, std::string _path);
void findDivergentPartsPropToExtMom_3loop(Diagram _diag, std::string _path);

void factorTauOutsideIntoOverallFactor_proptoExtMom_1loop(FinalIntegral &_fin);
void factorTauOutsideIntoOverallFactor_proptoExtMom_2loop(FinalIntegral &_fin);
void factorTauOutsideIntoOverallFactor_proptoExtMom_3loop(FinalIntegral &_fin);

void fixGamma_proptoExtMom_1loop(FinalIntegral &_fin);
void fixGamma_proptoExtMom_2loop(FinalIntegral &_fin);
void fixGamma_proptoExtMom_3loop(FinalIntegral &_fin);

void factorD0IntoOverallFactor2_proptoExtMom(FinalIntegral &_fin);
void combineNumFactWithI_proptoExtMom(FinalIntegral &_fin);
//===========================================================================

#endif /* SECTORDECOMPOSITION_HPP_ */
