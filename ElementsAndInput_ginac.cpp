/*
 * ElementsAndInput.cpp
 *
 *  Created on: Aug 7, 2023
 *      Author: M. Kecer
 */

#include "ElementsAndInput_ginac.hpp"
#include "TimeOrderingsAndTimeCuts_ginac.hpp"

//===========================================================================
//===========================================================================
// class Vertex - Implementation

// constructors
Vertex::Vertex() {
	this->vertName = "";
	this->vertType = "";
}

Vertex::Vertex(std::string _vertName) {
	this->vertName = _vertName;
	this->vertType = "";
}

Vertex::Vertex(std::string _vertName, std::string _type) {
	this->vertName = _vertName;
	this->vertType = _type;
}

// getters
std::string Vertex::getVertName() const {
	return this->vertName;
}

std::string Vertex::getVertType() const {
	return this->vertType;
}

// setters
void Vertex::setVertName(std::string _name) {
	this->vertName = _name;
	return;
}
void Vertex::setVertType(std::string _type) {
	this->vertType = _type;
	return;
}

// print
std::string Vertex::info() const {
	std::stringstream s;
	s << "{" << this->vertName << ", " << this->getVertType() << "}";
	return s.str();
}

void Vertex::print() const {
	std::cout << "{" << this->vertName << ", " << this->getVertType() << "}";
	return;
}

// overloading operators == compare vertices
bool Vertex::operator==(const Vertex &rhs) const {

	if (this->vertName.compare(rhs.vertName) != 0) {
		return false;
	}
	if (this->vertType.compare(rhs.vertType) != 0) {
		std::cout
				<< "Careful! Two of your vertices have same name but different type.";
		return false;
	}

	return true;
}

// overloading operators != compare vertices
bool Vertex::operator!=(const Vertex &rhs) const {
	return !(*this == rhs);
}

// overloading assignment operator = for vertices
Vertex& Vertex::operator=(const Vertex &rhs) {
	if (this == &rhs) {      // Same object?
		return *this;
	}

	this->vertName = rhs.vertName;
	this->vertType = rhs.vertType;
	return *this;
}

//===========================================================================
// class Propagator - Implementation

// constructors
Propagator::Propagator() {
	this->startVert = Vertex();
	this->endVert = Vertex();
	this->momentum = 0; // default value
	this->propType = "";
}

Propagator::Propagator(std::string _end, std::string _start, std::string _mom) :
		startVert(_start), endVert(_end) {
	//GiNaC::symbol momentum(_mom);
	//this->momentum = momentum;
	this->momentum = GiNaC::symbol(_mom);
	this->propType = "";
}

Propagator::Propagator(std::string _end, std::string _start, std::string _mom,
		std::string _type) :
		startVert(_start), endVert(_end) {
	//GiNaC::symbol momentum(_mom);
	//this->momentum = momentum;
	this->momentum = GiNaC::symbol(_mom);
	this->propType = _type;
}

Propagator::Propagator(Vertex _end, Vertex _start, GiNaC::ex _mom) {
	//TODO - do you need the & - for symbol probably yes, for ex ???
	this->endVert = _end;
	this->startVert = _start;
	this->momentum = _mom;
	this->propType = "";
}

Propagator::Propagator(Vertex _end, Vertex _start, GiNaC::ex _mom,
		std::string _type) {
	//TODO - do you need the & - for symbol probably yes, for ex ???
	this->endVert = _end;
	this->startVert = _start;
	this->momentum = _mom;
	this->propType = _type;
}

// getters
//std::string Propagator::getStartVert() const {
//	return this->startVert;
//}
Vertex Propagator::getStartVert() const {
	return this->startVert;
}

Vertex Propagator::getEndVert() const {
	return this->endVert;
}

GiNaC::ex Propagator::getMomentum() const {
	return this->momentum;
}

std::string Propagator::getPropType() const {
	return this->propType;
}

// setters
void Propagator::setStartVert(Vertex _start) {
	this->startVert = _start;
	return;
}
void Propagator::setEndVert(Vertex _end) {
	this->endVert = _end;
	return;
}
void Propagator::setMomentum(GiNaC::ex _mom) {
	this->momentum = _mom;
	return;
}
void Propagator::setPropType(std::string _type) {
	this->propType = _type;
	return;
}

// printer
std::string Propagator::info() const {
	std::stringstream s;
	s << "{" << this->endVert.getVertName() << ", "
			<< this->startVert.getVertName() << ", " << this->momentum << ", "
			<< this->getPropType() << "}";
	return s.str();
}

void Propagator::print() const {
	std::cout << "{" << this->endVert.getVertName() << ", "
			<< this->startVert.getVertName() << ", " << this->momentum << ", "
			<< this->getPropType() << "}";
	return;
}

// overloading == operator
bool Propagator::operator==(const Propagator &rhs) const {

	if (this->endVert.getVertName().compare(rhs.endVert.getVertName()) != 0) {
		return false;
	}

	if ((this->startVert.getVertName()).compare(rhs.startVert.getVertName())
			!= 0) {
		return false;
	}

	if (this->momentum != rhs.momentum) {
		return false;
	}

	if (this->propType != rhs.propType) {
		return false;
	}
	return true;
}

// overloading != operator
bool Propagator::operator!=(const Propagator &rhs) const {
	return !(*this == rhs);
}

// overloading assignment operator = for propagators
Propagator& Propagator::operator=(const Propagator &rhs) {
	if (this == &rhs) {      // Same object?
		return *this;
	}

	this->endVert = rhs.endVert;
	this->startVert = rhs.startVert;
	this->momentum = rhs.momentum;
	this->propType = rhs.propType;
	return *this;
}

//===========================================================================
//===========================================================================
// class Diagram - Implementation

// constructors
Diagram::Diagram() {
	this->name = "";
	this->extPropags = { };
	this->intPropags = { };
	this->vertices = { };
	this->extPoint = { };
	this->numOfLoops = 0;
	this->intMomenta = { };
	this->extMomenta = { };
	this->tau = GiNaC::symbol("t");
	this->u_0 = GiNaC::symbol("u0");
	this->D_0 = GiNaC::symbol("D0");
	this->symmetryFactor = 1;
	this->factorsFromVertsAndProps = 1;
	this->d = GiNaC::symbol("d");
	this->createdNumeratorByDerivative = 1;
	this->indicesOfPropsWhichContrToNumerator = { };
	this->y1 = GiNaC::symbol("y1");
	this->y2 = GiNaC::symbol("y2");
}

Diagram::Diagram(std::string _name, std::vector<Propagator> _propags,
		std::vector<Vertex> _verts, int _symmetryFactor) {
	this->name = _name;
	this->extPropags = findExtPropags(_propags);
	this->intPropags = findIntPropags(_propags);
	this->vertices = findVertices(_verts);
	this->extPoint = findExtPoint(_verts);
	this->numOfLoops = findNumOfLoopsPhi3();
	this->tau = GiNaC::symbol("t");
	this->u_0 = GiNaC::symbol("u0");
	this->D_0 = GiNaC::symbol("D0");
	this->symmetryFactor = _symmetryFactor;
	this->factorsFromVertsAndProps = findFactorFromVertsAndProps(*this);
	this->d = GiNaC::symbol("d");
	//this->momenta = defSymbolMomenta();
	this->createdNumeratorByDerivative = 1;
	this->indicesOfPropsWhichContrToNumerator = { };
	this->y1 = GiNaC::symbol("y1");
	this->y2 = GiNaC::symbol("y2");
}

Diagram::Diagram(std::string _name, std::vector<Propagator> _propags,
		std::vector<Vertex> _verts, int _numOfLoops, int _symmetryFactor) {
	this->name = _name;
	this->extPropags = findExtPropags(_propags);
	this->intPropags = findIntPropags(_propags);
	this->vertices = findVertices(_verts);
	this->extPoint = findExtPoint(_verts);
	this->numOfLoops = _numOfLoops;
	this->tau = GiNaC::symbol("t");
	this->u_0 = GiNaC::symbol("u0");
	this->D_0 = GiNaC::symbol("D0");
	this->symmetryFactor = _symmetryFactor;
	this->factorsFromVertsAndProps = findFactorFromVertsAndProps(*this);
	this->d = GiNaC::symbol("d");
	//this->momenta = defSymbolMomenta();
	this->createdNumeratorByDerivative = 1;
	this->indicesOfPropsWhichContrToNumerator = { };
	this->y1 = GiNaC::symbol("y1");
	this->y2 = GiNaC::symbol("y2");
}

Diagram::Diagram(std::string _nickel, int _symmetryFactor) {
	// constructs diagram from nickel indices, assigns general momenta
	// int mom: q0, ..., qn (for every internal leg distinct symbol)
	// ext mom: s0, ..., sm (for every external leg distinct symbol)
	std::vector<Vertex> verts = findNodesFromNickel(_nickel);
	std::vector<Propagator> props = findEdgesFromNickel(_nickel, verts);
	this->name = _nickel;
	this->extPoint = findExtPoint(verts);
	this->vertices = assignVertexTypes(findVertices(verts), props);
	this->extPropags = addVertTypesToProps(findExtPropags(props),
			this->vertices);
	this->intPropags = addVertTypesToProps(findIntPropags(props),
			this->vertices);
	this->numOfLoops = findNumOfLoopsPhi3();
	this->extMomenta = addGeneralExtMomenta(*this, "n"); // assigns non-zero ext. momenta
	this->intMomenta = addGeneralIntMomenta(*this);
	this->tau = GiNaC::symbol("t");
	this->u_0 = GiNaC::symbol("u0");
	this->D_0 = GiNaC::symbol("D0");
	this->symmetryFactor = _symmetryFactor;
	this->factorsFromVertsAndProps = findFactorFromVertsAndProps(*this);
	this->d = GiNaC::symbol("d");
	this->createdNumeratorByDerivative = 1;
	this->indicesOfPropsWhichContrToNumerator = { };
	this->y1 = GiNaC::symbol("y1");
	this->y2 = GiNaC::symbol("y2");
}

Diagram::Diagram(std::string _nickel, std::string _IPMomRouting,
		std::string _zeroExtMomenta, int _symmetryFactor) {
	// constructs diagram from nickel indices
	// Argument _nickel - nickel indices for the diagram (e.g. "e12|e3|33|:0P_mP_pP|0p_Pp|mP_pP|"
	//
	// Argument _zeroExtMomenta - if "y" - sets ext momenta to 0, if "n" assigns ext mom according
	// to either generalMomRouting or IPMomRouting
	//
	// Argument _IPMomRouting - if "n" - assigns general momenta, if "y" then
	// assigns IP momentum routing - uses spanning trees to assign standard momentum routing where loop momenta "k0, ..., kn"
	// are assigned to "mP" type vertices, ext. momenta are assigned as "p0"
	// and all the other propagators have momenta set based on momentum conservation
	//
	std::vector<Vertex> verts = findNodesFromNickel(_nickel);
	std::vector<Propagator> props = findEdgesFromNickel(_nickel, verts);
	this->name = _nickel;
	this->symmetryFactor = _symmetryFactor;
	this->extPoint = findExtPoint(verts);
	this->vertices = assignVertexTypes(findVertices(verts), props);
	this->extPropags = addVertTypesToProps(findExtPropags(props),
			this->vertices);
	this->intPropags = addVertTypesToProps(findIntPropags(props),
			this->vertices);
	this->numOfLoops = findNumOfLoopsPhi3();
	this->tau = GiNaC::symbol("t");
	this->u_0 = GiNaC::symbol("u0");
	this->D_0 = GiNaC::symbol("D0");
	this->factorsFromVertsAndProps = findFactorFromVertsAndProps(*this);
	this->d = GiNaC::symbol("d");
	this->createdNumeratorByDerivative = 1;
	this->indicesOfPropsWhichContrToNumerator = { };
	this->y1 = GiNaC::symbol("y1");
	this->y2 = GiNaC::symbol("y2");

	//this->print();
	// if IPMomRouting "n"
	if (_IPMomRouting.compare("n") == 0) {
		if (_zeroExtMomenta.compare("y") == 0) {
			this->extMomenta = addGeneralExtMomenta(*this, "y");
			this->intMomenta = addGeneralIntMomenta(*this);
		} else {
			if (_zeroExtMomenta.compare("n") == 0) {
				this->extMomenta = addGeneralExtMomenta(*this, "n");
				this->intMomenta = addGeneralIntMomenta(*this);
			} else {
				std::cout
						<< "Something went wrong in the constructor. Probably wrong input."
						<< std::endl;
			}
		}

	} else {
		// if IPMomRouting "y"
		if (_IPMomRouting.compare("y") == 0) {
			if (_zeroExtMomenta.compare("y") == 0) {
				this->extMomenta = addGeneralExtMomenta(*this, "y"); // this needs to be done first
				this->intMomenta = addGeneralIntMomenta(*this); // this needs to be done first
				// now ext before int - see functions

				//this->print();
				this->extMomenta = assignExtMom_IPMomRouting(*this, "y");
				this->intMomenta = assignIntMom_IPMomRouting(*this);

			} else {
				if (_zeroExtMomenta.compare("n") == 0) {
					this->extMomenta = addGeneralExtMomenta(*this, "n"); // this needs to be done first
					this->intMomenta = addGeneralIntMomenta(*this); // this needs to be done first

					//this->print();
					this->extMomenta = assignExtMom_IPMomRouting(*this, "n");
					this->intMomenta = assignIntMom_IPMomRouting(*this);
				} else {
					std::cout
							<< "Something went wrong in the constructor. Probably wrong input."
							<< std::endl;
				}
			}

		} else {
			std::cout
					<< "Something went wrong in the constructor. Probably wrong input."
					<< std::endl;
		}
	}

}

// getters
std::string Diagram::getName() const {
	return this->name;
}

std::vector<Propagator> Diagram::getIntPropags() const {
	return this->intPropags;
}

std::vector<Propagator> Diagram::getExtPropags() const {
	return this->extPropags;
}

std::vector<Vertex> Diagram::getVertices() const {
	return this->vertices;
}

int Diagram::getNumOfLoops() const {
	return this->numOfLoops;
}

std::vector<GiNaC::symbol> Diagram::getIntMomenta() const {
	return this->intMomenta;
}

std::vector<GiNaC::symbol> Diagram::getExtMomenta() const {
	return this->extMomenta;
}

Vertex Diagram::getExtPoint() const {
	return this->extPoint;
}

Propagator Diagram::getIntPropagAtIndex(int _index) const {
	return this->intPropags.at(_index);
}

Propagator Diagram::getExtPropagAtIndex(int _index) const {
	return this->extPropags.at(_index);
}
Vertex Diagram::getVertexAtIndex(int _index) const {
	return this->vertices.at(_index);
}
GiNaC::symbol Diagram::getTau() const {
	return this->tau;
}
GiNaC::symbol Diagram::getU_0() const {
	return this->u_0;
}

GiNaC::symbol Diagram::getD_0() const {
	return this->D_0;
}

int Diagram::getSymmetryFactor() const {
	return this->symmetryFactor;
}

GiNaC::ex Diagram::getFactorsFromVertsAndProps() const {
	return this->factorsFromVertsAndProps;
}
GiNaC::symbol Diagram::getD() const {
	return this->d;
}
GiNaC::ex Diagram::getCreatedNumeratorByDerivative() const {
	return this->createdNumeratorByDerivative;
}
std::vector<int> Diagram::getIndicesOfPropsWhichContrToNumerator() const{
	return this->indicesOfPropsWhichContrToNumerator;
}
GiNaC::symbol Diagram::getY1() const{
	return this->y1;
}
GiNaC::symbol Diagram::getY2() const{
	return this->y2;
}
// setters
void Diagram::setName(std::string _name) {
	this->name = _name;
	return;
}
void Diagram::setExtPropags(std::vector<Propagator> _props) {
	this->extPropags = _props;
	return;
}
void Diagram::setIntPropags(std::vector<Propagator> _props) {
	this->intPropags = _props;
	return;
}
void Diagram::setVertices(std::vector<Vertex> _verts) {
	this->vertices = _verts;
	return;
}
void Diagram::setNumOfLoops(int _num) {
	this->numOfLoops = _num;
	return;
}
void Diagram::setIntMomenta(std::vector<GiNaC::symbol> _momenta) {
	this->intMomenta = _momenta;
	return;
}
void Diagram::setExtMomenta(std::vector<GiNaC::symbol> _momenta) {
	this->extMomenta = _momenta;
	return;
}
void Diagram::setExtPoint(Vertex _extPoint) {
	this->extPoint = _extPoint;
	return;
}
void Diagram::setIntPropagAtIndex(Propagator _propToSet, int _index) {
	this->intPropags.at(_index) = _propToSet;
	return;
}
void Diagram::setExtPropagAtIndex(Propagator _propToSet, int _index) {
	this->extPropags.at(_index) = _propToSet;
	return;
}
void Diagram::setVertexAtIndex(Vertex _vertToSet, int _index) {
	this->vertices.at(_index) = _vertToSet;
	return;
}
void Diagram::setTau(GiNaC::symbol _tau) {
	this->tau = _tau;
	return;
}

void Diagram::setU_0(GiNaC::symbol _u_0) {
	this->u_0 = _u_0;
	return;
}

void Diagram::setD_0(GiNaC::symbol _D_0) {
	this->D_0 = _D_0;
	return;
}

void Diagram::setSymmetryFactor(int _symmetryFactor) {
	this->symmetryFactor = _symmetryFactor;
	return;
}

void Diagram::setFactorsFromVertsAndProps(GiNaC::ex _factor) {
	this->factorsFromVertsAndProps = _factor;
	return;
}

void Diagram::setD(GiNaC::symbol _d) {
	this->d = _d;
	return;
}
void Diagram::setCreatedNumeratorByDerivative(GiNaC::ex _numerator) {
	this->createdNumeratorByDerivative = _numerator;
	return;
}
void Diagram::setIndicesOfPropsWhichContrToNumerator(std::vector<int> _indices){
	this->indicesOfPropsWhichContrToNumerator = _indices;
	return;
}
void Diagram::setY1(GiNaC::symbol _y1){
	this->y1 = _y1;
	return;
}
void Diagram::setY2(GiNaC::symbol _y2){
	this->y2 = _y2;
}

// print
void Diagram::print() const {
	std::cout << "Diagram: " << this->name << "\n";

	std::cout << "Number of Loops: " << this->numOfLoops;
	std::cout << "\n";

	std::cout << "Vertices: ";
	for (int i = 0; i < this->vertices.size(); i++) {
		this->vertices.at(i).print();
		if (i != this->vertices.size() - 1) {
			std::cout << ", ";
		}
	}

	std::cout << "\n";
	std::cout << "External point designation: ";
	this->extPoint.print();
	std::cout << "\n";

	std::cout << "Extermal propagators: ";
	for (int i = 0; i < this->extPropags.size(); i++) {
		this->extPropags.at(i).print();
		std::cout << ", ";
	}
	std::cout << "\n";

	std::cout << "Internal propagators: ";
	for (int i = 0; i < this->intPropags.size(); i++) {
		this->intPropags.at(i).print();
		if (i != this->intPropags.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "\n";

	std::cout << "Internal momenta: {";
	for (int i = 0; i < this->intMomenta.size(); i++) {
		std::cout << this->intMomenta.at(i);
		if (i != this->intMomenta.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}" << "\n";

	std::cout << "External momenta: {";
	for (int i = 0; i < this->extMomenta.size(); i++) {
		std::cout << this->extMomenta.at(i);
		if (i != this->extMomenta.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "}" << "\n";

	std::cout << "Symmetry factor: ";
	std::cout << this->symmetryFactor << "\n";
	//std::cout << "\n";

	std::cout << "Factor from vertices and mP propags: ";
	std::cout << this->factorsFromVertsAndProps << "\n";
	//std::cout << "\n";

	std::cout << "Factor in numerator created by p^2 derivative: ";
	std::cout << this->createdNumeratorByDerivative << "\n";
	//std::cout << "\n";
	std::cout << "Indices of propagators contributing to numerator created by p^2 derivative: {";
	for (int i = 0; i < this->indicesOfPropsWhichContrToNumerator.size(); i++) {
			std::cout << this->indicesOfPropsWhichContrToNumerator.at(i);
			if (i != this->indicesOfPropsWhichContrToNumerator.size() - 1) {
				std::cout << ", ";
			}
		}
	std::cout << "}" << "\n";
	std::cout << "\n";
	return;
}

Diagram& Diagram::operator=(const Diagram &rhs) {
	if (this == &rhs) {      // Same object?
		return *this;
	}

	this->name = rhs.name;
	this->extPropags = rhs.extPropags;
	this->intPropags = rhs.intPropags;
	this->vertices = rhs.vertices;
	this->extPoint = rhs.extPoint; // overkill but still, all diagrams have the same ext point {"e","0"}
	this->numOfLoops = rhs.numOfLoops;
	this->tau = rhs.tau;
	this->u_0 = rhs.u_0;
	this->D_0 = rhs.D_0;
	this->intMomenta = rhs.intMomenta; // todo - check if this will not go baserk - symbols in ginac
	this->extMomenta = rhs.extMomenta;
	this->symmetryFactor = rhs.symmetryFactor;
	this->factorsFromVertsAndProps = rhs.factorsFromVertsAndProps;
	this->d = rhs.d;
	this->createdNumeratorByDerivative = rhs.createdNumeratorByDerivative;
	this->indicesOfPropsWhichContrToNumerator = rhs.indicesOfPropsWhichContrToNumerator;
	this->y1 = rhs.y1;
	this->y2 = rhs.y2;
	return *this;
}

//===========================================================================

//===========================================================================
// functions to be used in diagram constructors
std::vector<Propagator> findExtPropags(std::vector<Propagator> _propToSearch) {
	std::vector<Propagator> extPropags;
	/* Function takes vector<Propagator> and returns only those propags from it
	 * that are external (start or end in external point "E")*/

	for (int i = 0; i < _propToSearch.size(); i++) {
		if (_propToSearch.at(i).getEndVert().getVertName().compare("e") == 0) {
			extPropags.push_back(_propToSearch.at(i));
		} else {
			if (_propToSearch.at(i).getStartVert().getVertName().compare("e")
					== 0) {
				extPropags.push_back(_propToSearch.at(i));
			}
		}

	}
	return extPropags;
}

std::vector<Propagator> findIntPropags(std::vector<Propagator> _propToSearch) {
	/* Function takes vector<Propagator> and returns only those propags from it
	 * that are internal (don't have start or endpoint "E")*/
	std::vector<Propagator> intPropags;

	for (int i = 0; i < _propToSearch.size(); i++) {

		if (_propToSearch.at(i).getEndVert().getVertName().compare("e") != 0
				&& _propToSearch.at(i).getStartVert().getVertName().compare("e")
						!= 0) {
			intPropags.push_back(_propToSearch.at(i));
		}

	}
	return intPropags;
}

int Diagram::findNumOfLoopsPhi3() {
// Function determines number of loops. So far it can only do phi3
// type of vertices

// initialize
	int numOfPropags = this->intPropags.size(); // number of all propags
	std::vector<Vertex> vertices = this->getVertices();
	std::vector<Vertex> phi3Vertices;
	std::vector<Vertex> phi2Vertices;

	for (int i = 0; i < vertices.size(); i++) {
		if (vertices.at(i).getVertType().length() == 3) {
			// phi3 theory type of vertex
			phi3Vertices.push_back(vertices.at(i));
		} else {
			if (vertices.at(i).getVertType().length() == 2) {
				phi2Vertices.push_back(vertices.at(i));
				numOfPropags -= 1;
			} else {
				if (vertices.at(i).getVertType().length() == 0) {
					std::cout
							<< "You didn't specify type of vertices. - findNumOfLoopsPhi3()"
							<< std::endl;
					return -1;
				} else {
					std::cout
							<< "The diagram has non Phi3 vertices. I don't know yet. -findNumOfLoopsPhi3()"
							<< std::endl;
					return -1;
				}
			}
		}
	}

	return (numOfPropags - phi3Vertices.size() + 1);
}

std::vector<Vertex> findVertices(std::vector<Vertex> _vertsToSearch) {
// Function takes all vector<Vertex> _vertsToSearch and returns
// only those vertices that are not external point (vertex {"e", "0"})
// there should always be only one in any given diagram
// if there is multiple -> problems

	std::vector<Vertex> vertices;
	for (int i = 0; i < _vertsToSearch.size(); i++) {
		if (_vertsToSearch.at(i).getVertName().compare("e") != 0) {
			vertices.push_back(_vertsToSearch.at(i));
		}
	}
	return vertices;
}

Vertex findExtPoint(std::vector<Vertex> _vertsToSearch) {
// Function takes all vector<Vertex> _vertsToSearch and returns
// only the external point (vertex {"E","0"})
// there should always be only one in any given diagram
// if there is multiple -> problems

	Vertex extPoint;
	int numOfExtPoints = 0;
	for (int i = 0; i < _vertsToSearch.size(); i++) {
		if (_vertsToSearch.at(i).getVertName().compare("e") == 0) {
			numOfExtPoints += 1;
			extPoint = _vertsToSearch.at(i);
		}
	}
	if (numOfExtPoints > 1) {
		std::cout << "More than one ext. point -> wrong input! -findExtPoint()";
		throw("Problem");
	}

	return extPoint;
}

GiNaC::ex findFactorFromVertsAndProps(Diagram _diag) {
	// Function takes propagators (internal suffice) and vertices (these have to already have assigned
	// types in order for function to work) and assigns factor +D_0 u_0 for every Ppp vert., -D_0 u_0 for every Ppm vert.
	// and D_0 for every Pm propagator

	GiNaC::ex result = 1;
	std::vector<Propagator> intProps = _diag.getIntPropags();
	std::vector<Vertex> verts = _diag.getVertices();

	// find factor from verts.
	Vertex currentVert;
	for (int i = 0; i < verts.size(); i++) {
		currentVert = verts.at(i);
		if (currentVert.getVertType().compare("Ppm") == 0) {
			result *= (-1) * _diag.getU_0() * _diag.getD_0();
		} else {
			if (currentVert.getVertType().compare("PPp") == 0) {
				result *= (+1) * _diag.getU_0() * _diag.getD_0();
			} else {
				std::cout
						<< "You are probably trying to find factor from vertices and propags (u_0, D_0) after the"
								"props have been split (mP -> mM, MP). This might not work. see findFactorFromVertsAndProps(), class Diagram."
						<< std::endl;
			}
		}
	}

	// find factor from propags.
	Propagator currentProp;
	for (int i = 0; i < intProps.size(); i++) {
		currentProp = intProps.at(i);
		if (currentProp.getPropType().compare("mP") == 0) {
			result *= _diag.getD_0();
		}
	}

	return result;
}
//===========================================================================

//===========================================================================
// functions to make program compatible with previous code
std::vector<Propagator> makeDiagramsPropagatorsFromStrings(
		std::vector<std::string> _input) {

	/* Function takes input in form  { "AX(0)" , "XB(k)" , "AD(k)" , ... }
	 * where every vector contains strings, every one of which represents propag.
	 * in the format "startVert endVert momentum"
	 *
	 * TODO - it is stupid format - make it better - right now it only works
	 * for vertex names out of one char
	 *
	 * return: std::vector<Propagator> that gets fetched into diagram constructor
	 * */

// declare variables
	std::vector<Propagator> propagators;
	std::string start;
	std::string end;
	std::string mom;

	for (int i = 0; i < _input.size(); i++) {
		// initialize/reset aux variables
		end = "";
		start = "";
		mom = "";

		end = _input.at(i).at(0);
		start = _input.at(i).at(1);

		for (int j = 2; j < _input.at(i).size(); j++) {
			mom += _input.at(i).at(j);
		}

		Propagator current(end, start, mom);
		propagators.push_back(current);
	}

	return propagators;
}

std::vector<Vertex> makeDiagramsVerticesFromString(std::string _input) {
	/* Function takes input in form  "ABCDEFXYZ"
	 * where every char represents vertex
	 *
	 * TODO - it is stupid format - make it better - right now it only works
	 * for vertex names out of one char
	 *
	 * return: std::vector<Vertex> that gets fetched into diagram constructor
	 * */

// declare variables
	std::vector<Vertex> vertices;
	std::string aux;

// main cycle
	for (int i = 0; i < _input.size(); i++) {
		aux = _input.at(i);
		Vertex current(aux);
		vertices.push_back(current);
	}

	return vertices;
}
//===========================================================================

//===========================================================================
// Code for constructor of Diagram from Nickel indices

std::vector<std::vector<std::string>> divideNickel(std::string _nickelInput) {
// Auxiliary function. It takes nickel index notation e.g. "e12|e3|33|:eP_mP_pP|ep_Pp|mP_pP|"
// and divides it to vector<vector<string>> { {e12, e3, 33}, {eP_mP_pP, ep_Pp, mP_pP} }
// Arguments - nickel index diagram description (string format)
// returns - vector<vector<string>> - diagram description divided into more useful form

	std::vector<std::vector<std::string>> result;
	std::vector<std::string> topologies;
	std::vector<std::string> fields;
	std::string topologyEntry;
	std::string fieldEntry;
	int newEntry = 0;
	int topToFieldTransition = 0;

	for (int i = 0; i < _nickelInput.size(); i++) {

		if (_nickelInput.at(i) == ':') {
			topToFieldTransition = 1;
			continue;
		}

		if (_nickelInput.at(i) == '|') {
			newEntry = 1;
		}

		if (topToFieldTransition == 0) {
			if (newEntry == 0) {
				topologyEntry += _nickelInput.at(i);
			} else {
				// if newEntry == 1
				topologies.push_back(topologyEntry);
				topologyEntry = "";
				newEntry = 0;
				continue;
			}
		} else {
			// if topToFieldTransition == 1
			if (newEntry == 0) {
				fieldEntry += _nickelInput.at(i);
			} else {
				// if newEntry == 1
				fields.push_back(fieldEntry);
				fieldEntry = "";
				newEntry = 0;
				continue;
			}

		}

	}
//test
//std::cout<<topologies.size();
//std::cout<<_nickelInput;
//for(int i = 0; i < 5; i++){
//	std::cout<<topologies.at(i)<<std::endl;
//	std::cout<<fields.at(i)<<std::endl;
//}
	result.push_back(topologies);
	result.push_back(fields);
	return result;
}

std::vector<Vertex> findNodesFromNickel(std::string _nickelInput) {
// Auxiliary function used in constructor of propagator
// creates and returns a vector<Vertex> of diagrams vertices
// NOTE: also external

//std::vector<Vertex> test;
	std::vector<Vertex> result;
	Vertex current;
	int alreadyHaveExt = 0;
	std::string newVertName;

// create vertex with name "0"
	result.push_back(Vertex("0"));

	for (int i = 0; i < _nickelInput.size(); i++) {

		if (_nickelInput.at(i) == ':') {
			return cleanRepeatedVertices(result);
		}
		if (_nickelInput.at(i) == '|') {
			continue;
		}

		if (_nickelInput.at(i) == 'e') {
			// if letter in nickel notation is external leg
			if (alreadyHaveExt == 0) {
				current.setVertName("e");
				current.setVertType("0");
				result.push_back(current);
				alreadyHaveExt = 1;
			}
		} else {
			if (_nickelInput.at(i) >= 48 && _nickelInput.at(i) <= 57) {
				// if current letter in index notation is number of a vertex
				newVertName = _nickelInput.at(i); // implicit change of type char -> string
				current.setVertName(newVertName);
				current.setVertType("");
				result.push_back(current);

			} else {
				std::cout
						<< "Wrong input! - Nickel indices are fishy. -findNodesFromNickel()"
						<< std::endl;
			}
		}
	}

	return cleanRepeatedVertices(result);
}

std::vector<Vertex> cleanRepeatedVertices(std::vector<Vertex> _vertices) {
// Aux function.
// Argument - vector<Vertex> in which there can be repeated entries (e.g. {a,a,b,c,d})
// return - vector<Vertex> cleaned of repetitions (e.g. {a,b,c,d})
//
// NOTE: there can be repeated propags. inside diagram. (Not really if momenta are set
// but at zero momenta {1,2,0,pP}{1,2,0,pP}) there can.
// The propags are the same and yet there are two distinct instances of
// it creating a loop. So you should first always set some momenta.
// but there can only be one vertex of each inside diagram.

	std::vector<Vertex> result;
	Vertex current;
	int repeated = 0;

	for (int i = 0; i < _vertices.size(); i++) {
		current = _vertices.at(i);

		for (int j = i; j < _vertices.size(); j++) {

			if (j == i) {
				continue;
			}

			if (_vertices.at(i) == _vertices.at(j)) {
				repeated = 1;
			}
		}

		if (repeated == 0) {
			result.push_back(current);
		}
		repeated = 0;
	}
	return result;
}

Vertex pickVertexByNameInsideVector(std::string _name,
		std::vector<Vertex> _input) {
// aux function -> picks out and returns a vertex from vector<Vertex>
// that matches the name _name
// if such vertex is not found returns vertex {-1,-1}

	Vertex vert("-1", "-1");

	for (int i = 0; i < _input.size(); i++) {
		if (_input.at(i).getVertName().compare(_name) == 0) {
			return _input.at(i);
		}
	}

	return vert;
}

int findIndexOfVertexInsideVector(Vertex _vert, std::vector<Vertex> _input) {
// function takes as input a vertex and returns at which index
// it is stored inside vector
	int index = 0;

	for (int i = 0; i < _input.size(); i++) {
		if (_vert == _input.at(i)) {
			return i;
		}
	}
	return -1;
}

int findIndexOfPropInsideVector(Propagator _prop,
		std::vector<Propagator> _input) {
	int index = 0;
	for (int i = 0; i < _input.size(); i++) {
		if (_prop == _input.at(i)) {
			return i;
		}
	}
	return -1;
}

std::vector<Propagator> findEdgesFromNickel(std::string _nickel,
		std::vector<Vertex> _verts) {
// The function returns vector<Propagator> of propagators that
// make up diagram encoded by nickel input _nickel
// NOTE: vertices from nickel notation must be already built
// for this to work. It takes these as arguments.
// Therefore use findVerticesFromNickel(std::string _nickelInput) first

	std::vector<std::vector<std::string>> dividedNickel = divideNickel(_nickel);
	std::vector<std::vector<std::string>> dividedFields = divideDividedNickel(
			dividedNickel);
	std::vector<Propagator> result;
	Propagator currentProp;
	char aux;
	std::string endVertString;
	std::string startVertString;

	if (dividedNickel.size() != 2
			|| dividedNickel.at(0).size() != dividedNickel.at(1).size()) {
		std::cout << "Wrong input in divided nickel.";
	}

	for (int i = 0; i < dividedNickel.at(0).size(); i++) {
		aux = 48 + i;
		endVertString = aux;
		currentProp.setEndVert(
				pickVertexByNameInsideVector(endVertString, _verts));
		for (int j = 0; j < dividedNickel.at(0).at(i).size(); j++) {
			startVertString = dividedNickel.at(0).at(i).at(j);
			currentProp.setStartVert(
					pickVertexByNameInsideVector(startVertString, _verts));
			currentProp.setPropType(dividedFields.at(i).at(j));
			result.push_back(currentProp);
		}
	}

//return result;
	return repairCausality(result);
}

std::vector<std::vector<std::string>> divideDividedNickel(
		std::vector<std::vector<std::string>> _divNickel) {
// divideDividedNickel does the following
// dividedNickel looks like: { { e12, e3, 33 }, { 0P_mP_pP, 0p_Pp, mP_pP } }
// what this function does is it divides the second entry of this vector to
// { {0P,mP,pP}, {0p, Pp}, {mP, pP} } -> this is returned

	std::vector<std::vector<std::string>> fields;
	std::vector<std::string> foo; // foo is say {eP,mP,pP} out of { eP_mP_pP, ep_Pp, mP_pP }
	std::string aux;
	std::string result;
//int endOfOneEntry = 0;

	for (int k = 0; k < _divNickel.at(1).size(); k++) {
		//endOfOneEntry = 0;
		for (int m = 0; m < _divNickel.at(1).at(k).size(); m++) {

			if (_divNickel.at(1).at(k).at(m) == '_') {
				//endOfOneEntry = 1;
				foo.push_back(aux);
				aux = "";
				continue;
			}

			//if (endOfOneEntry == 0) {
			aux += _divNickel.at(1).at(k).at(m);
			//}
		}
		foo.push_back(aux);
		aux = "";
		//fields.push_back(foo);
		fields.push_back(foo);
		foo.clear();

	}

	return fields;
}

std::vector<Propagator> repairCausality(std::vector<Propagator> _propags) {
// Function goes through all propagators in vector<Propagator> _propags
// and if it finds propagator where causality is messed up (endVertex is sooner in time then
// startVertex) this is fixed, i.e. "Pm" type -> "mP", "Pp -> pP" and so it is with ext. props

	Vertex aux;

	for (int i = 0; i < _propags.size(); i++) {

		if (_propags.at(i).getPropType().compare("Pm") == 0) {
			aux = _propags.at(i).getEndVert();
			_propags.at(i).setEndVert(_propags.at(i).getStartVert());
			_propags.at(i).setStartVert(aux);
			_propags.at(i).setPropType("mP");
		}
		if (_propags.at(i).getPropType().compare("Pp") == 0) {
			aux = _propags.at(i).getEndVert();
			_propags.at(i).setEndVert(_propags.at(i).getStartVert());
			_propags.at(i).setStartVert(aux);
			_propags.at(i).setPropType("pP");
		}
		// ======================
		// SOMETHING SEEMS FISHY ABOUT NICKEL INDICES (FIELD DESCR.) FOR EXT. POINTS
		// ======================
		//if (_propags.at(i).getPropType().compare("0p") == 0) {
		//	aux = _propags.at(i).getEndVert();
		//	_propags.at(i).setEndVert(_propags.at(i).getStartVert());
		//	_propags.at(i).setStartVert(aux);
		//	_propags.at(i).setPropType("p0");
		//}
		//if (_propags.at(i).getPropType().compare("P0") == 0) {
		//	aux = _propags.at(i).getEndVert();
		//	_propags.at(i).setEndVert(_propags.at(i).getStartVert());
		//	_propags.at(i).setStartVert(aux);
		//	_propags.at(i).setPropType("0P");
		//}
		// ======================

		// ======================
		// IT SEEMS AS IF FIELDS IN NICKEL ARE ALWAYS WRITTEN AS "0P" OR "0P"
		// REGARDLES OF ACCORDANCE WITH TOPOLOGY -> SEE OUR ARTICLE
		// TODO - check Kompaniets
		// ======================
		if (_propags.at(i).getPropType().compare("0p") == 0) {
			if (_propags.at(i).getEndVert().getVertName().compare("e") == 0) {
				aux = _propags.at(i).getEndVert();
				_propags.at(i).setEndVert(_propags.at(i).getStartVert());
				_propags.at(i).setStartVert(aux);
				_propags.at(i).setPropType("p0");
			} else {
				if (_propags.at(i).getStartVert().getVertName().compare("e")
						== 0) {
					_propags.at(i).setPropType("p0");
				}
			}

		}
		if (_propags.at(i).getPropType().compare("0P") == 0) {
			if (_propags.at(i).getEndVert().getVertName().compare("e") != 0) {
				aux = _propags.at(i).getEndVert();
				_propags.at(i).setEndVert(_propags.at(i).getStartVert());
				_propags.at(i).setStartVert(aux);
				_propags.at(i).setPropType("0P");
			}
		}
	}

	return _propags;
}

void mySwap(std::string &_one, std::string &_two) {
	std::string aux = _one;
	_one = _two;
	_two = aux;
	return;
}

std::vector<Vertex> assignVertexTypes(std::vector<Vertex> _vertsWithoutType,
		std::vector<Propagator> _allEdges) {
// Function assigns types to vertices without a set type
// Argument - vector<Vertex> of vertices without set type (but already set names) (-nodes of graph)
// Argument - vector<Propagator> of all propags internal and external (with types of props)
// return - vector<Vertex> of vertices with the type determined from prop. types and set

// NOTE: this should be done on diagram on the level of "P", "p", and "m" fields
// not on the MM types after the split! -> this is because of chosen implementation

// initialize
	std::vector<Vertex> vertsWithTypes;
	Vertex current;
	int numOfFields_m = 0;
	int numOfFields_p = 0;
	int numOfFields_P = 0;

// for all vertices
	for (int i = 0; i < _vertsWithoutType.size(); i++) {
		current = _vertsWithoutType.at(i);

		// for all propagators - determine the type of vertex
		for (int j = 0; j < _allEdges.size(); j++) {
			if (_allEdges.at(j).getStartVert() == current) {
				if (_allEdges.at(j).getPropType().at(1) == 'm') {
					numOfFields_m += 1;
				}
				if (_allEdges.at(j).getPropType().at(1) == 'p') {
					numOfFields_p += 1;
				}
				if (_allEdges.at(j).getPropType().at(1) == 'P') {
					numOfFields_P += 1;
				}
			}

			if (_allEdges.at(j).getEndVert() == current) {
				if (_allEdges.at(j).getPropType().at(0) == 'm') {
					numOfFields_m += 1;
				}
				if (_allEdges.at(j).getPropType().at(0) == 'p') {
					numOfFields_p += 1;
				}
				if (_allEdges.at(j).getPropType().at(0) == 'P') {
					numOfFields_P += 1;
				}
			}
		}

		// assign type of vertex in standard format "PPp", and "Ppm"
		if (numOfFields_m == 1 && numOfFields_p == 1 && numOfFields_P == 1) {
			current.setVertType("Ppm");
		} else {
			if (numOfFields_m == 0 && numOfFields_p == 1
					&& numOfFields_P == 2) {
				current.setVertType("PPp");
			} else {
				std::cout
						<< "Wrong input. You have vertex types that are different from PPp and Ppm. -assignVertexTypes()"
						<< std::endl;
			}
		}

		// reset auxiliary variables
		numOfFields_P = 0;
		numOfFields_p = 0;
		numOfFields_m = 0;

		// push_back the vertex with assigned type to result
		vertsWithTypes.push_back(current);
	}

	return vertsWithTypes;
}

std::vector<Propagator> addVertTypesToProps(
		std::vector<Propagator> _propsWOVertType,
		std::vector<Vertex> _vertsWithType) {

// Function is used in constructor of propagators.
// If we already found nodes of graph and their types (togheter make up vertices)
// and we have edges of graph maybe with propagator type, but as startVert and endVert
// there are still just nodes without type -> add type of vertex into propagators.
	std::vector<Propagator> result;
	Vertex currentVert;
	Propagator currentProp;

	for (int i = 0; i < _propsWOVertType.size(); i++) {
		currentProp = _propsWOVertType.at(i);

		for (int j = 0; j < _vertsWithType.size(); j++) {
			if (_vertsWithType.at(j).getVertName().compare(
					currentProp.getEndVert().getVertName()) == 0) {
				currentProp.setEndVert(_vertsWithType.at(j));
			}

			if (_vertsWithType.at(j).getVertName().compare(
					currentProp.getStartVert().getVertName()) == 0) {
				currentProp.setStartVert(_vertsWithType.at(j));
			}

		}
		result.push_back(currentProp);
	}

	return result;
}

std::vector<GiNaC::symbol> addGeneralIntMomenta(Diagram &_diag) {

// Function is used in constructor of Diagram.
// It assigns a general ginac momentum to every int. propagator
// and saves these symbols to Diagram attribute - intMomenta.
// Internal momenta are designated by letter q.

	std::vector<GiNaC::symbol> result;
	std::stringstream name;
	Propagator propToSet; //help variable

	for (int i = 0; i < _diag.getIntPropags().size(); i++) {
		propToSet = _diag.getIntPropags().at(i);
		name.str("");
		name << "q" << i; // implicitné pretypovanie
		GiNaC::symbol current(name.str());
		// set new momentum
		propToSet.setMomentum(current);
		// set this as new intPropag
		_diag.setIntPropagAtIndex(propToSet, i);
		// save the symbol for newly assigned momentum
		result.push_back(current);
	}
	return result;
}

std::vector<GiNaC::symbol> addGeneralExtMomenta(Diagram &_diag,
		std::string _zeroExtMomenta) {

// Argument _zeroExtMomenta - if "y" - all ext. momenta are set to 0
// if "n" - general ext. momenta will be assigned
//
// Function is used in constructor of Diagram.
// It assigns a general ginac momentum to every ext. propagator
// and saves these symbols to Diagram attribute - extMomenta.
// general ext. moment are designated by letter s

	std::vector<GiNaC::symbol> result;
	std::stringstream name;
	Propagator propToSet; //help variable

// if ext momenta are to be non-zero
	if (_zeroExtMomenta.compare("n") == 0) {
		for (int i = 0; i < _diag.getExtPropags().size(); i++) {
			propToSet = _diag.getExtPropags().at(i);
			name.str("");
			name << "s" << i; // implicitné pretypovanie
			GiNaC::symbol current(name.str());
			// set new momentum
			propToSet.setMomentum(current);
			// set this as new intPropag
			_diag.setExtPropagAtIndex(propToSet, i);
			// save the symbol for newly assigned momentum
			result.push_back(current);
		}
	} else {
		// if ext momenta are to be zero
		if (_zeroExtMomenta.compare("y") == 0) {
			for (int i = 0; i < _diag.getExtPropags().size(); i++) {
				propToSet = _diag.getExtPropags().at(i);
				// set new momentum
				propToSet.setMomentum(0);
				// set this as new intPropag
				_diag.setExtPropagAtIndex(propToSet, i);
			}
		}
	}
	return result;
}

std::vector<GiNaC::symbol> assignExtMom_IPMomRouting(Diagram &_diag,
		std::string _zeroExtMomenta) {
// TODO test
// assign external momenta
// so far it works for 2 point functions

// argument _zeroExtMomenta - possible inputs "y" - diag. is at zero
// ext. momenta or "n" - non-zero ext. momenta are kept.
// return - vector of ginac symobls corresponding to momenta in my diagram _diag.
//
// NOTE: function not only returns the vector of extMomenta but also sets extMomenta
// inside _diag (notice &) (to propagators)
//
// NOTE: at this point the function can only handle 2 point functions.

	std::vector<GiNaC::symbol> result = { };
	Propagator propToSet; // help variable
// if 2 point
	if (_zeroExtMomenta.compare("n") == 0) {
		// if ext. momenta are to be non-zero
		if (_diag.getExtPropags().size() == 2) {
			GiNaC::symbol p0("p0");
			result.push_back(p0);
			for (int i = 0; i < _diag.getExtPropags().size(); i++) {
				propToSet = _diag.getExtPropags().at(i);
				propToSet.setMomentum(p0);
				_diag.setExtPropagAtIndex(propToSet, i);
				//usedProps.push_back(_diag.getExtPropags().at(i));
			}
		}

		if (_diag.getExtPropags().size() == 3) {
			GiNaC::symbol p0("p0");
			GiNaC::symbol p1("p1");
			result.push_back(p0);
			result.push_back(p1);
			for (int i = 0; i < _diag.getExtPropags().size(); i++) {
				propToSet = _diag.getExtPropags().at(i);

				if (i == 0) {
					propToSet.setMomentum(p0);
					_diag.setExtPropagAtIndex(propToSet, i);
				}

				if (i == 1) {
					propToSet.setMomentum(p1);
					_diag.setExtPropagAtIndex(propToSet, i);
				}

				if (i == 2) {
					propToSet.setMomentum(p0 + p1);
					_diag.setExtPropagAtIndex(propToSet, i);
				}
				//usedProps.push_back(_diag.getExtPropags().at(i));
			}
		}
	} else {
		// if ext. momenta are to be zero
		//if (_diag.getExtPropags().size() == 2) {
		for (int i = 0; i < _diag.getExtPropags().size(); i++) {
			propToSet = _diag.getExtPropags().at(i);
			propToSet.setMomentum(0);
			_diag.setExtPropagAtIndex(propToSet, i);
			//usedProps.push_back(_diag.getExtPropags().at(i));
		}
		//}
	}
// if 3 point
// exactly the same thing is done

	return result;
}

std::vector<GiNaC::symbol> assignIntMom_IPMomRouting(Diagram &_diag) {
// Function finds and asigns a momentum routing to diagram _diag.
// What is refered to as IP momentum routing is nothing else than
// standard momentum routing, such that loop momenta are assigned
// to propagators type mP. The function uses spanning trees to find
// such std. routing.
//
// ASSUMES EXT MOMENTA ARE ALREADY SET!!!
//
// return - vector of ginac symobls corresponding to momenta in my diagram _diag.
//
// NOTE: function not only returns the vector of intMomenta but also sets intMomenta
// inside _diag (notice &)
//
// NOTE: at this point the function can only handle 2 point functions.

	Propagator propToSet; // help variable
	std::vector<GiNaC::symbol> result;
	std::string typeOfCurrentDelProp;
	std::stringstream name;
	int aux = 0; // help variable
	std::vector<Propagator> usedProps;
	std::vector<Propagator> newProps;

// assign loop momenta to the "mP" type propags
	for (int i = 0; i < _diag.getIntPropags().size(); i++) {
		propToSet = _diag.getIntPropags().at(i);

		if (_diag.getIntPropags().at(i).getPropType().compare("mP") == 0) {
			name.str("");
			name << "k" << aux; // implicitné pretypovanie
			aux++;
			GiNaC::symbol current(name.str());
			result.push_back(current);
			propToSet.setMomentum(current);
			_diag.setIntPropagAtIndex(propToSet, i);
			usedProps.push_back(_diag.getIntPropags().at(i));
		}

	}

// function assumes external momenta are already set, so ext propags should be
// added as used

	for (int i = 0; i < _diag.getExtPropags().size(); i++) {
		usedProps.push_back(_diag.getExtPropags().at(i));
	}

// assign momenta to the remaining propagators using momentum conservation
// at every vertex
	newProps = assignMomByConservation(_diag, usedProps);

// at this point the propagators with set momenta are kept in newProps
// there are both intProps and extProps

// set new intProps (extProps are already set - see beggining of code)
	_diag.setIntPropags(findIntPropags(newProps));

//test
	if (doesMomConsHold(_diag)) {
		return result;
	} else {
		std::cout
				<< "Something went wrong. Momentum conservation does not hold."
				<< std::endl;
		return {};
	}

	return result;
}

std::vector<Propagator> assignMomByConservation(Diagram &_diag,
		std::vector<Propagator> &_usedProps) {
// Function is used in addIPMomRouting.
// If I have propags. which already have set momenta in _usedProps
// and there is as many props as needed to fully determine the remaining
// propagators momenta in _diag -> assign these remaining propagagtors momenta
// by using momentum conservation at every vertex.
//
// NOTE: right now it works only for vertices type phi3.
// It won't really work well if phi2 type verts are used - do before splitting
// propagators, i.e. before employing dividePropsType_mP() or
// addVertCorrespToFrequencyDeriv() methods.

	std::vector<int> numOfTimesVertsAreUsed(_diag.getVertices().size(), 0);
	std::cout << std::endl;
	Propagator currentProp;
	std::vector<Propagator> propPool;
	GiNaC::ex kirch = 0;
	std::vector<Propagator> allProps = _diag.getExtPropags();
	appendVectors(allProps, _diag.getIntPropags());
	bool isUsed;
	bool isInflowing;
	bool isNotUsedInflowing;
	Vertex currentVert;

// just a little check
	for (int i = 0; i < _diag.getVertices().size(); i++) {
		if (_diag.getVertices().at(i).getVertType().size() != 3) {
			std::cout
					<< "Right now I cannot do diagrams that have non phi3 vertices. -assignMomByConservation()"
					<< std::endl;
			return {};
		}
	}

// repeat the procedure until momenta of all the propagators are determined
	while (_usedProps.size()
			!= _diag.getIntPropags().size() + _diag.getExtPropags().size()) {

		// basically in phi3 type theory, if at some vertex there connect
		// two propagators which already have their momentum set and one
		// that has momentum yet to be set -> this is fully determined by
		// conservation of momentum.

		// so find such vertex where 2 of the connecting propags have already set
		// momenta - first go through all such props with already set momenta
		for (int i = 0; i < _usedProps.size(); i++) {

			// don't count e - that we dont care about
			if (!_usedProps.at(i).getStartVert().getVertName().compare("e")
					== 0) {
				// count usage of the vertex as ++
				numOfTimesVertsAreUsed.at(
						findIndexOfVertexInsideVector(
								_usedProps.at(i).getStartVert(),
								_diag.getVertices()))++;}

				// don't count e - that we dont care about
			if (!_usedProps.at(i).getEndVert().getVertName().compare("e")
					== 0) {
				// count usage of the vertex as ++
				numOfTimesVertsAreUsed.at(
						findIndexOfVertexInsideVector(
								_usedProps.at(i).getEndVert(),
								_diag.getVertices()))++;}
			}

			// if you found such vertex where only one connecting leg is to have set momentum
		for (int i = 0; i < numOfTimesVertsAreUsed.size(); i++) {
			if (numOfTimesVertsAreUsed.at(i) == 2) {

				// set such vertex as currentVertex
				currentVert = _diag.getVertices().at(i);

				// find all props that connect to it
				propPool = getPropsContainingVert(_diag.getVertices().at(i),
						allProps);
				break;
			}
		}

		// reset for next iteration
		numOfTimesVertsAreUsed = std::vector<int>(_diag.getVertices().size(),
				0);

		isUsed = false;
		kirch = 0;
		int indexOfNotUsed = 0;
		int indexInAllProps = 0;

		// now go through all the props in propPool
		// find if given prop has already set momentum
		for (int i = 0; i < propPool.size(); i++) {
			isUsed = false;
			for (int j = 0; j < _usedProps.size(); j++) {
				if (propPool.at(i) == _usedProps.at(j)) {
					isUsed = true;
				}
			}

			// find if the propagator is inflowing to vert or outflowing
			if (propPool.at(i).getEndVert() == currentVert) {
				isInflowing = true;
			} else {
				isInflowing = false;
			}

			if (isUsed) {
				// if it has set momentum and is inflowing - add its momentum
				if (isInflowing) {
					kirch += propPool.at(i).getMomentum();
					continue;
				} else {
					// has set momentum but is outflowing - subtract its
					kirch -= propPool.at(i).getMomentum();
					continue;
				}
			} else {
				// if the propag doesn't have momentum set yet, get its index
				// in propPool
				indexOfNotUsed = i;
				isNotUsedInflowing = isInflowing; // is it inflowing to vertex??
			}
		}

		// reset bool variable befoer next iteration
		isUsed = false; // i am not sure if i am not doing this redundantly - see start of for cycle

		// find where this to-have-momentum-set propagator is in vector allProps
		indexInAllProps = findIndexOfPropInsideVector(
				propPool.at(indexOfNotUsed), allProps);

		if (isNotUsedInflowing) {
			// if it is inflowing - set momentum
			propPool.at(indexOfNotUsed).setMomentum(-kirch);
		} else {
			// if it is outflowing - set momentum
			propPool.at(indexOfNotUsed).setMomentum(kirch);
		}

		// update usedProps in order not to set the same propags. momentum multiple times
		_usedProps.push_back(propPool.at(indexOfNotUsed));

		// upadte allProps - the propagator has now momentum set by mom. cons.
		allProps.at(indexInAllProps) = propPool.at(indexOfNotUsed);
		indexOfNotUsed = 0; //TODO - do i need this??

		// reset variable before next iteration
		propPool.clear();
	}

// return vector with propagators that have already set momenta
// note that here in return there will be both internal and external propags.
	return _usedProps;
}

std::vector<Propagator> getPropsContainingVert(Vertex _vert,
		std::vector<Propagator> _propsToSearch) {
// Function returns all vector of all the propagators that contain vertex
// from argument.

	std::vector<Propagator> result;
	std::vector<Propagator> aux;
	result = getPropagatorsWithEndVertex(_vert, _propsToSearch);
	appendVectors(result, getPropagatorsWithStartVertex(_vert, _propsToSearch));
	return result;
}

bool doesMomConsHold(Diagram _diag) {
// Tester function. It checks if momentum conservation holds
// at every vertex of the diagram.

	std::vector<Propagator> propPool;
	std::vector<Propagator> allProps = _diag.getExtPropags();
	appendVectors(allProps, _diag.getIntPropags());
	Vertex currentVert;
	GiNaC::ex sum = 0;
//bool isInflowing;

// go through all vertices
	for (int i = 0; i < _diag.getVertices().size(); i++) {
		currentVert = _diag.getVertices().at(i);
		propPool = getPropsContainingVert(currentVert, allProps);
		// propPool holds all the propags that start or end in given vertex

		for (int j = 0; j < propPool.size(); j++) {
			// just a test
			if (propPool.size() != 3) {
				std::cout << "Careful. You have non phi3 type vertices."
						<< std::endl;
			}

			if (propPool.at(j).getEndVert() == currentVert) {
				// if the momentum is inflowing
				//isInflowing = true;
				sum += propPool.at(j).getMomentum();
			} else {
				// if the momentum is outflowing
				//isInflowing = false;
				sum -= propPool.at(j).getMomentum();
			}

			//sum += propPool.at(j).getMomentum();
		}

		if (sum == 0) {
			continue;
		} else {
			return false;
		}
		sum = 0;
	}

	return true;
}

//===========================================================================

//===========================================================================
// Code for dividing propags type mP
void dividePropsType_mP(Diagram &_diag) {
// In IP the fields m actually work as addidtional vertices and they should be counted
// when trying to find all possible time orderings.
// This function implements this concept. Propagators type mP are split
// into two propagators and new vertex is created. The vertex type is designated as
// MM while two new propagators will be types pM and MP. (type pM has zero momentum - see article)
// Newly created vertices are given capital letter names starting with A.

	Vertex createdVert;
	Propagator split1;
	Propagator split2;
	std::vector<Propagator> newProps = { };
	std::vector<Vertex> newVerts = { };
	Propagator currentProp;

	int nameNum = 65;
	char name;
	std::string strName = "";

// go through all internal propagators
	for (int i = 0; i < _diag.getIntPropags().size(); i++) {
		currentProp = _diag.getIntPropags().at(i);
		// i want to only divide mP propagators
		if (!_diag.getIntPropags().at(i).getPropType().compare("mP") == 0) {
			continue;
		}

		// when i have mP propagator a create a name for new vertex
		name = nameNum; // implicitné pretypovanie
		nameNum++;
		strName = "";
		strName += name;

		// I create new vertex with the name
		createdVert = Vertex(strName, "MM");

		// the propagator is split into two
		split1 = Propagator(currentProp.getEndVert(), createdVert, 0, "pM");
		split2 = Propagator(createdVert, currentProp.getStartVert(),
				currentProp.getMomentum(), "MP");

		// one of the newly created propagators directly replaces the old one in diagram
		_diag.setIntPropagAtIndex(split2, i);
		newProps.push_back(split1);
		newVerts.push_back(createdVert);
	}

// the others (split1's) are appended
	appendVectors(newProps, _diag.getIntPropags());
	appendVectors(newVerts, _diag.getVertices());
	_diag.setIntPropags(newProps);
	_diag.setVertices(newVerts);

	return;
}

//===========================================================================

//===========================================================================
// Code for getting part proportional to ext. frequency
std::vector<Diagram> addVertCorrespToFrequencyDeriv(Diagram _parentDiag) {
// In IP we are mainly interested in finding out exponent "z", for which
// diverging part of two point functions proportional to ext. frequency needs
// to be calculated. This can be found by performing derivative w.r.t. ext. frequency.
// On graph level in (p,t) representation this corresponds to attaching new unit
// vertex on some internal edge carring ext. frequency. (all such insertions must
// be found). This function takes implements this concept, and returns vector
// of diagrams, each of which has insertion on one of the internal legs
// that carries ext. frequency). Such results can be further looked at to find
// time versions and etc.

// NOTE: at this point it works for 2 point functions only and IPMomRouting
// is assumed (_parentDiagram already has split mP props).

	std::vector<Diagram> result;
	Diagram created = _parentDiag;
	std::stringstream createdDiagName;
	Vertex createdVert;
	Propagator split1;
	Propagator split2;
	std::vector<Propagator> newProps = { };
	std::vector<Vertex> newVerts = { };
	Propagator currentProp;

// where there is ext. momentum there is ext. frequency which we are after
// i.e. the new vertex can only be inserted to internal leg carrying ext. momentum

// now we only programmed IPMomRouting for 2 point functions, so there is only one
// external momentum p0 (in IPMomRouting)
// TODO - more general case
	GiNaC::symbol extMom = created.getExtMomenta().at(0);

// this will be name of newly created vertex (only one needs to be created - one frequency derivative)
	std::string name = "f";

// just some tests to catch possible problems
	if (_parentDiag.getExtPropags().size() > 2) {
		std::cout
				<< "Right now I cannot do diagrams of 3 point and more point functions. -addVertCorrespToFrequencyInt()"
				<< std::endl;
		return result;
	}

// naturally the function can be used only when diagram has non-zero ext. mom. (and ext. freq.)
	for (int i = 0; i < _parentDiag.getExtPropags().size(); i++) {
		if (_parentDiag.getExtPropags().at(i).getMomentum() == 0) {
			std::cout
					<< "Illegal actions. Diagram has ext. momenta set to zero (freq's also)."
							"Therefore derivative w.r.t. ext. frequency gives zero. -addVertCorrsepToFrequencyDeriv()"
					<< std::endl;
			return result;
		}
	}

// since only 2 point functions are allowed and IPMomRouting must be used
	if (_parentDiag.getExtMomenta().size() != 1) {
		std::cout
				<< "Illegal actions. You probably didn't use IPMomRouting which is required"
						"-addVertCorrsepToFrequencyDeriv()" << std::endl;
		return result;
	}

// back to code
// go through all internal propagators
	for (int i = 0; i < _parentDiag.getIntPropags().size(); i++) {
		created = _parentDiag;

		currentProp = created.getIntPropags().at(i);

		// i want to only divide paropagators with ext. momentum (frequency)
		if (!currentProp.getMomentum().has(extMom)) {
			continue;
		}

		// create new vertex
		createdVert = Vertex(name, "pP");

		// split the propagator into two
		split1 = Propagator(currentProp.getEndVert(), createdVert,
				currentProp.getMomentum(), "pP");
		split2 = Propagator(createdVert, currentProp.getStartVert(),
				currentProp.getMomentum(), "pP");

		// one of the newly created props directly replaces old propagator in diagram
		created.setIntPropagAtIndex(split2, i);
		newProps.push_back(split1);
		newVerts.push_back(createdVert);

		// the rest are appended
		appendVectors(newProps, created.getIntPropags());
		appendVectors(newVerts, created.getVertices());
		created.setIntPropags(newProps);
		created.setVertices(newVerts);

		// give returend diagram some different name to differentiate it from parent
		createdDiagName.str("");
		createdDiagName << "fq_" << result.size() << "_"
				<< _parentDiag.getName();
		created.setName(createdDiagName.str());

		// if it is proportional to ext frequency it will not be to ext. momentum
		// therefore ext. mom may be set to zero
		result.push_back(setExtMomToZero(created));

		// reset variables for next iteration
		newProps.clear();
		newVerts.clear();
	}

	return result;
}

Diagram setExtMomToZero(Diagram _inputDiag) {
// Auxiliary function that sets external momentum to zero in all its propagators.
// NOTE: so far it only works for 2 point functions - only one ext. momentum

	Diagram result = _inputDiag;
	Propagator currentProp;

// now we only programmed IPMomRouting for 2 point functions, so there is only one
// external momentum p0 (in IPMomRouting)
// TODO - more general case
	GiNaC::symbol extMom = result.getExtMomenta().at(0);

	GiNaC::ex currentMom;
//std::vector<GiNaC::symbol> newMomenta;

// test - catching possible problems and wrong input
	if (_inputDiag.getExtPropags().size() > 2) {
		std::cout
				<< "Right now I cannot do diagrams of 3 point and more point functions. -setExtMomToZero()"
				<< std::endl;
		return result;
	}

// clear ext. momentum from ext. props
	for (int i = 0; i < _inputDiag.getExtPropags().size(); i++) {
		currentProp = result.getExtPropags().at(i);
		currentMom = currentProp.getMomentum();
		currentProp.setMomentum(currentMom.subs(extMom == 0));
		result.setExtPropagAtIndex(currentProp, i);
	}

// clear ext. momentum from int. props
	for (int i = 0; i < _inputDiag.getIntPropags().size(); i++) {
		currentProp = result.getIntPropags().at(i);
		currentMom = currentProp.getMomentum();
		currentProp.setMomentum(currentMom.subs(extMom == 0));
		result.setIntPropagAtIndex(currentProp, i);
	}

// now delete ext. momentum from diagram attribute vector<ginac::symbol> momenta
	result.setExtMomenta( { });

	return result;
}

//===========================================================================

//===========================================================================
// Code for getting part proportional to tau - TODO
std::vector<Diagram> addVertCorrespToTauDeriv(Diagram _parentDiag) {
// Here we deal with question of finding divergent part of two point functions
// proportional to tau. The function is analogical to finding part propto ext. freq.
//
// This can all be done by performing derivative w.r.t. tau
// On graph level in (p,t) representation this corresponds to attaching new unit
// vertex on some internal edge carring tau. (all such insertions must
// be found). This function then implements this concept, and returns vector
// of diagrams, each of which has insertion on one of the internal legs
// that carries tau). Such results can be further looked at to find
// time versions and etc.

// NOTE: at this point it works for 2 point functions only and IPMomRouting
// is assumed (_parentDiagram already has split mP props).

	std::vector<Diagram> result;
	Diagram created = _parentDiag;
	std::stringstream createdDiagName;
	Vertex createdVert;
	Propagator split1;
	Propagator split2;
	std::vector<Propagator> newProps = { };
	std::vector<Vertex> newVerts = { };
	Propagator currentProp;

// where there is ext. momentum there is ext. frequency which we are after
// i.e. the new vertex can only be inserted to internal leg carrying ext. momentum

// now we only programmed IPMomRouting for 2 point functions, so there is only one
// external momentum p0 (in IPMomRouting)
// TODO - more general case
//GiNaC::symbol extMom = created.getExtMomenta().at(0);

// this will be name of newly created vertex (only one needs to be created - one frequency derivative)
	std::string name = "f";

// just some tests to catch possible problems
	if (_parentDiag.getExtPropags().size() > 2) {
		std::cout
				<< "Right now I cannot do diagrams of 3 point and more point functions. -addVertCorrespToFrequencyInt()"
				<< std::endl;
		return result;
	}

// back to code
// go through all internal propagators
	for (int i = 0; i < _parentDiag.getIntPropags().size(); i++) {
		created = _parentDiag;

		currentProp = created.getIntPropags().at(i);

		// i want to only divide paropagators with tau
		if (currentProp.getPropType().compare("MM") == 0) {
			continue;
		}

		if (currentProp.getPropType().compare("pP") == 0) {
			// if we are dividing pP propagator

			// create new vertex
			createdVert = Vertex(name, "pP");

			// split the propagator into two
			split1 = Propagator(currentProp.getEndVert(), createdVert,
					currentProp.getMomentum(), "pP");
			split2 = Propagator(createdVert, currentProp.getStartVert(),
					currentProp.getMomentum(), "pP");

			// one of the newly created props directly replaces old propagator in diagram
			created.setIntPropagAtIndex(split2, i);
			newProps.push_back(split1);
			newVerts.push_back(createdVert);

			// the rest are appended
			appendVectors(newProps, created.getIntPropags());
			appendVectors(newVerts, created.getVertices());
			created.setIntPropags(newProps);
			created.setVertices(newVerts);

			// give returend diagram some different name to differentiate it from parent
			createdDiagName.str("");
			createdDiagName << "tau_" << result.size() << "_"
					<< _parentDiag.getName();
			created.setName(createdDiagName.str());

			// if it is proportional to tau it will not be to ext. momentum/frequency
			// therefore ext. mom may be set to zero
			result.push_back(setExtMomToZero(created));

			// reset variables for next iteration
			newProps.clear();
			newVerts.clear();
		}

		if (currentProp.getPropType().compare("MP") == 0) {
			// if we are dividing MP propagator

			// create new vertex
			createdVert = Vertex(name, "pP");

			// split the propagator into two
			// TODO .....
			split1 = Propagator(currentProp.getEndVert(), createdVert,
					currentProp.getMomentum(), "MP");
			split2 = Propagator(createdVert, currentProp.getStartVert(),
					currentProp.getMomentum(), "pP");

			// one of the newly created props directly replaces old propagator in diagram
			created.setIntPropagAtIndex(split2, i);
			newProps.push_back(split1);
			newVerts.push_back(createdVert);

			// the rest are appended
			appendVectors(newProps, created.getIntPropags());
			appendVectors(newVerts, created.getVertices());
			created.setIntPropags(newProps);
			created.setVertices(newVerts);

			// give returend diagram some different name to differentiate it from parent
			createdDiagName.str("");
			createdDiagName << "tau_" << result.size() << "_"
					<< _parentDiag.getName();
			created.setName(createdDiagName.str());

			// if it is proportional to tau it will not be to ext. momentum/frequency
			// therefore ext. mom may be set to zero
			result.push_back(setExtMomToZero(created));

			// reset variables for next iteration
			newProps.clear();
			newVerts.clear();
		}

	}

	return result;
}

//===========================================================================

//===========================================================================
// Code for getting part proportional to ext. momentum (p^2)

std::vector<Diagram> addVertsCorrespToPDeriv(Diagram &_parentDiag) {
	// Function combines effects of the functions addVertType1, addVertsType2
	// and addVertsType3_crossTerms
	std::vector<Diagram> result;
	std::vector<Diagram> res1 = addVertType1(_parentDiag);
	std::vector<Diagram> res2 = addVertsType2(_parentDiag);
	std::vector<Diagram> res3 = addVertsType3_crossTerms(_parentDiag);

	for (int i = 0; i < res1.size(); i++) {
		result.push_back(res1.at(i));
	}

	for (int i = 0; i < res2.size(); i++) {
		result.push_back(res2.at(i));
	}

	for (int i = 0; i < res3.size(); i++) {
		result.push_back(res3.at(i));
	}

	return result;
}

std::vector<Diagram> addVertType1(Diagram &_parentDiag) {
	// Vert type 1 when differentiating w.r.t p^2 is the one analogous to
	// tau derivative (second derivative acts on the factor of momenta in numerator
	// brought about by first derivative)
	std::vector<Diagram> result;
	Diagram created = _parentDiag;
	std::stringstream createdDiagName;
	Vertex createdVert;
	Propagator split1;
	Propagator split2;
	std::vector<Propagator> newProps = { };
	std::vector<Vertex> newVerts = { };
	Propagator currentProp;
	GiNaC::symbol D0 = _parentDiag.getD_0();

	// i.e. the new vertex can only be inserted to internal leg carrying ext. momentum

	// now we only programmed IPMomRouting for 2 point functions, so there is only one
	// external momentum p0 (in IPMomRouting)

	//GiNaC::symbol extMom = created.getExtMomenta().at(0);

	// this will be name of newly created vertex (only one needs to be created)
	std::string name = "f";

	// just some tests to catch possible problems
	if (_parentDiag.getExtPropags().size() > 2) {
		std::cout
				<< "Right now I cannot do diagrams of 3 point and more point functions. -addVertCorrespToFrequencyInt()"
				<< std::endl;
		return result;
	}

	// naturally the function can be used only when diagram has non-zero ext. mom. (and ext. freq.)
	for (int i = 0; i < _parentDiag.getExtPropags().size(); i++) {
		if (_parentDiag.getExtPropags().at(i).getMomentum() == 0) {
			std::cout
					<< "Illegal actions. Diagram has ext. momenta set to zero (freq's also)."
							"Therefore derivative w.r.t. ext. mom gives zero. -addVertCorrsepToMomDeriv()"
					<< std::endl;
			return result;
		}
	}

	// since only 2 point functions are allowed and IPMomRouting must be used
	if (_parentDiag.getExtMomenta().size() != 1) {
		std::cout
				<< "Illegal actions. You probably didn't use IPMomRouting which is required"
						"-addVertCorrsepToMomDeriv()" << std::endl;
		return result;
	}

	GiNaC::symbol extMom = created.getExtMomenta().at(0);
	// back to code
	// go through all internal propagators
	for (int i = 0; i < _parentDiag.getIntPropags().size(); i++) {
		created = _parentDiag;

		currentProp = created.getIntPropags().at(i);

		// i want to only divide paropagators with p
		if (!currentProp.getMomentum().has(extMom)) {
			continue;
		}

		// create new vertex
		createdVert = Vertex(name, "pP");

		// split the propagator into two
		split1 = Propagator(currentProp.getEndVert(), createdVert,
				currentProp.getMomentum(), "pP");
		split2 = Propagator(createdVert, currentProp.getStartVert(),
				currentProp.getMomentum(), "pP");

		// one of the newly created props directly replaces old propagator in diagram
		created.setIntPropagAtIndex(split2, i);
		newProps.push_back(split1);
		newVerts.push_back(createdVert);

		// the rest are appended
		appendVectors(newProps, created.getIntPropags());
		appendVectors(newVerts, created.getVertices());
		created.setIntPropags(newProps);
		created.setVertices(newVerts);

		// give returend diagram some different name to differentiate it from parent
		createdDiagName.str("");
		createdDiagName << "type1_" << result.size() << "_"
				<< _parentDiag.getName();
		created.setName(createdDiagName.str());

		// you have to multiply this by factor -D_0 - because of derivatives
		GiNaC::ex newFactor = created.getFactorsFromVertsAndProps();
		//GiNaC::ex newFactor = created.getFactorsFromVertsAndProps() * (-D0);
		created.setFactorsFromVertsAndProps(newFactor);

		//GiNaC::ex createdNumerator = (-D0*D0);
		GiNaC::ex createdNumerator = (-D0);
		created.setCreatedNumeratorByDerivative(createdNumerator);
		// if it is proportional to tau it will not be to ext. momentum/frequency
		// therefore ext. mom may be set to zero
		result.push_back(setExtMomToZero(created));

		// reset variables for next iteration
		newProps.clear();
		newVerts.clear();

	}

	return result;
}

std::vector<Diagram> addVertsType2(Diagram &_parentDiag) {
	// Vert type 2 when differentiating w.r.t p^2 is the one where
	// second derivative acts again on the same propagator as first
	// this brings about factor of momentum of the prop. squared into numerator
	// and 2 unit vertices
	std::vector<Diagram> result;
	Diagram created = _parentDiag;
	std::stringstream createdDiagName;
	Vertex createdVert1;
	Vertex createdVert2;
	Propagator split1;
	Propagator split2;
	Propagator split3;
	std::vector<Propagator> newProps = { };
	std::vector<Vertex> newVerts = { };
	Propagator currentProp;
	GiNaC::symbol D0 = _parentDiag.getD_0();
	GiNaC::symbol d = _parentDiag.getD();

	// i.e. the new vertex can only be inserted to internal leg carrying ext. momentum

	// now we only programmed IPMomRouting for 2 point functions, so there is only one
	// external momentum p0 (in IPMomRouting)

	//GiNaC::symbol extMom = created.getExtMomenta().at(0);

	// this will be name of newly created vertex
	std::string name1 = "f";
	std::string name2 = "g";

	// just some tests to catch possible problems
	if (_parentDiag.getExtPropags().size() > 2) {
		std::cout
				<< "Right now I cannot do diagrams of 3 point and more point functions. -addVertCorrespToFrequencyInt()"
				<< std::endl;
		return result;
	}

	// naturally the function can be used only when diagram has non-zero ext. mom. (and ext. freq.)
	for (int i = 0; i < _parentDiag.getExtPropags().size(); i++) {
		if (_parentDiag.getExtPropags().at(i).getMomentum() == 0) {
			std::cout
					<< "Illegal actions. Diagram has ext. momenta set to zero (freq's also)."
							"Therefore derivative w.r.t. ext. mom gives zero. -addVertCorrsepToMomDeriv()"
					<< std::endl;
			return result;
		}
	}

	// since only 2 point functions are allowed and IPMomRouting must be used
	if (_parentDiag.getExtMomenta().size() != 1) {
		std::cout
				<< "Illegal actions. You probably didn't use IPMomRouting which is required"
						"-addVertCorrsepToMomDeriv()" << std::endl;
		return result;
	}

	GiNaC::symbol extMom = created.getExtMomenta().at(0);
	// back to code
	// go through all internal propagators
	for (int i = 0; i < _parentDiag.getIntPropags().size(); i++) {
		created = _parentDiag;

		currentProp = created.getIntPropags().at(i);

		// i want to only divide paropagators with p
		if (!currentProp.getMomentum().has(extMom)) {
			continue;
		}

		// create new vertex
		createdVert1 = Vertex(name1, "pP");
		createdVert2 = Vertex(name2, "pP");
		// split the propagator into two
		split1 = Propagator(currentProp.getEndVert(), createdVert1,
				currentProp.getMomentum(), "pP");
		split2 = Propagator(createdVert1, createdVert2,
				currentProp.getMomentum(), "pP");
		split3 = Propagator(createdVert2, currentProp.getStartVert(),
				currentProp.getMomentum(), "pP");

		// one of the newly created props directly replaces old propagator in diagram
		created.setIntPropagAtIndex(split2, i);
		newProps.push_back(split1);
		newProps.push_back(split3);
		newVerts.push_back(createdVert1);
		newVerts.push_back(createdVert2);

		// the rest are appended
		appendVectors(newProps, created.getIntPropags());
		appendVectors(newVerts, created.getVertices());
		created.setIntPropags(newProps);
		created.setVertices(newVerts);

		// give returend diagram some different name to differentiate it from parent
		createdDiagName.str("");
		createdDiagName << "type2_" << result.size() << "_"
				<< _parentDiag.getName();
		created.setName(createdDiagName.str());

		// you have to multiply this by factor (+4/d)*D_0^2 - because of derivatives
		//GiNaC::ex newFactor = created.getFactorsFromVertsAndProps();
		//GiNaC::ex newFactor = created.getFactorsFromVertsAndProps()
		//		* (4 * GiNaC::pow(D0, 2) / d);
		//created.setFactorsFromVertsAndProps(newFactor);

		// there is factor of momentum squared of the given propagator that appears
		// in numerator and factor of (4 * GiNaC::pow(D0, 2) / d)
		// the former is not written - it is taken care of later thanks to fact
		// that it is type 2
		// the latter is explicitly written here and later added to I in FinalIntegral

		//GiNaC::ex createdNumerator = currentProp.getMomentum().subs(
		//		extMom == 0);
		//createdNumerator = GiNaC::pow(createdNumerator, 2);
		GiNaC::ex createdNumerator = (4 * GiNaC::pow(D0, 2) / d);
		created.setCreatedNumeratorByDerivative(createdNumerator);
		std::vector<int> indices = {};
		//created.print();
		//split2.print();
		//std::cout<<findIndexOfPropInsideVector(split2, created.getIntPropags());
		indices.push_back(findIndexOfPropInsideVector(split2, created.getIntPropags()));
		//std::cout<<indices.at(0);
		created.setIndicesOfPropsWhichContrToNumerator(indices);
		//created.print();
		// if it is proportional to tau it will not be to ext. momentum/frequency
		// therefore ext. mom may be set to zero
		result.push_back(setExtMomToZero(created));

		// reset variables for next iteration
		newProps.clear();
		newVerts.clear();

	}
	return result;
}

std::vector<Diagram> addVertsType3_crossTerms(Diagram &_parentDiag) {
	// Cross terms - one derivative acts on one propagator carrying ext. frequency
	// second on another propag. carrying ext. frequency - unit vertex added to both
	// proppagators and factor with dot product of resp. momenta appears in numerator
	std::vector<Diagram> result;
	Diagram created = _parentDiag;
	std::stringstream createdDiagName;
	Vertex createdVert1;
	Vertex createdVert2;
	Propagator split1;
	Propagator split2;
	Propagator split3;
	Propagator split4;
	std::vector<Propagator> newProps = { };
	std::vector<Vertex> newVerts = { };
	Propagator currentProp1;
	Propagator currentProp2;
	GiNaC::symbol D0 = _parentDiag.getD_0();
	GiNaC::symbol d = _parentDiag.getD();

	// i.e. the new vertex can only be inserted to internal leg carrying ext. momentum

	// now we only programmed IPMomRouting for 2 point functions, so there is only one
	// external momentum p0 (in IPMomRouting)

	//GiNaC::symbol extMom = created.getExtMomenta().at(0);

	// this will be name of newly created vertex
	std::string name1 = "f";
	std::string name2 = "g";

	// just some tests to catch possible problems
	if (_parentDiag.getExtPropags().size() > 2) {
		std::cout
				<< "Right now I cannot do diagrams of 3 point and more point functions. -addVertCorrespToFrequencyInt()"
				<< std::endl;
		return result;
	}

	// naturally the function can be used only when diagram has non-zero ext. mom. (and ext. freq.)
	for (int i = 0; i < _parentDiag.getExtPropags().size(); i++) {
		if (_parentDiag.getExtPropags().at(i).getMomentum() == 0) {
			std::cout
					<< "Illegal actions. Diagram has ext. momenta set to zero (freq's also)."
							"Therefore derivative w.r.t. ext. frequency gives zero. -addVertCorrsepToFrequencyDeriv()"
					<< std::endl;
			return result;
		}
	}

	// since only 2 point functions are allowed and IPMomRouting must be used
	if (_parentDiag.getExtMomenta().size() != 1) {
		std::cout
				<< "Illegal actions. You probably didn't use IPMomRouting which is required"
						"-addVertCorrsepToFrequencyDeriv()" << std::endl;
		return result;
	}

	GiNaC::symbol extMom = created.getExtMomenta().at(0);
	// back to code
	// go through all pairs of internal propagators
	for (int i = 0; i < _parentDiag.getIntPropags().size(); i++) {
		for (int j = 0; j < _parentDiag.getIntPropags().size(); j++) {
			created = _parentDiag;

			// if i==j then we have vertsType2 - already done
			if (i == j) {
				continue;
			}

			currentProp1 = created.getIntPropags().at(i);
			currentProp2 = created.getIntPropags().at(j);

			// i want to only divide paropagators with p
			if (!currentProp1.getMomentum().has(extMom)) {
				continue;
			}
			if (!currentProp2.getMomentum().has(extMom)) {
				continue;
			}

			// create new vertex
			createdVert1 = Vertex(name1, "pP");
			createdVert2 = Vertex(name2, "pP");
			// split the propagator into two
			split1 = Propagator(currentProp1.getEndVert(), createdVert1,
					currentProp1.getMomentum(), "pP");
			split2 = Propagator(createdVert1, currentProp1.getStartVert(),
					currentProp1.getMomentum(), "pP");
			split3 = Propagator(currentProp2.getEndVert(), createdVert2,
					currentProp2.getMomentum(), "pP");
			split4 = Propagator(createdVert2, currentProp2.getStartVert(),
					currentProp2.getMomentum(), "pP");

			// one of the newly created props directly replaces old propagator in diagram
			created.setIntPropagAtIndex(split2, i);
			created.setIntPropagAtIndex(split4, j);
			newProps.push_back(split1);
			newProps.push_back(split3);
			newVerts.push_back(createdVert1);
			newVerts.push_back(createdVert2);

			// the rest are appended
			appendVectors(newProps, created.getIntPropags());
			appendVectors(newVerts, created.getVertices());
			created.setIntPropags(newProps);
			created.setVertices(newVerts);

			// give returend diagram some different name to differentiate it from parent
			createdDiagName.str("");
			createdDiagName << "type3_" << result.size() << "_"
					<< _parentDiag.getName();
			created.setName(createdDiagName.str());



			GiNaC::ex newFactor = created.getFactorsFromVertsAndProps();
			//GiNaC::ex newFactor = created.getFactorsFromVertsAndProps()
			//		* (2 * d * GiNaC::pow(D0, 2));
			created.setFactorsFromVertsAndProps(newFactor);

			// there is factor of dot product of momenta of the given propagators that are diff'ed
			// in numerator and factor of (+2*d)*D_0^2
			// the former is not written - it is taken care of later thanks to fact
			// that it is type 3
			// the latter is explicitly written here and later added to I in FinalIntegral

			GiNaC::ex createdNumerator = (2 * d * GiNaC::pow(D0, 2));
			created.setCreatedNumeratorByDerivative(createdNumerator);
			std::vector<int> indices = {};
			indices.push_back(findIndexOfPropInsideVector(split2, created.getIntPropags()));
			indices.push_back(findIndexOfPropInsideVector(split4, created.getIntPropags()));
			created.setIndicesOfPropsWhichContrToNumerator(indices);
			// if it is proportional to tau it will not be to ext. momentum/frequency
			// therefore ext. mom may be set to zero
			result.push_back(setExtMomToZero(created));

			// reset variables for next iteration
			newProps.clear();
			newVerts.clear();

		}
	}
	return result;
}

//===========================================================================
