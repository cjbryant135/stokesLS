/*
 *  units.h
 *  levelset
 *
 *  Created by David Chopp on 3/31/10.
 *  Copyright 2010 Northwestern University. All rights reserved.
 *
 */

#ifndef UNITS_H
#define UNITS_H

#include <string>
#include <vector>
#include <map>
#include <list>

namespace levelset {
	
class Units {
private:
	
	static const char primitivechar = '!';
	std::map<std::string,std::string> unittable;
	typedef struct {
#if 0
		std::vector<std::string> numerator;
		std::vector<std::string> denominator;
#else
		std::list<std::string> numerator;
		std::list<std::string> denominator;		
#endif
		double factor;
	} unittype;
	std::vector<std::pair<std::string,std::string> > prefixtable;
	
public:
	
	Units(void) {ReadUnits();}
	double Convert(const double amt, const std::string have, const std::string want);
	
private:
	
	void ReadUnits(void);
	bool AddUnit(unittype* theunit, std::string toadd, bool flip);
	std::string LookUpUnit(const std::string unit);
	int ReduceProduct(unittype* theunit, bool flip);
	int ReduceUnit(unittype* theunit);
	void CancelUnits(unittype* unit);
	int CompleteReduce(unittype* unit);
	
	double Reduce(const std::string expression);
};

}
#endif
