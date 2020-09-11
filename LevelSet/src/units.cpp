/*
 * units.c   Copyright (c) 1993 by Adrian Mariano (adrian@cam.cornell.edu)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * Disclaimer:  This software is provided by the author "as is".  The author
 * shall not be liable for any damages caused in any way by this software.
 *
 * I would appreciate (though I do not require) receiving a copy of any
 * improvements you might make to this program.
 */

/* C++ version by David Chopp */

#ifndef lint
static const char rcsid[] =
    "$FreeBSD: src/usr.bin/units/units.c,v 1.10 2002/07/28 16:23:28 dwmalone Exp $";
#endif /* not lint */

#include <ctype.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

//#define       _PATH_UNITSLIB  "/usr/share/misc/units.lib"

#define VERSION "1.0"

#ifndef UNITSFILE
#define UNITSFILE _PATH_UNITSLIB
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <ctype.h>
#include "units.h"

namespace levelset {
        
    std::list<std::string> tokenize(const std::string line, const std::string delim)
    {
        std::list<std::string> result;
        std::string::size_type first = line.find_first_not_of(delim);
        while (first != std::string::npos) {
            std::string::size_type last = line.find_first_of(delim, first);
            result.push_back(line.substr(first, last-first));
            first = line.find_first_not_of(delim, last);
        }
        return result;
    }

    void Units::ReadUnits(void)
    {
        std::string buffer;
        char* ufile;
        ufile = getenv("UNITSFILE");
        std::ifstream unitfile;
        if (ufile) {
            unitfile.open(ufile);
        } else { 
            unitfile.open(UNITSFILE);
        }
        if (!unitfile.is_open()) {
            std::cerr << "Unable to locate units file.  Set the environment variable UNITSFILE to the correct value.\n";
            exit(1);
        }
        while (getline(unitfile, buffer)) {
            if (buffer.size() > 0 && buffer[0] != '/') {
                std::list<std::string> toks = tokenize(buffer, " \n\t");
                if (toks.size() > 0) {
                    if (toks.front()[toks.front().size()-1] == '-')  
                        prefixtable.push_back(std::make_pair(toks.front().substr(0,toks.front().size()-1),toks.back()));
                    else {
                        std::list<std::string>::iterator i = toks.begin();
                        ++i;
                        std::string foo = buffer.substr(buffer.find(*i));
                        unittable.insert(std::make_pair(toks.front(),buffer.substr(buffer.find(*i,toks.front().size()))));
                    }
                }
            }
        }
#if 0   
        for (std::map<std::string,std::string>::iterator i=unittable.begin(); i!=unittable.end(); ++i) 
            std::cout << i->first << " --> " << i->second << '\n';
#endif
    }

#if 0
    bool Units::AddUnit(unittype* theunit, std::string toadd, bool flip)
    {
        std::string scratch = toadd;
        
        // remove minus signs that aren't part of an exponential notation
        std::string::size_type index = scratch.find_first_of("-");
        while (index != std::string::npos) {
            if (index == 0 || tolower(scratch[index-1]) != 'e' 
                || !(isdigit(scratch[index+1]) || scratch[index+1]=='.') ) 
                scratch[index] = ' ';
            index = scratch.find_first_of("-",index+1);
        }
        index = scratch.find_first_of("/");
        std::vector<std::string> toks;
        if (index == std::string::npos) 
            toks = tokenize(scratch, " *\t\n");
        else { 
            toks = tokenize(scratch.substr(0,index), " *\t\n/");
            toks.push_back(std::string("/"));
            std::vector<std::string> toks2 = tokenize(scratch.substr(index+1), " *\t\n/");
            toks.insert(toks.end(), toks2.begin(), toks2.end());
        }

        int doingtop = 1;
        for (std::vector<std::string>::iterator i=toks.begin(); i!=toks.end() && doingtop >=0; ++i) {
            if (*i == std::string("/")) 
                --doingtop;
            else {
                if (isnumber(*i->c_str())) {
                    if (doingtop ^ flip) 
                        theunit->factor *= atof(i->c_str());
                    else
                        theunit->factor /= atof(i->c_str());
                } else {
                    int repeat = 1;
                    if (isdigit(*(i->rbegin()))) {
                        repeat = *(i->rbegin()) - '0';
                        *i = i->substr(0,i->size()-1);
                    }
                    for (; repeat; --repeat) {
                        if (doingtop ^ flip)
                            theunit->numerator.push_back(*i);
                        else
                            theunit->denominator.push_back(*i);
                    }
                }
            }
        }
        return 0;
    }

    std::string Units::LookUpUnit(const std::string unit)
    {
        std::map<std::string,std::string>::iterator i = unittable.find(unit);
        if (i != unittable.end()) {
            return i->second;
        } else {
            std::string tunit = unit;
            if (*tunit.rbegin() == '^') {
                tunit.erase(tunit.size()-1);
                i = unittable.find(tunit);
                if (i != unittable.end()) {
                    return tunit;
                }
            }
            if (*tunit.rbegin() == 's') {
                tunit.erase(tunit.size()-1);
                i = unittable.find(tunit);
                if (i != unittable.end()) {
                    return tunit;
                } else if (*tunit.rbegin() == 'e') {
                    tunit.erase(tunit.size()-1);
                    i = unittable.find(tunit);
                    if (i != unittable.end()) {
                        return tunit;
                    }
                }
            }
            for (std::map<std::string,std::string>::iterator j = prefixtable.begin(); j != prefixtable.end(); ++j) {
                if (j->first == unit.substr(0,j->first.size())) {
                    if ((unit.size() <= j->first.size()) || LookUpUnit(unit.substr(j->first.size())) != "") {
                        return j->first+std::string(" ")+unit.substr(j->first.size());
                    }
                }
            }
        }
        return std::string("");
    }

    int Units::ReduceProduct(unittype* theunit, bool flip)
    {
        int didsomething = 2;
        std::vector<std::string> product = flip ? theunit->denominator : theunit->numerator;
        for (std::vector<std::string>::iterator i = product.begin(); i != product.end(); ++i) {
            bool done = false;
            while (i->size() > 0 && !done) {
                std::string toadd = LookUpUnit(*i);
                if (toadd == "") {
                    std::cerr << "Unknown unit " << *i << ". terminating.\n";
                    exit(1);
                } else {
                    if (toadd.find(primitivechar) != std::string::npos) 
                        done = true;
                    else {
                        didsomething = 1;
                        *i = std::string("");
                        AddUnit(theunit, toadd, flip);
                    }
                }
            }
        }
        return didsomething;
    }

    int Units::ReduceUnit(unittype* theunit)
    {
        int ret = 1;
        while (ret & 1) {
            ret = ReduceProduct(theunit, 0) | ReduceProduct(theunit, 1);
            if (ret & 4)
                return 1;
        }
        return 0;
    }

    void Units::CancelUnits(unittype* unit)
    {
        for (std::vector<std::string>::iterator i=unit->numerator.begin(); i!=unit->numerator.end(); ++i) {
            std::vector<std::string>::iterator j;
            for (j=unit->denominator.begin(); j!=unit->denominator.end() && *i != *j; ++j) ;
            if (j != unit->denominator.end()) {
                *i = std::string("");
                *j = std::string("");
            }
        }
    }

    int Units::CompleteReduce(unittype* unit)
    {
        if (ReduceUnit(unit))
            return 1;
        CancelUnits(unit);
        return 0;
    }
#endif

    double Units::Reduce(const std::string expression)
    {
        std::string::size_type index;
        double factor = 1.;
        if ((index = expression.find_first_of('/')) != std::string::npos) {
            factor = Reduce(expression.substr(0,index))/Reduce(expression.substr(index+1));
        } else {
            std::list<std::string> toklist;
            toklist = tokenize(expression, " \t\n");
            std::list<std::string>::iterator i = toklist.begin();
            std::list<std::string>::iterator j;
            while (i != toklist.end()) {
#if 0
                // debug by listing tok list
                std::cout << "toklist:\n";
                for (j=toklist.begin(); j!=toklist.end(); ++j) 
                    std::cout << *j << '\n';
                std::cout << '\n';
#endif
                if ((*i)[0] == '!') {
                    ++i;
                } else {
                    std::map<std::string,std::string>::iterator k = unittable.find(*i);
                    if (k != unittable.end()) {
                        *i = k->second;
                    } else {
                        if ((index = i->find_first_of(" -\t\n")) != std::string::npos) {
                            j = i;
                            toklist.insert(i,i->substr(0,index));
                            --j;
                            toklist.insert(i,i->substr(index+1));
                            toklist.erase(i);
                            i = j;
                        } else {
                            if (i->find_first_of("/") != std::string::npos) {
                                factor *= Reduce(*i);
                                i = toklist.erase(i);
                            } else {
                                if (isdigit((*i)[0])) {  // don't forget about the prefixes
                                    if (index = i->find_first_of('|') != std::string::npos) {
                                        factor *= atof(i->substr(0,index).c_str());
                                        factor /= atof(i->substr(index+1).c_str());
                                    } else
                                        factor *= atof(i->c_str());
                                    j = i;
                                    ++j;
                                    toklist.erase(i);
                                    i = j;
                                } else {
                                    if (isdigit((*i)[i->size()-1])) {
                                        int repeat = (*i)[i->size()-1] - '0';
                                        j = i;
                                        for (int n=0; n<repeat; ++n) {
                                            toklist.insert(i,i->substr(0,i->size()-1));
                                            if (n==0)
                                                --j;
                                        }
                                        toklist.erase(i);
                                        i = j;
                                    } else {
                                        if ((*i)[i->size()-1] == '^') {
                                            *i = i->substr(0,i->size()-1);
                                        } else {
                                            if ((*i)[i->size()-1] == 's') {
                                                *i = i->substr(0,i->size()-1);
                                            } else {
                                                if ((*i)[i->size()-1] == 'e') {
                                                    *i = i->substr(0,i->size()-1);
                                                } else {
                                                    std::vector<std::pair<std::string,std::string> >::iterator pk;
                                                    bool found = false;
                                                    for (pk=prefixtable.begin(); pk!=prefixtable.end() && !found; ++pk) {
                                                        //                                                                                              std::string foo = pk->first;
                                                        //                                                                                              std::string bar = i->substr(0,pk->first.size());
                                                        found = pk->first == i->substr(0,pk->first.size());
                                                        if (found) {
                                                            if (isdigit((pk->second)[0])) {
                                                                *i = i->substr(pk->first.size());
                                                                factor *= atof(pk->second.c_str());
                                                            } else {
                                                                *i = pk->second + i->substr(pk->first.size());
                                                            }
                                                        }
                                                    }
                                                    if (!found) {
                                                        std::cerr << "Unknown unit: " << *i << '\n';
                                                        exit(1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return factor;
    }
                                
                                
                                
    double Units::Convert(const double amt, const std::string have, const std::string want)
    {
#if 0
        unittype haveunit;
        haveunit.factor = 1.;
        unittype wantunit;
        wantunit.factor = 1.;
        AddUnit(&haveunit, have, 0);
        AddUnit(&wantunit, want, 0);
        CompleteReduce(&haveunit);
        CompleteReduce(&wantunit);
        return amt*haveunit.factor/wantunit.factor;
#else
        return amt*Reduce(have)/Reduce(want);
#endif
    }

}
