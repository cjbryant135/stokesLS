/*
 *  dat2xml.cpp
 *  levelset
 *
 *  Created by David Chopp on 3/25/10.
 *  Copyright 2010 Northwestern University. All rights reserved.
 *
 */


#include "inputparams.h"
#include "units.h"
#include <iostream>

using namespace levelset;

static Units UConverter;

int main(int argc, char* argv[])
{
    InputParams oldfile(argv[1]);
        
    // read control
        
    double d1 = oldfile.GetDoubleParamWithUnits("control/startOfSimulation", "hours");
    double d2 = oldfile.GetDoubleParamWithUnits("control/endOfSimulation", "min");
    double d3 = oldfile.GetDoubleParamWithUnits("control/timeStep", "sec");
    double d4 = oldfile.GetDoubleParamWithUnits("control/outputPeriod", "years");
    std::cout << "Start at time " << d1 << " hours.\n";
    std::cout << "End at time " << d2 << " minutes.\n";
    std::cout << "Each time step is " << d3 << " seconds long.\n";
    std::cout << "Print every interval of length " << d4 << " years.\n";
        
    // read computationDomain
        
    int nx = oldfile.GetIntParam("computationDomain/grid/nX");  // note: I added this funcionality after we talked
    int ny = oldfile.GetIntParam("computationDomain/grid/nY");
    double resx = oldfile.GetDoubleParamWithUnits("computationDomain/grid/resolutionX", "m");
    double resy = oldfile.GetDoubleParamWithUnits("computationDomain/grid/resolutionY", "yards");
    std::cout << "The domain has dimensions " << nx << "x" << ny << '\n';
    std::cout << "The x resolution is " << resx << " meters per grid cell\n";
    std::cout << "The y resolution is " << resy << " yards per grid cell\n";
        
    IPElement e1 = oldfile.GetElement("computationDomain/initfunc");
    std::cout << "Initial domain is the union of:\n";
    std::multimap<std::string,IPElement> iflist = e1.GetChildren();
    for (std::multimap<std::string,IPElement>::iterator i=iflist.begin(); i!=iflist.end(); ++i) {
        std::cout << i->first << ": ";
        double cx = i->second.GetAttribute("centerX");
        double cy = i->second.GetAttribute("centerY");
        double r = i->second.GetAttribute("radius");
        std::cout << "(" << UConverter.Convert(cx,i->second.GetAttribute("unit"),"inches") << ',' 
                  << UConverter.Convert(cy,i->second.GetAttribute("unit"),"inches") << ") and radius "
                  << UConverter.Convert(r,i->second.GetAttribute("unit"),"inches") << " in inches.\n";
    }
        
    std::vector<IPElement> bclist = oldfile.GetParamList("computationDomain/boundaryCondition");
    for (std::vector<IPElement>::iterator i = bclist.begin(); i!=bclist.end(); ++i) {
        std::cout << "For " << std::string(i->GetAttribute("location")) 
                  << " the bc's are " << std::string(i->GetAttribute("class"));
        if (std::string(i->GetAttribute("class")) == "Periodic") 
            std::cout << " with " << std::string(i->GetAttribute("with"));
        std::cout << '\n';
    }
        
    // read biofilm 
    // am going to skip a few things here onlyl because they are similar to the above examples.
        
    double bl = oldfile.GetDoubleParamWithUnits("biofilm/boundaryLayer", "angstroms");
    std::cout << "The boundary layer is " << bl << " angstroms.\n";
    double sa = oldfile.GetDoubleParamWithUnits("biofilm/specificArea", "1/m");
    std::cout << "The specific area is " << sa << " per meter.\n";
    bool detach = oldfile.GetBooleanParam("biofilm/maintainVolumeFractionSum");
    std::string bm = oldfile.GetStringParam("biofilm/maintainVolumeFractionSum/biomass");
    std::cout << "maintainVolumeFractionSum is " << (detach ? "on" : "off") 
              << " and reference biomass is " << bm << '\n';
        
    //std::cout << foo << '\n';
//      std::vector<Element> alist = oldfile.GetParamList("biomass");
//      for (std::vector<Element>::iterator i=alist.begin(); i!=alist.end(); ++i)
//              std::cout << std::string(i->GetAttribute("name")) << '\n';
//      oldfile.WriteXML("deleteme.xml");
    return 0;
}
