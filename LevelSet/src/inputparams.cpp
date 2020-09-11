#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include "tinyxml.h"
#include "inputparams.h"
#include "units.h"

namespace levelset {

    int Value::index_ = 0;

    template <class T>
    std::string toString(const T n)
    {
        std::ostringstream s;
        s << n;
        return s.str();
    }

    std::string toString(const bool b)
    {
        if (b) 
            return std::string("true");
        else
            return std::string("false");
    }

    static Units UConverter;

    void Value::SetUnits(const std::string u)
    {
        if (units != "") {
            double newval = UConverter.Convert(double(*this), units, u);
            s = toString(newval);
        }
        units = u;
    }

    void IPElement::AddAttribute(const std::string label, const std::string value)
    {
        attrib.insert(std::make_pair(label, value));
    }

    Value IPElement::GetAttribute(const std::string label) 
    {
        Value v = attrib[label];
        return v;
    }

    bool IPElement::IsAttribute(const std::string label) const
    {
        return attrib.find(label) != attrib.end();
    }

    Value IPElement::GetValue(const std::string unit) 
    {
        if (unit == "") {
            return GetAttribute("!value");
        } else {
            Value val = GetAttribute("!value");
            val.SetUnits(unit);
            return val;
        }
    }

    std::multimap<std::string,IPElement> IPElement::GetChildren(std::string label)
    {
        if (label == std::string("!all")) {
            return child;
        } else {
            std::multimap<std::string,IPElement> subset;
            std::pair<std::multimap<std::string,IPElement>::iterator,std::multimap<std::string,IPElement>::iterator> rng;
            rng = child.equal_range(label);
            subset.insert(rng.first,rng.second);
            return subset;
        }
    }

    void IPElement::AddChild(const IPElement& elem)
    {
        child.insert(std::make_pair(elem.name,elem));
    }

    IPElement IPElement::GetChild(const std::string label) const
    {
        std::multimap<std::string,IPElement>::const_iterator i = child.find(label);
        if (i == child.end()) {
            std::cerr << "Error: unable to locate parameter " << label << '\n';
            exit(1);
        }
        return i->second;
    }

    std::ostream& operator<<(std::ostream& s, const IPElement& elem)
    {
        s << elem.name << '\n';
        for (std::map<std::string,Value>::const_iterator i=elem.attrib.begin(); i!=elem.attrib.end(); ++i)
            s << i->first << " = " << std::string(i->second) << '\n';
        for (std::multimap<std::string,IPElement>::const_iterator i=elem.child.begin(); i!=elem.child.end(); ++i)
            s << i->second << '\n';
        return s;
    }

    InputParams::InputParams(char* filename) : fname(filename), maketemplate(false), bname(filename) 
    {
        std::string fname = filename;
        isXML = !fname.substr(fname.length()-3).compare("xml");
        ReadFile(filename);
    }

    InputParams::InputParams(int nargs, char* argv[])
        : fname(argv[1]), maketemplate(false), bname(argv[1])
    {
        std::string fname = argv[1];
        isXML = !fname.substr(fname.length()-3).compare("xml");
        ReadFile(argv[1]);
        for (int i=2; i<nargs; ++i)
            AddParam(argv[i]);
    }

    void InputParams::ExtractNode(TiXmlNode* node, IPElement& elem) 
    {
        switch(node->Type()) {
        case TiXmlNode::TINYXML_DOCUMENT:
        {
            std::cout << "Document:\n" << node->Value() << "\n";
            for (TiXmlNode* cnode = node->FirstChild(); cnode != NULL; cnode = cnode->NextSibling()) 
                ExtractNode(cnode, root);
        }
        break;
        case TiXmlNode::TINYXML_ELEMENT: 
        {
            //std::cout << "IPElement:\n" << node->Value() << "\n";
            TiXmlElement* enode = node->ToElement();
            std::string newlabel;
            if (std::string(node->Value()) == std::string("param")) {
                enode->QueryStringAttribute("name",&newlabel);
            } else {
                newlabel = node->Value();
            }
            IPElement newelem(newlabel);
            for (TiXmlAttribute* pAtt = enode->FirstAttribute(); pAtt != NULL; pAtt = pAtt->Next()) {
                newelem.AddAttribute(pAtt->Name(),pAtt->ValueStr());
                //std::cout << "  " << pAtt->Name() << " = " << pAtt->ValueStr() << '\n';
            }
            for (TiXmlNode* cnode = enode->FirstChild(); cnode != NULL; cnode = cnode->NextSibling()) 
                ExtractNode(cnode, newelem);
            elem.AddChild(newelem);
        }
        break;
        case TiXmlNode::TINYXML_TEXT:
            elem.AddAttribute("!value",node->Parent()->ToElement()->GetText());
            //std::cout << "    Value = " << node->Parent()->ToElement()->GetText() << '\n';
            break;
        default:
            break;
        }
    }

    void InputParams::ReadFile(const char* fname)
    {
        if (isXML) {
            xmldoc.LoadFile(fname);
            ExtractNode(&xmldoc, root);

            // bvm: commented this out 7/23/10 because the 'version 1.0' code is meant to
            // be a simple list of parameters, but in reality the version # refers to the
            // version of xml so this isn't quite right. For now, just ignore the xml file
            // type that is just a list of numbers

//              TiXmlDeclaration* declaration = xmldoc.FirstChild()->ToDeclaration();
//              std::string foo = declaration->Version();
//              if (declaration->Version() == std::string("1.0")) {
//                      TiXmlHandle docHandle(&xmldoc);
//                      TiXmlHandle h1 = docHandle.FirstChildElement("Parameters");
//                      TiXmlElement* h2 = h1.FirstChildElement("param").ToElement();
//                      while (h2 != NULL) {
//                              std::string key;
//                              h2->QueryStringAttribute("name",&key);
//                              Value val;
//                              val = h2->GetText();
//                              ipmap::iterator iter = values.find(key);
//                              if (iter == values.end()) {
//                                      values.insert(std::make_pair(key,val));
//                              } else {
//                                      iter->second = val;
//                              }
//                              h2 = h2->NextSiblingElement();
//                      }
//              } else {
//                      ExtractNode(&xmldoc, root);
//              //      std::cout << "\n\nDatabase:\n" << root << '\n';
//              }
        } else {
            std::ifstream file(fname);
            char buffer[256];
            std::string key;
            Value val;
                
            if (file.good()) {
                while (!file.eof()) {
                    file.getline(buffer, 256);
                    if (buffer[0] != '#') {
                        key = buffer;
                        int n = key.find(":");
                        if (n != std::string::npos) {
                            val = key.substr(n+1);
                            key.erase(n);
                            ipmap::iterator iter = values.find(key);
                            if (iter == values.end()) {
                                values.insert(std::make_pair(key,val));
                            }
                            else {
                                iter->second = val;
                            }
                        }
                    }
                }
            } else {
                std::cout << "File " << fname << " does not exist, create it? (y/n) ";
                char c;
                std::cin >> c;
                maketemplate = c == 'y';
            }   
        }
    }

    void InputParams::AddParam(const char* buffer)
    {
        std::string key = buffer;
        int n = key.find(":");
        if (n == std::string::npos)
            n = key.find("=");
        if (n != std::string::npos) {
            Value val;
            val = key.substr(n+1);
            key.erase(n);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                updates.insert(std::make_pair(key,val));
            }
            else {
                iter->second = val;
            }
        }
    }

    Value InputParams::GetParamValue(const char* name)
    {
        Value *val;
        
        if (isXML) {
            std::string key(name);
            std::string::size_type index = key.find_first_of('/');
            IPElement elem = root;
            while (index != std::string::npos) {
                elem = elem.GetChild(key.substr(0,index));
                key = key.substr(index+1);
                index = key.find_first_of('/');
            }
            if (elem.IsAttribute(key)) {
                return Value(elem.GetAttribute(key));
            }
            elem = elem.GetChild(key);
            return Value(elem.GetAttribute("!value"));
        } else {
            if (!maketemplate) {
                std::string key(name);
                ipmap::iterator iter = updates.find(key);
                if (iter == updates.end()) {
                    ipmap::iterator iter2 = values.find(key);
                    if (iter2 != values.end()) {
                        val = &(iter2->second);
                    }
                    else {
                        std::cerr << "Error: parameter " << key << " not found.\n";
                        exit(1);
                    }
                }
                else {
                    ipmap::iterator iter2 = values.find(key);
                    if (iter2 != values.end()) {
                        iter2->second = iter->second;
                        val = &(iter->second);
                    }
                    else {
                        std::cerr << "Error: parameter " << key << " not found.\n";
                        exit(1);
                    }
                }
            }
            else {
                std::string key(name);
                ipmap::iterator iter = updates.find(key);
                if (iter == updates.end()) {
                    std::cout << name << ':';
                    std::string ans;
                    std::cin >> ans;
                    Value v(ans);
                    values.insert(std::make_pair(key,v));
                    ipmap::iterator iter2 = values.find(key);
                    val = &(iter2->second);
                }
                else {
                    values.insert(*iter);
                    ipmap::iterator iter2 = values.find(key);
                    val = &(iter2->second);
                }
            }
        }
        return *val;
    }

    double InputParams::GetDoubleParamWithUnits(const char* name, const char* unit)
    {   
        if (isXML) {
            std::string key(name);
            std::string::size_type index = key.find_first_of('/');
            IPElement elem = root;
            while (index != std::string::npos) {
                elem = elem.GetChild(key.substr(0,index));
                key = key.substr(index+1);
                index = key.find_first_of('/');
            }
            Value val;
            if (elem.IsAttribute(key)) {
                val = Value(elem.GetAttribute(key));
            } else {
                elem = elem.GetChild(key);
                val = elem.GetAttribute("!value");
            }
            if (elem.IsAttribute("unit")) 
                val.SetUnits(elem.GetAttribute("unit"));
            val.SetUnits(unit);
            return double(val);
        } else {
            std::cerr << "GetParamWithUnits is only available for XML files.\n";
            exit(1);
        }
    }

    Value InputParams::GetParamValue(const char* name, const std::string defalt, 
                                     const char* desc)
    {
        Value *val;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    val = &(iter2->second);
                }
                else {
                    Value v(defalt);
                    if (desc != NULL) v.SetDescription(desc);
                    values.insert(std::make_pair(key,v));
                    iter2 = values.find(key);
                    val = &(iter2->second);
                }
            }
            else {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    iter2->second = iter->second;
                    val = &(iter->second);
                }
                else {
                    std::cerr << "Error: parameter " << key << " not found.\n";
                    exit(1);
                }
            }  
        }
        else {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                if (desc != NULL) std::cout << "# " << desc << '\n';
                std::cout << name << ':';
                std::string ans;
                std::cin >> ans;
                std::string key(name);
                if (ans.size() == 0) ans = defalt;
                Value v(ans);
                if (desc != NULL) v.SetDescription(desc);
                values.insert(std::make_pair(key,v));
                ipmap::iterator iter2 = values.find(key);
                val = &(iter2->second);
            }
            else {
                std::cout << name << ':' << std::string(iter->second) << '\n';
                values.insert(*iter);
                val = &(iter->second);
            }
        }
        return *val;
    }

    std::vector<IPElement> InputParams::GetParamList(const char* name)
    {
        if (isXML) {
            std::string key(name);
            std::string::size_type index = key.find_first_of('/');
            IPElement elem = root;
            while (index != std::string::npos) {
                elem = elem.GetChild(key.substr(0,index));
                key = key.substr(index+1);
                index = key.find_first_of('/');
            }
            std::multimap<std::string,IPElement> elist = elem.GetChildren(key);
            std::vector<IPElement> elems;
            for (std::multimap<std::string,IPElement>::const_iterator i=elist.begin(); i!=elist.end(); ++i)
                elems.push_back(i->second);
            return elems;
        } else {
            std::cerr << "Error: GetParamList valid only for XML input files.\n";
            exit(1);
        }
    }

    IPElement InputParams::GetElement(const char* name)
    {
        if (isXML) {
            std::string key(name);
            std::string::size_type index = key.find_first_of('/');
            IPElement elem = root;
            while (index != std::string::npos) {
                elem = elem.GetChild(key.substr(0,index));
                key = key.substr(index+1);
                index = key.find_first_of('/');
            }
            return elem.GetChild(key);
        } else {
            std::cerr << "Error: GetParamList valid only for XML input files.\n";
            exit(1);
        }
    }

    template <class T>
    T InputParams::GetParam(const char* name, const T defalt, const char* desc)
    {
        Value val = GetParamValue(name, toString(defalt), desc);
        T ans = T(val);
        //std::cout << "Parameter: " << name << " = " << ans << '\n';
        return ans;
    }

    int InputParams::GetIntParam(const char* name)
    {
        Value val = GetParamValue(name);
        int ans = int(val);
        //std::cout << "Parameter: " << name << " = " << ans << '\n';
        return ans;
    }

    int InputParams::GetIntParam(const char* name, const int defalt, 
                                 const char* desc)
    {
        return GetParam(name, defalt, desc);
    }

    double InputParams::GetDoubleParam(const char* name)
    {
        Value val = GetParamValue(name);
        double ans = double(val);
        //std::cout << "Parameter: " << name << " = " << ans << '\n';
        return ans;
    }

    double InputParams::GetDoubleParam(const char* name, const double defalt, 
                                       const char* desc)
    {
        return GetParam(name, defalt, desc);
    }

    void InputParams::GetCharParam(const char* name, char* ans)
    {
        Value val = GetParamValue(name);
        strcpy(ans,std::string(val).c_str());
        //std::cout << "Parameter: " << name << " = " << ans << '\n';
    }

    void InputParams::GetCharParam(const char* name, char* ans, 
                                   const char* defalt, const char* desc)
    {
        Value val = GetParamValue(name, std::string(defalt), desc);
        strcpy(ans,std::string(val).c_str());
        //std::cout << "Parameter: " << name << " = " << ans << '\n';
    }

    std::string InputParams::GetStringParam(const char* name, const char* defalt, const char* desc)
    {
        std::string df;
        if (defalt == NULL) {
            df = std::string("");
        } else {
            df = defalt;
        }
        Value val = GetParamValue(name /*, df, desc*/);  
        std::string ans = std::string(val);
        //std::cout << "Parameter: " << name << " = " << ans << '\n';
        return ans;
    }

    bool InputParams::GetBooleanParam(const char* name)
    {
        Value val = GetParamValue(name);
        bool ans = val.isTrue();
        //std::cout << "Parameter: " << name << " = " << (ans ? "true" : "false") << '\n';
        return ans;
    }

    bool InputParams::GetBooleanParam(const char* name, const bool defalt, 
                                      const char* desc)
    {
        Value val = GetParamValue(name, toString(defalt), desc);
        bool ans = val.isTrue();
        //std::cout << "Parameter: " << name << " = " << (ans ? "true" : "false") << '\n';
        return ans;
    }


#if 0
    int InputParams::GetIntParam(const char* name)
    {
        int ans;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    ans = int(iter2->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }
                else {
                    std::cerr << "Error: parameter " << key << " not found.\n";
                    exit(1);
                }
            }
            else {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    iter2->second = iter->second;
                    ans = int(iter->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }
                else {
                    std::cerr << "Error: parameter " << key << " not found.\n";
                    exit(1);
                }
            }
        }
        else {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                std::cout << name << ':';
                std::cin >> ans;
                std::string key(name);
                Value val(toString(ans),values.size());
                values.insert(std::make_pair(key,val));
            }
            else {
                std::cout << name << ':' << int(iter->second) << '\n';
                values.insert(*iter);
                ans = int(iter->second);
            }
        }
        
        return ans;
    }

    int InputParams::GetIntParam(const char* name, const int defalt, const char* desc)
    {
        int ans;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    ans = int(iter2->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }
                else {
                    ans = defalt;
                    Value val(toString(defalt),values.size());
                    val.SetDescription(desc);
                    values.insert(std::make_pair(key,val));
                }
            }
            else {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    iter2->second = iter->second;
                    ans = int(iter->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }
                else {
                    std::cerr << "Error: parameter " << key << " not found.\n";
                    exit(1);
                }
            }  
        }
        else {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                std::cout << "# " << desc << '\n';
                std::cout << name << ':';
                std::cin >> ans;
                std::string key(name);
                Value val(toString(ans),values.size());
                val.SetDescription(desc);
                values.insert(std::make_pair(key,val));
            }
            else {
                std::cout << name << ':' << int(iter->second) << '\n';
                values.insert(*iter);
                ans = int(iter->second);
            }
        }
        
        return ans;
    }

    double InputParams::GetDoubleParam(const char* name)
    {
        double ans;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    ans = double(iter2->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }
                else {
                    std::cerr << "Error: parameter " << key << " not found.\n";
                    exit(1);
                }
            }
            else {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    iter2->second = iter->second;
                    ans = double(iter->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }
                else {
                    std::cerr << "Error: parameter " << key << " not found.\n";
                    exit(1);
                }
            }
        }
        else {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                std::cout << name << ':';
                std::cin >> ans;
                std::string key(name);
                Value val(toString(ans),values.size());
                values.insert(std::make_pair(key,val));
            }
            else {
                ans = double(iter->second);
                std::cout << name << ':' << ans << '\n';
                values.insert(*iter);
            }
        }
        
        return ans;
    }

    double InputParams::GetDoubleParam(const char* name, const double defalt, 
                                       const char* desc)
    {
        double ans;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = updates.find(key);
            if (iter == updates.end()) {
                ipmap::iterator iter2 = values.find(key);
                if (iter2 != values.end()) {
                    ans = double(iter2->second);
                    std::cout << "Parameter: " << name << " = " << ans << '\n';
                }    
                else {
                    ans = defalt;
                    Value val(toString(defalt),values.size());
                    val.SetDescription(desc);
                    values.insert(std::make_pair(key,val));
                }
            }
            else {
                        
            }
        }
        else {
            std::cout << "# " << desc << '\n';
            std::cout << name << ':';
            std::cin >> ans;
            std::string key(name);
            Value val(toString(ans),values.size());
            val.SetDescription(desc);
            values.insert(std::make_pair(key,val));
        }
        
        return ans;
    }

    void InputParams::GetCharParam(const char* name, char* ans)
    {
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = values.find(key);
            if (iter != values.end()) {
                strcpy(ans,std::string(iter->second).c_str());
                std::cout << "Parameter: " << name << " = " << ans << '\n';
            }
            else {
                std::cerr << "Error: parameter " << key << " not found.\n";
                exit(1);
            }
        }
        else {
            std::cout << name << ':';
            std::cin >> ans;
            std::string key(name);
            Value val(values.size());
            val = std::string(ans);
            values.insert(std::make_pair(key,val));
        }
    }

    void InputParams::GetCharParam(const char* name, char* ans, const char* defalt, const char* desc)
    {
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = values.find(key);
            if (iter != values.end()) {
                strcpy(ans, std::string(iter->second).c_str());
                std::cout << "Parameter: " << name << " = " << ans << '\n';
            }
            else {
                strcpy(ans,defalt);
                std::string foo(defalt);
                Value val(foo, values.size());
                val.SetDescription(desc);
                values.insert(std::make_pair(key,val));
            }
        }
        else {
            std::cout << "# " << desc << '\n';
            std::cout << name << ':';
            std::cin >> ans;
            std::string key(name);
            Value val(toString(ans), values.size());
            val.SetDescription(desc);
            values.insert(std::make_pair(key,val));
        }
    }

    bool InputParams::GetBooleanParam(const char* name)
    {
        bool ans;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = values.find(key);
            if (iter != values.end()) {
                ans = iter->second.isTrue();
                std::cout << "Parameter: " << name << " = " 
                          << (ans ? "true" : "false") << '\n';
            }
            else {
                std::cerr << "Error: parameter " << key << " not found.\n";
                exit(1);
            }
        }
        else {
            char buffer[256];
            std::cout << name << ':';
            std::cin >> buffer;
            std::string key(name);
            std::string foo(buffer);
            Value val(foo, values.size());
            values.insert(std::make_pair(key,val));
            ans = val.isTrue();
        }
        
        return ans;
    }

    bool InputParams::GetBooleanParam(const char* name, const bool defalt, const char* desc)
    {
        bool ans;
        
        if (!maketemplate) {
            std::string key(name);
            ipmap::iterator iter = values.find(key);
            if (iter != values.end()) {
                ans = iter->second.isTrue();
                std::cout << "Parameter: " << name << " = " 
                          << (ans ? "true" : "false") << '\n';
            }
            else {
                ans = defalt;
                Value val(defalt ? std::string("true") : std::string("false"), 
                          values.size());
                val.SetDescription(desc);
                values.insert(std::make_pair(key,val));
            }
        }
        else {
            std::cout << "# " << desc << '\n';
            std::cout << name << ':';
            char buffer[256];
            std::cin >> buffer;
            std::string key(name);
            std::string foo(buffer);
            Value val(foo, values.size());
            val.SetDescription(desc);
            values.insert(std::make_pair(key,val));
            ans = val.isTrue();
        }
        
        return ans;
    }
#endif

    void InputParams::CloseTemplate(void)
    {
        if (maketemplate) {
            if (!isXML) {
                std::ofstream s(fname.c_str());
                std::map<int,std::pair<std::string,Value> > values2;
                for (ipmap::iterator i=values.begin(); i!=values.end(); ++i) 
                    values2.insert(std::make_pair(i->second.Index(),
                                                  std::make_pair(i->first,i->second)));
                std::map<int,std::pair<std::string,Value> >::iterator j;
                for (j=values2.begin(); j!=values2.end(); ++j) {
                    if (j->second.second.Description().length() > 0) 
                        s << "# " << j->second.second.Description() << '\n';
                    s << j->second.first << ": " << std::string(j->second.second) << '\n';
                }
                s.close();
            } else {
                WriteXML(fname.c_str());
            }
        }
    }

    void InputParams::WriteXML(const char* name)
    {
        std::string oname;
        if (name == NULL) {
            oname = fname;
            if (!fname.substr(fname.length()-3).compare("dat")) {
                oname = fname.substr(0,fname.length()-3)+"xml";
            } else {
                oname = fname+".xml";
            }
        } else {
            oname = name;
        }
        if (!isXML) {
            xmldoc.SetValue(oname);
            TiXmlElement* root = xmldoc.RootElement();
            std::map<int,std::pair<std::string,Value> >::iterator j;
            for (ipmap::iterator j=values.begin(); j!=values.end(); ++j) {
                TiXmlElement foo("param");
                foo.SetAttribute("name",j->first);
                TiXmlText data(j->second);
                foo.InsertEndChild(data);
                root->InsertEndChild(foo);
            }
                
        }
        xmldoc.SaveFile(oname);
    }

}


