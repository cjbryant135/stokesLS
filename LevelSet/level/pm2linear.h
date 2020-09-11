//
//  pm2linear.hpp
//  LevelSet
//
//  Created by David Chopp on 7/23/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#ifndef pm2linear_hpp
#define pm2linear_hpp

#include "pm2boundary.h"

namespace levelset {
    
    class PM2_LinearBdry : public PM2_Boundary
    {
    public:
        
        PM2_LinearBdry(const int wxlo = 1, const int wxhi = 1,
                       const int wylo = 1, const int wyhi = 1)
        : PM2_Boundary(wxlo, wxhi, wylo, wyhi) {}
        void Apply(const int k);
        void Apply(const int i, const int j, const int k);
    };
    
}

#endif /* pm2linear_hpp */
