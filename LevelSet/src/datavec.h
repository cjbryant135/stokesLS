#ifndef __DATAVEC_H__
#define __DATAVEC_H__

#include <iostream>
#include "datablock.h"
#include "defs.h"

namespace levelset {
        
    template <class T>
        class DataVec;
// Global I/O operators
    template <class T>
        std::ostream& operator<<(std::ostream& s, DataVec<T>& v);
    template <class T>
        std::istream& operator>>(std::istream& s, DataVec<T>& v);

    template <class T>
        class DataVec {
    protected:
   
        DataBlock<T>*    data;
        int              length;

    public:

        // Constructors
   
    DataVec(void) : data(NULL), length(0)  {;}
        DataVec(const int n, const T* array = NULL);
        DataVec(const DataVec<T>& vec);
   
        ~DataVec(void) {FreeRef(data);} 
   
        // Methods
   
        DataVec<T> Copy(void) const;
        DataVec<T> Ref(void) const;
        int Length(void) const {return length;}
        void Resize(const int n);
        T* Ptr(void) {return data->Ptr();} 
   
        // operators
   
        DataVec<T>& operator=(const DataVec<T>& vec);
        DataVec<T>& operator+=(const DataVec<T>& vec);
        inline T& operator[](const unsigned int i)  {return (*data)[i];}
        inline T operator[](const unsigned int i) const {return (*data)[i];}
        char operator==(const DataVec<T>& v1);
        char operator!=(const DataVec<T>& v1) {return !(*this==v1);}
        friend std::ostream& operator<<<>(std::ostream& s, DataVec<T>& v);
        friend std::istream& operator>><>(std::istream& s, DataVec<T>& v);
    };


// Global comparison operators
    template <class T>
        char operator==(const DataVec<T>& v1, const DataVec<T>& v2);

    template <class T>
        inline char operator!=(const DataVec<T>& v1, const DataVec<T>& v2) {return !(v1 == v2);}

}

#ifndef NO_INCLUDE_CC_FILES
#include "datavec.cc"
#endif // NO_INCLUDE_CC_FILES

#endif // __DATAVEC_H__

