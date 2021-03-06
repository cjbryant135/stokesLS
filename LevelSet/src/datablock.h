#ifndef __DATABLOCK_H__
#define __DATABLOCK_H__

#include "reference.h"
#include "defs.h"

namespace levelset {
/// 
    template <class T>
        class DataBlock : public Reference {
    private:

        unsigned int length;
        T*  data;
//   int refs;

    public:

        /** @name Constructors and Destructors */
        //@{
        /// Default constructor leaves data array null and length 0
    DataBlock(void) : Reference(), length(0), data(NULL) {;}
        /** Construct a vector of length n and allocate space 
         *      with all values initially 0
         */
        DataBlock(const int n, const T* array = NULL);
        /// Destructor disposes of allocated memory
        ~DataBlock(void);
        //@}
        
        /** @name Methods */
        //@{
        /// Get the length of the vector
        inline unsigned int Length(void) const { return length; }
        /// Make a unique copy of the data vector
        DataBlock<T>* Copy(void);
        /// Truncate or allocate unitialized entries
        void Resize(const unsigned int n);
        /// Get the actual data pointer
        T* Ptr(void) {return data;}
        //@}

        /** @name Operators */
        //@{
        /// Access elements of the vector, no bounds checking
        inline T operator[](int i) const {return data[i];}
        ///
        inline T& operator[](int i) {return data[i];}
        //@}

    };

}

#ifndef NO_INCLUDE_CC_FILE
#include "datablock.cc"
#endif

#endif // __DATABLOCK_H__
