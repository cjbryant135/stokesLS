#ifndef __DATABLOCK_CC__
#define __DATABLOCK_CC__

namespace levelset {

template <class T>
DataBlock<T>::DataBlock(const int size, const T* array) :
      Reference(), length(size)
{
   data = new T[size];
   if (array) {
      int i;
      for (i=0; i<size; ++i) data[i] = array[i];
   }
}

template <class T>
DataBlock<T>::~DataBlock(void)
{
   if (data) delete[] data;
}

template <class T>
DataBlock<T>* DataBlock<T>::Copy(void)
{
   if (data) {
      DataBlock<T>* newblock = new DataBlock<T>(length);
      for (int i=0; i<length; ++i)
         newblock->data[i] = data[i];
      return newblock;
   } else {
      DataBlock<T>* newblock = new DataBlock<T>;
      return newblock;
   }
}


template <class T>
void DataBlock<T>::Resize(const unsigned int n)
{
   if (n > length) {
      T* temp = new T[n];
      if (data) {
         unsigned int i;
         for (i=0; i<length; ++i) temp[i] = data[i];
         delete[] data;
      }
      data = temp;
      length = n;
   }
}

}

#endif // __DATABLOCK_CC__
