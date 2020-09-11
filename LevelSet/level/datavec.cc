#include "datavec.h"

namespace levelset {
	
template <class T>
DataVec<T>::DataVec(const int n, const T* array)
{
   data = new DataBlock<T>(n, array);
   data->AddRef();
   length = n;
}

template <class T>
DataVec<T>::DataVec(const DataVec<T>& vec)
{
   data = vec.data;
   length = vec.length;
   data->AddRef();
}

template <class T>
DataVec<T> DataVec<T>::Copy(void) const
{
   DataVec<T> newvec;
   newvec.data = data->Copy();
   newvec.length = length;
   newvec.data->AddRef();
   return newvec;
}

template <class T>
DataVec<T> DataVec<T>::Ref(void) const
{
   DataVec<T> newvec;
   newvec.data = data;
   newvec.length = length;
   newvec.data->AddRef();
   return newvec;
}

template <class T>
void DataVec<T>::Resize(const int n)
{
   if (data) {
      if (n > length) 
         (*data).Resize(n);
   }
   else {
      if (n > 0) {
         data = new DataBlock<T>(n);
         data->AddRef();
      }
   }
   length = n;
}

template <class T>
DataVec<T>& DataVec<T>::operator=(const DataVec<T>& vec)
{
   if (data) FreeRef(data);
   data = vec.data;
   data->AddRef();
   length = vec.length;
   return *this;
}

template <class T>
DataVec<T>& DataVec<T>::operator+=(const DataVec<T>& vec)
{
   int len = length;
   Resize(length+vec.Length());
   for (int i=0; i<vec.Length(); ++i) (*data)[len+i] = vec[i];
   return *this;
}

template <class T>
char DataVec<T>::operator==(const DataVec<T>& v1)
{
   if (length == v1.length) {
      for (int i=0; i<length; ++i)
         if ((*this)[i] != v1[i]) return 0x00;
      return 0x01;
   } else {
      return 0x00;
   }
}

template <class T>
std::ostream& operator<<(std::ostream& s, DataVec<T>& v)
{
   s << v.length << std::endl;
   for (int i = 0; i < v.length; ++i)
      s << v[i] << std::endl;
   return s;
}

template <class T>
std::istream& operator>>(std::istream& s, DataVec<T>& v)
{
   int i;
   s >> i;
   T* data = new T[i];
   for (int j = 0; j < i; ++j)
      s >> data[j];
   DataVec<T> temp(i, data);
   v = temp;
   return s;
}
	

}



