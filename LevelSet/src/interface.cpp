#include "interface.h"
#include "numtrait.h"
#ifdef USE_AVS
#include "avsfield.h"
#endif
#include <sstream>
#include <vector>

namespace levelset {

    int Interface::LoopCount(void) 
    {
        char* mark = new char[length];
        int i;
        for (i=0; i<length; ++i) mark[i] = 0x00;
        int count = 0;
        int start;
        for (start=0; start<length && mark[start]; ++start) ;
        while (start < length) {
            int here = start;
            for (here=start; here!=-1 && !mark[here]; here=(*data)[here].nabor[SRight])
                mark[here] = 0x01;
            if (here == -1) {
                int last;
                for (here=start, last=here; here!=-1; here=(*data)[here].nabor[SLeft]) {
                    mark[here] = 0x01;
                    last = here;
                }
                mark[last] = 0x02;
            } else {
                mark[here] = 0x02;
            }
            ++count;
            for (start=0; start<length && mark[start]; ++start) ;
        }
        if (loopstart) delete[] loopstart;
        loopstart = new int[count];
        for (start=0, i=0; start<length; ++start)
            if (mark[start] == 0x02) loopstart[i++] = start;
        delete[] mark;
        return count;
    }

    int Interface::LoopLength(const int loop) const
    {
        int count;
        int here;
        for (count=0, here=loopstart[loop]; 
             here!=-1 && (count==0 || here!=loopstart[loop]); 
             ++count, here=(*data)[here].nabor[SRight]) ;
        return count;
    }

    void Interface::InterpToEnds(const int si, const int sl, const int sh)
    {
        int nbr;
        double len[2];
   
        for (int i=0; i<length; ++i) {
            len[0] = (*data)[i].Length();
            if ((nbr=(*data)[i].nabor[SLeft]) != -1) {
                len[1] = (*data)[nbr].Length();
                (*data)[i][sl] = (len[0]*(*data)[nbr][si]+len[1]*(*data)[i][si])
                    /(len[0]+len[1]);
            }
            else {
                (*data)[i][sl] = (*data)[i][si];
            }
            if ((nbr=(*data)[i].nabor[SRight]) != -1) {
                len[1] = (*data)[nbr].Length();
                (*data)[i][sh] = (len[0]*(*data)[nbr][si]+len[1]*(*data)[i][si])
                    /(len[0]+len[1]);
            }
            else {
                (*data)[i][sh] = (*data)[i][si];
            }
        }
    }

    void Interface::GetEndpoints(const int loop, int& n, double*& x, double*& y) const
    {
        n = LoopLength(loop);
        char c = Closed(loop);
        if (!c) ++n;
        x = new double[n];
        y = new double[n];
        
        int here;
        int i;
        for (i=0, here=loopstart[loop]; here!=-1 && (i==0 || here!=loopstart[loop]); 
             ++i, here=(*data)[here].nabor[SRight]) {
            x[i] = (*data)[here].X1();
            y[i] = (*data)[here].Y1();
            if (!c && (*data)[here].nabor[SRight]==-1) {
                x[i+1] = (*data)[here].X2();
                y[i+1] = (*data)[here].Y2();
            }
        }
    }

    void Interface::GetEndpoints(int& n, double*& x, double*& y) const
    {
        x = new double[length+1];
        y = new double[length+1];
        
        n=0;
        for (int i=0; i<length; ++i) {
            x[n] = (*data)[i].X1();
            y[n++] = (*data)[i].Y1();
            if ((*data)[i].nabor[SRight]==-1) {
                x[n] = (*data)[i].X2();
                y[n] = (*data)[i].Y2();
            }
        }
    }

    void Interface::GetEndpoints(std::vector<Path>& p) const
    {
        p.clear();
        bool *mark = new bool[length];
        for (int i=0; i<length; ++i)
            mark[i] = false;
        int n = 0;
        
        while (n < length) {
            // start a path from first unmarked segment
            int j;
            for (j=0; j<length && mark[j]; ++j) ;
            Path q;
            q.Clear();
            int k = j;
            q.x.insert(q.x.begin(),(*data)[k].X1());
            q.y.insert(q.y.begin(),(*data)[k].Y1());
            q.x.push_back((*data)[k].X2());
            q.y.push_back((*data)[k].Y2());
            mark[k] = true; ++n;
            k = (*data)[j].nabor[SLeft];
            while (k != j && k >= 0) {
                q.x.insert(q.x.begin(),(*data)[k].X1());
                q.y.insert(q.y.begin(),(*data)[k].Y1());
                mark[k] = true; ++n;
                k = (*data)[k].nabor[SLeft];
            }
            if (k == j) {
                q.x.insert(q.x.begin(),(*data)[k].X1());
                q.y.insert(q.y.begin(),(*data)[k].Y1());
            } else {
                k = (*data)[j].nabor[SRight];
                while (k >= 0) {
                    q.x.push_back((*data)[k].X2());
                    q.y.push_back((*data)[k].Y2());
                    mark[k] = true; ++n;
                    k = (*data)[j].nabor[SRight];
                }
            }
            p.push_back(q);
        }
    }

    void Interface::GetEndpoints(const int loop, int& n, DataVec<double>& x, DataVec<double>& y) const
    {
        n = LoopLength(loop);
        char c = Closed(loop);
        if (!c) ++n;
        x.Resize(n);
        y.Resize(n);
        
        int here;
        int i;
        for (i=0, here=loopstart[loop]; here!=-1 && (i==0 || here!=loopstart[loop]); 
             ++i, here=(*data)[here].nabor[SRight]) {
            x[i] = (*data)[here].X1();
            y[i] = (*data)[here].Y1();
            if (!c && (*data)[here].nabor[SRight]==-1) {
                x[i+1] = (*data)[here].X2();
                y[i+1] = (*data)[here].Y2();
            }
        }
    }

    char Interface::Closed(const int loop) const
    {
        return (*data)[loopstart[loop]].nabor[SRight] != -1
            && (*data)[loopstart[loop]].nabor[SLeft] != -1;
    }
   
    char Interface::Verify(const int sindex) const
    {
        char answer = true;
   
        for (int ii=0; ii<length; ++ii)
            answer = answer && (*data)[ii].Verify(sindex);
        return answer;
    }

    char Interface::Verify(void) const
    {
        char* mark = new char[length];
        int i;
        char ans = true;
        for (i=0; i<length; ++i) mark[i] = 0x00;
        int start = 0;
        while (start < length) {
            int here;
            for (here=start; here!=-1 && !mark[here]; here=(*data)[here].nabor[SRight])
                mark[here] = 0x01;
            if (here == -1) {
                int last;
                for (here=(*data)[start].nabor[SLeft], last=start; here!=-1 && !mark[here]; 
                     here=(*data)[here].nabor[SLeft]) {
                    mark[here] = 0x01;
                    last = here;
                }
                if (here != -1) {
                    std::cout << "Error in interface with start = " << start << std::endl;
                    ans = false;
                }
                mark[last] = 0x02;
            } else {
                if (here != start) {
                    std::cout << "Error in interface with start = " << start << std::endl;
                    ans = false;
                }
                mark[here] = 0x02;
            }
            for (start=0; start<length && mark[start]; ++start) ;
        }
        delete[] mark;
        return ans;
    }

#ifndef NO_GRAPHICS
    void Interface::Plot(PlotWindow2D& port, const int k,
                         const char showvel, const double scale) const
    {
        for (int ii=0; ii<length; ++ii)
            (*data)[ii].Plot(port, k, showvel, scale);
    }

    void Interface::LabelPlot(PlotWindow2D& port) const
    {
        for (int ii=0; ii<length; ++ii) {
            (*data)[ii].Plot(port, 0);
            std::ostringstream s;
            s << ii << std::ends;
            port.Text(((*data)[ii].X1()+(*data)[ii].X2())/2,
                      ((*data)[ii].Y1()+(*data)[ii].Y2())/2, s.str().c_str());
        }
    }
#endif

#ifdef USE_AVS
    ostream& AVS_Out(ostream& s, Interface* m, const char* comment)
    {
        int loopcnt = m->LoopCount();
        int dim = m->Length() + loopcnt - 1;
        NumericTrait<double> t;
        char datatype[20] = "double";
        AVSWriteHead(s, comment, 1, &dim, 2, 1, datatype, IrregularField,
                     NULL, NULL, NULL, NULL, NULL, NULL);
        int ii;
        for (ii=0; ii<dim; ++ii)
            AVSWriteBinary(s, double(1.0));
        for (ii=0; ii<loopcnt; ++ii) {
            int looplen = m->LoopLength(ii);
            Segment seg;
            int j;
            for (j=0, seg=m->LoopStart(ii); j<looplen; ++j, seg=(*m)[seg.nabor[SRight]])
                AVSWriteBinary(s, seg.x[0]);
            AVSWriteBinary(s, seg.x[0]);
        }
        for (ii=0; ii<loopcnt; ++ii) {
            int looplen = m->LoopLength(ii);
            Segment seg;
            int j;
            for (j=0, seg=m->LoopStart(ii); j<looplen; ++j, seg=(*m)[seg.nabor[SRight]])
                AVSWriteBinary(s, seg.y[0]);
            AVSWriteBinary(s, seg.y[0]);
        }
        return s;
    }
#endif

    int Interface::GetLoopIndices(const int loop, int*& ind) const
    {
        int len = LoopLength(loop);
        ind = new int[len];
        int i;
        for (i=0, ind[0] = loopstart[loop]; i<len-1; ++i) 
            ind[i+1] = (*data)[ind[i]].nabor[SRight];
        return len;
    }

}
