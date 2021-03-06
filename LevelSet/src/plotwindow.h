#ifndef H_PLOTWINDOW_CLASS
#define H_PLOTWINDOW_CLASS

#include <X11/Xos.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

#ifndef true
#define true 0x01
#define false 0x00
#endif

#include <fstream>

namespace levelset {
        
    struct xpcb {
        Display *dpy;
        int     screen;
        Window  win;
        GC      gc;
        int     width;
        int     height;
        double  xmax, ymax;
        double  xmin, ymin;
    };
//#ifndef WWidth
//#define WWidth 800
//#endif // WWidth
//#ifndef WHeight
//#define WHeight 800
//#endif // WHeight
#ifndef FullScale
#define FullScale 10000.0
#endif // FullScale
#define xflushWaiting 0
#define NumberOfColors 80
#define DOTSIZE 6

    class PlotWindow {

        struct xpcb    *itsWindow;
        struct xpcb    *tempWindow;
        XFontStruct    *StandardFont;
        double         theScale;
        int            trans[2];
        std::ofstream*      PSFile;
        int            shade;
        int            filecnt;
        char           basename[50];
        char           pause_at_exit;
        int                             WWidth;
        int                             WHeight;

    protected:
   
        // Plotting Primitives
   
        void   PlotSegment(int x0, int y0, int x1, int y1, 
                           unsigned short cshade = 0);
        void   PlotText(int x0, int y0, const char* text);
        void   PlotTriangle(int x0, int y0, int x1, int y1, int x2, int y2,
                            unsigned short cshade);
        void   PlotPolygon(int num, int* xp, int* yp, unsigned short cshade);
        void   PlotDot(int x0, int y0, int color);
        void   PlotOval(int x0, int y0, int w, int h);
        void     PlotFilledOval(int x0, int y0, int w, int h);
        void     Clear(int x0, int y0, int w, int h);
   
        //   Utilities
   
        inline int   Round(double x)   {return int(x+0.5);};
   
    public:

        //PlotWindow(const int w = 800, const int h = 800) 
        //              : pause_at_exit(true), WWidth(w), WHeight(h)  {;}
    PlotWindow(bool yes_window = true, bool use_pause = true, const int w = 800, const int h = 800)
        : pause_at_exit(use_pause), WWidth(w), WHeight(h) {Initialize(yes_window);}
        ~PlotWindow(void);
        void Initialize(bool yes_window);
   
        // Window Hiding/Viewing command
   
        void   ShowPlotWindow(bool yes);
   
        //   General graphics commands
   
        void   SetColor(int);
        void   SetWidth(int);
        void   Clear(const bool pause = false);

        enum WinPos {LL=0, LR, UL, UR};
        void   PlotText(WinPos wp, const char* s);

        // Postscript commands
   
        bool StartNewFile(const char* name);
        bool StartNewFile(void);
        void EndFile(void);
        void Comment(char* s);
        void Comment(std::istream& s);
   
        // Image commands
        
        XImage* GetImage(void);
        void SaveToTiff(char* fname);
        void SaveToTiff(const char* bname, const int num, const int digits = 6);
    };

}
#endif

